'''
Created on Jan 21, 2014

@author: Jason Merritt
'''

# Code to run chemostats/continuous-culture devices, using
# LabJacks (via LabJackPython) for system control and Pt Grey
# cameras with the FlyCapture SDK (by executing a C program)
# for imaging, with external PID and image segmentation
# Python scripts.
#
# Ultimately (with image segmentation script) outputs 4 CSV
# files: a "chemostat" file with experiment data and measurements
# taken by LabJack, a "dataset" file with data on cell image
# segmentation, a "clump dataset" file with data on clump image
# segmentation, and a "counts" file with a very rough unfiltered
# estimate of counts from each image according to data from image
# segmentation.
#
# Also generates a subfolder for images, with images saved by
# C program.

import u3 # LabJackPython import
import math
import numpy as np
import time
import traceback
import subprocess
import sys
import datetime
import os

# Import PID class for temperature control
from PID import PID


# Import thermometer and system number calibration file
#calibration_file = open('calibration.sys','r')
#calibration = np.genfromtxt(calibration_file,delimiter=',')
#calibration_file.close()
systemNumber=1#int(calibration[0]+0.1)

# Create a folder with today's date, and generate a data output file (without overwriting)
path = str(datetime.date.today())
datefolder=path
imgfolder=os.path.join(path,"imgs_S"+str(systemNumber))
filename = 'chemostat_S'+str(systemNumber)
filetype = '.csv'
try: 
    os.makedirs(path)
except OSError:
    if not os.path.isdir(path):
        raise

path = os.path.join(path,str(datetime.date.today())+'_' + filename)
if os.path.isfile(path + filetype):
    testNum = 0
    while os.path.isfile(path + '_' + str(testNum) + filetype):
        testNum+=1
    path=path+'_' + str(testNum)

config_file = sys.argv[0]

###########################
# START CONTROL VARIABLES #
###########################
f = open(config_file, "r")
thirtyDegVolt= float(f.readline())#calibration[1] # SHOULD NOT TYPICALLY BE CHANGED: Voltage from calibration file corresponding to 30 degrees Celsius
degVoltScale= float(f.readline())#calibration[2] # SHOULD NOT TYPICALLY BE CHANGED: Thermometer degree-volt scaling near 30 degrees Celsius

a_uM = float(f.readline())
b_uM = float(f.readline())
c_uM = float(f.readline())

# Temperature/pump cycle format: 2-column array controlling the temperature
# and flow cycles of the chemostat. Each row corresponds to one step in the
# cycle. The first element of each row is a temperature in Celsius (for tempCycle)
# or a nominal volume of media (in mL) to pass through each hour (for pumpCycle).
# The second element of each row is the duration of the step in minutes. For example,
# in pumpCycle, a row of [4.0,24*60] corresponds to attempting to pass 4 mL/hour for
# 24 hours. Note that for pumpCycle this is only a target value, and the actual value
# is set later in the code based on this number. A value of nominalPumpFlowRate*60 in
# the first column of pumpCycle corresponds to maximum allowed flow rate (minimized
# later in the code). Inflow pumps cannot be on for more than 52 seconds per minute.
#
# At the end of tempCycle the cycle simply repeats from the beginning, so that
# the standard tempCycle of [[30,24*60]] simply loops at the end of each day, keeping
# a target of 30 degrees Celsius. At the end of a pumpCycle, the control loop
# returns to the row indicated by pumpCycleReturnIndex.
###########################
tempCycle=[[30,24]]##*60]] # SHOULD NOT TYPICALLY BE CHANGED: Controls temperature
LEDCycle = [[100,10]]##*60]]

if len(sys.argv) > 1:
    os.system(sys.argv[1])
else :
    def temp_f(s):
        return(30)
        
    def LED_f(s):
        if s % 10 == 0:
            print(math.floor(s/10)*10)
        return(math.floor(s/10)*10)
        #return(100)
        # if s < 30:
            # return((10/3)*s)
        # elif s < 45:
            # return(100)
        # elif s < 75:
            # return(100-(10/3)*(s-45))
        # else:
            # return(0)
###########################
#  END CONTROL VARIABLES  #
###########################

for i in tempCycle:
    i[1]=i[1]*60 # Temperature cycle control code functions at the level of seconds, not minutes

# Thermometer pin on LabJack
tempPin=0

# OD measurement pins on LabJack
IRPin=19
diodePin=7 # NOTE: Not directly used anywhere, but the real-world value of this pin
           # IS relevant for readOD(): 'AIN7' refers to pin 7 in that context

# PWM control pins on LabJack for Peltier control (MD)
pwmAPin=5
pwmBPin=6
pwmPin=4

# Microscope and algae LED control pins on LabJack
algaePin=15


# The value pwmAPin should have to control temperature correctly
# Should not be changed
positiveDirection = 1

tempAvgRange = 20 # Number of temperature readings to average over each minute
cycleLength = 1.0 # Length of temperature PID cycle: 1 second

# PID control constants
Kp = 600.0
Ki = 8.0
Kd = 11250.0

# These correspond to PID variables; PID wasn't running before T=0,
# so we just initialize both to 0
lastOutput = 0
lastCycle = 0

# Return the thermometer voltage corresponding to a temperature in Celsius
def voltsToDegree(numVolts):
    return (numVolts-thirtyDegVolt)/degVoltScale+30

# Return the temperature in Celsius corresponding to a thermometer voltage
def degreesToVolts(numDegrees):
    return (numDegrees-30)*degVoltScale+thirtyDegVolt

# Initial temperature goal based on tempCycle
goal = degreesToVolts(tempCycle[0][0])

def uMToLED(uM):
    LED_out = int(a_uM*uM**3+b_uM*uM**2+c_uM*uM+d_uM)
    if LED_out > 255:
        LED_out = 255
    elif LED_out < 0:
        LED_out = 0
    return(LED_out)

# Print header rows for data file
f = open(path+filetype,'w')
#f.write('System ' + str(systemNumber) + ', Temp Cycle: ' + str(tempCycle)+'\n')
f.write('Time (m),Photodiode voltage (V), Photodiode stdev (V), Temp. Setpoint (C), Temp. (C), Thermometer voltage (V), Thermometer stdev (V), Output, Output stdev (V),  Temp. Loop Index')#, Kp = ' + str(Kp) + ",Ki = " + str(Ki) + ",Kd = " + str(Kd) + '\n')
f.close()

# Initialize LabJack (uses LabJackPython)
d = u3.U3()
d.getCalibrationData()
# Initialize LabJack pin types and pin values
d.configIO(FIOAnalog = 129)
d.setDOState(algaePin,1)
d.setDOState(IRPin,1)
d.configIO( NumberOfTimersEnabled = 1 )
baseValue = 65535.0

# Return a measurement of OD (in volts) with an estimate of uncertainty
def readOD():
    # Configure the LabJack for streaming data from OD photodiode
    d.streamConfig( NumChannels = 1, PChannels = [ 7 ], NChannels = [ 31 ], Resolution = 0, ScanFrequency = 100, SamplesPerPacket=24 )
    tempOD=-1.0
    stdev=0
    try:
        d.setDOState(algaePin,1)
        d.getFeedback([u3.DAC8(1,255)])
        d.setDOState(algaePin,0)
        
        d.setDOState(IRPin,0) # Turn on IR LED
        time.sleep(.05) # Briefly wait until after photodiode voltage spike
        d.streamStart() # Begin recording data.
        time.sleep(.25)
        d.setDOState(IRPin,1) # Turn off IR LED

        d.setDOState(algaePin,1)
        d.getFeedback([u3.DAC8(1,uMToLED(new_LED))])
        d.setDOState(algaePin,0)
        
        missed = 0
        dataCount = 0
        
        tempOD=-1.0
        stdev=0
        
        # Check the stream data
        for r in d.streamData():
            if r is not None:
                # If more than one dataset was returned, ignore
                # the later ones...
                if dataCount >= 1:
                    break
                # NOTE: The following errors triggering generally indicate a loose
                # USB cable or a broken LabJack.
                #
                # If there's a serious error, but an exception was not called,
                # record OD as -1.
                if r['errors'] != 0:
                    print( "Error: %s ; " % r['errors'])
                    tempOD = -1
                # Otherwise, if at least 15 readings were recorded...
                elif len(r['AIN7']) > 15:
                    # If readings were missed or we have the wrong number of
                    # data packets, still try to record OD, but record with
                    # negative value
                    if r['numPackets'] != d.packetsPerRequest:
                        print( "----- UNDERFLOW : %s : " % r['numPackets'])
                        tempOD = -sum(r['AIN7'][5:15])/10.0
                        stdev = -np.std(r['AIN7'][5:15])
                    elif r['missed'] != 0:
                        missed += r['missed']
                        print( "+++ Missed ", r['missed'])
                        tempOD = -sum(r['AIN7'][5:15])/10.0
                        stdev = -np.std(r['AIN7'][5:15])
                    else:
                        # If there were no problems, record OD based on
                        # measurements 5 through 14
                        tempOD=sum(r['AIN7'][5:15])/10.0
                        stdev = np.std(r['AIN7'][5:15])
                # If fewer than 15 readings were recorded, record OD as -1.
                else:
                    tempOD = -1
            else:
                # If no data was received, record OD as -1.
                print( "No data")
                tempOD = -1
            dataCount += 1
    except:
        print( "".join(i for i in traceback.format_exc()))
    finally:
        # Stop streaming data and return LabJack to normal operation.
        d.streamStop()
    return tempOD, stdev

# This variable stores the most recent OD value
lastOD = readOD()[0]
# Initial temperature value (un-averaged)
temp0 = d.getAIN(tempPin)

# Instantiate PID for temperature control
p = PID(newInput = temp0, newSetpoint = goal, newSampleTime=cycleLength, newKp=Kp, newKi=Ki, newKd=Kd)
# This variable stores the most recent suggested output from the PID
lastOutput = p.compute(temp0)

# See description of updatePWM() below
if lastOutput >= 0:
    d.setDOState(pwmAPin,positiveDirection)
    d.setDOState(pwmBPin,1-positiveDirection)
    d.getFeedback( u3.Timer0Config(TimerMode = 0, Value = (int)((1.0-(lastOutput/100.0))*baseValue)) )
else:
    d.setDOState(pwmBPin,positiveDirection)
    d.setDOState(pwmAPin,1-positiveDirection)
    d.getFeedback( u3.Timer0Config(TimerMode = 0, Value = -(int)((1.0-(lastOutput/100.0))*baseValue)) )

# Update PWM temperature control
def updatePWM(newOutput):
    if newOutput >= 0: # If suggested PID output is non-negative
        # Set pwmAPin and pwmBPin for heating
        d.setDOState(pwmAPin,positiveDirection)
        d.setDOState(pwmBPin,1-positiveDirection)
        # Set PWM duty cycle based on PID output magnitude (how hard Peltier should work)
        d.getFeedback( u3.Timer0( Value =  (int)((1.0-(newOutput/100.0))*baseValue), UpdateReset = True ) )
    else: 
        # Set pwmAPin and pwmBPin for cooling
        d.setDOState(pwmBPin,positiveDirection)
        d.setDOState(pwmAPin,1-positiveDirection)
        # Set PWM duty cycle based on PID output magnitude (how hard Peltier should work)
        d.getFeedback( u3.Timer0( Value =  -(int)((1.0-(newOutput/100.0))*baseValue), UpdateReset = True ) )

tempLoopIndex=0 # Start temperature control loop
nextTempSwitch=tempCycle[0][1] # Next second to trigger tempCycle step

LEDLoopIndex=0 # Start LED control loop
nextLEDSwitch=LEDCycle[0][1] # Next second to trigger LED step

initialtime=time.time() # Start time for experiment

# Read temperature as a voltage, averaging over multiple thermometer measurements
def readTemp():
    tempAvgTemp = 0
    for i in range(tempAvgRange):
        tempAvgTemp += d.getAIN(tempPin)
        d.getFeedback(u3.WaitShort(5))
    return (tempAvgTemp / tempAvgRange)


# Buffer some data in case of write error
tempAvgs = []
times=[]
dataBuffer=[]
bufferError=False
bufferedTemps=[]
bufferedOutputs=[]

###################################
# THIS IS THE ACTUAL CONTROL LOOP #
###################################


while 1==1: # Program only halted manually...
    tempAvg = readTemp() # Record temperature
    
    times.append(time.time()-initialtime) # Save times of sub-second-interval temperature measurements... this isn't really used for anything
    
    if not bufferError: # If no problems writing to data file, print some info to the console
        print( str(lastCycle) + ': ' + str(tempAvg) + " (goal: " + str(goal) + ") (output: " + str(lastOutput) + ") (od: " + str(lastOD) + ")")
    
    output = p.compute(tempAvg) # Ask PID what to do with temperature control based on thermometer reading
    tempAvgs.append(tempAvg) # Temporarily save temperature value
    
    if lastCycle < p.GetCycle(): # If it's been at least a second, the real business starts...
        if output != lastOutput: # Only bother messing with the LabJack if there's something new...
            updatePWM(output) # Tell the LabJack what the PID said to do for temperature control
        
        tempAvg = 0
        for i in range(len(times)):
            tempAvg += tempAvgs[i]
        tempAvg /= len(times) # Average over all temperature readings since last second
        
        bufferedTemps.append(tempAvg) # Temporarily save temperature averages from this minute
        bufferedOutputs.append(lastOutput) # Temporarily save PID output values from this minute
            
        if p.GetCycle() % 10 == 0: # If it's a new MINUTE, was 60
            lastOD, diodeStd = readOD() # Read the OD, once at the very beginning
            # Ignore this commented out bit...
            '''if not startPump:
                if lastOD < odThreshold:
                    thresholdODs+=1
                    if(thresholdODs>thresholdODsMax-1):
                        startPump=True
                else:
                    thresholdODs=0'''
            if(len(dataBuffer)==600): # If there's 10 hours of data in the buffer...
                dataBuffer.pop(0) # Drop the oldest minute's worth of data
            # Buffer all the data from this minute in case there's an error writing to data file
            ## clear all of this out
            dataBuffer.append([p.GetCycle()/60,lastOD,diodeStd,voltsToDegree(goal),voltsToDegree(np.average(bufferedTemps)),np.average(bufferedTemps),np.std(bufferedTemps),np.average(bufferedOutputs),np.std(bufferedOutputs),tempLoopIndex])
            # Clear the minute-scale buffers
            bufferedTemps=[]
            bufferedOutputs=[]
            
            try: # Try to save all the data in the buffer, or specifically this minute's data
                f = open(path+filetype,'a')
                for i in range(len(dataBuffer)):
                    f.write(str(dataBuffer[i][0])+','+str(dataBuffer[i][1])+','+str(dataBuffer[i][2])+','+str(dataBuffer[i][3])+','+str(dataBuffer[i][4])+','+str(dataBuffer[i][5])+','+str(dataBuffer[i][6])+','+str(dataBuffer[i][7])+','+str(dataBuffer[i][8])+','+str(dataBuffer[i][9])+'\n')
                dataBuffer=[]
                bufferError=False
                f.close()
            except: # We can't write to the data file (is it locked?)! Nothing's being saved (except images)!
                bufferError=True
                print( "Enerror printing to file! Please close any programs using it. Buffering up to 10 hours of data...")
        
        s = p.GetCycle()
        new_temp = temp_f(s)
        new_LED = LED_f(s)
        
        d.setDOState(algaePin,1)
        d.getFeedback([u3.DAC8(1,uMToLED(new_LED))])
        d.setDOState(algaePin,0)
        
        if new_temp != goal:
            p.setSetpoint(degreesToVolts(new_temp))
            goal = degreesToVolts(new_temp)
        


        # Just some book-keeping to ready for next minute of loop...
        lastOutput=output    
        lastCycle = p.GetCycle()
        tempAvgs = []
        times=[]
