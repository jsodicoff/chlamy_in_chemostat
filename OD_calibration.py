def temp_f(s):
    if s < 10*60:
      return(25)
    if s < 5*60:
      return(30)
    elif s < 2*5*60:
      return(35)
    elif s < 3*5*60:
      return(30)
    elif s < 4*5*60:
      return(40)
    elif s < 5*5*60:
      return(30)
    elif s < 6*5*60:
      return(40)
    else:
      return(40)
        
def LED_f(s):
    return(0)
