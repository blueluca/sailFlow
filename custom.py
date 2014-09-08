def profile(x, h, x, y):
#=====================================================================================
# The profile calculator
#
# inputs 
#   x  : point in the range 0 to 1 to calculate
#   m : camber percentage
#   p : camber position percentage
#   x : absolute position in x
#   y : absolute position in y
#=====================================================================================
    m = 0.15
    p = 0.5
    
    if (x < 0):
        return 0
    if (x < p):
        y = (m / (p * p)) * (2 * p * x - x * x);
    else:
        y = (m / ((1 - p) * (1 - p))) * ((1 - 2 * p) + 2 * p * x - x * x);
    # debug(("Profile return Z="+str(y))
    return y