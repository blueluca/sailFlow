import math

def profile(xo, h, absx, absy, miny, maxy):
#=====================================================================================
# The profile calculator
#
# inputs 
#   x  : point in the range 0 to 1 to calculate
#   h :  height in percentage
#   absx : absolute position in x
#   absy : absolute position in y
#   miny : min y
#   maxy : max y
#=====================================================================================
# NACA numbers mptt
# m max camber in percentage
# p pos camber in percentage
# t airfoil thickness
#
# NACA 2415

    m0 = 0.1
    p0 = 0.5
    c = 1
    t = 0.1
    error = 0.001
    maxcount = 100

    xu = 0
    x = xo-0.001
    
    if (x < 0):
        return 0

    m = m0
    p = p0    
    st = 0.4
    if h > st:
        p = p0+ 0.5*(1/(1-st)*(h-st))
        m = m0 - 0.099*(1/(1-st)*(h-st))
        
    print("p=",round(p,2)," m=",round(m,2))
 
#    m = 0.08
    count = 0
    previousErr = 1.1
    d = -0.001
    while True:
        count += 1 
        if (x < p):
            yc = (m / (p * p)) * (2 * p * x - x * x);
            yc1 = (2*m/p)*(p-x)
        else:
            yc = (m / ((1 - p) * (1 - p))) * ((1 - 2 * p) + 2 * p * x - x * x);
            yc1 = (2*m/((1-p)*(1-p)))*(p-x)
    
        yt = 5*t*c*(0.2969*math.sqrt(x/c)-0.1260*(x/c)-0.3516*math.pow(x/c,2)+0.2843*math.pow(x/c,3)-0.1015*math.pow(x/c,4));
        
        xu = x -yt * math.sin(math.atan(yc1))
        yu = yc + yt * math.cos(math.atan(yc1))
#        print("Delta xo-xu=", round(abs(xo-xu),4))
        if abs(xo-xu)>previousErr:
            d = -d
        previousErr = abs(xo-xu)
        if previousErr<error:
            break
        if count > maxcount:
            print("MAXCOUNT !")
            break
        x += d
        
#    print(("Profile return Z=",round(yu,4)," x=",round(xo,4)," xu=",round(xu,4)))
#    print("Exit Angle=",math.atan(yc1)*180/math.pi)
    return yu