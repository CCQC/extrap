import numpy as np
from scipy.optimize import curve_fit #if #variables = #datapts this produces an analytic fit yolo
from scipy.special import zeta
#fnclist = [exp,xm3,schwartz4,schwartz6,exp2,expGaus,RZ2,RZ3]

        
def exp(x,a,b,c): #exp 3 point fit (analytic), most often used for uncorrelated energy
    return a + b*np.exp(-c*x)
def exp2(x,a,b):
    return a + b*(x+1)*np.exp(-9*np.sqrt(x))
def expGauss(x,a,b,c):
    return a + b*np.exp(-(x-1)) + c*np.exp(-(x-1)**2)

def xm3(x,a,b): #x^-3 2 point fit (analytic), used for correlation energy
    return a + b/(x**3)
def xm5(x,a,b,c):
    return a + b/(x**3) + c/(x**5)

def schwartz4(x,a,b): 
    return a + b/(x + 0.5)**4
def schwartz6(x,a,b,c):
    return a + b/(x + 0.5)**4 + c/(x + 0.5)**6

#these function differently from the others, and (tml) cant be fit with scipy.optimize
def RZ2(xs,ys,pts=2):
    a = xs[-1]**4*(ys[-1] - ys[-2])
    E = ys[-1] + a*(zeta(4) - sum([l**-4 for l in range(1,xs[-1]+1)]))
    return (E, a), 0,  RZ2
def RZ3(xs,ys,pts=3):
    a = (-(xs[-1]-1)**6 * (ys[-2] - ys[-3]) + xs[-1]**6*(ys[-1] - ys[-2]))/(2*xs[-1] - 1)
    b = -a*xs[-1]**2 + xs[-1]**6*(ys[-1] - ys[-2])
    E = ys[-1] + a*(zeta(4) - sum([l**-4 for l in range(1,xs[-1]+1)])) +\
            b*(zeta(6) - sum([l**-6 for l in range(1,xs[-1]+1)]))
    return (E, a, b), 0, RZ3
lookup = \
{"exp":(exp,3),
 "exp2":(exp2,2),
 "expGauss":(expGauss,3),
 "xm3":(xm3,2),
 "xm5":(xm5,3),
 "schwartz4":(schwartz4,2),
 "schwartz6":(schwartz6,3),
 "RZ2":(RZ2,2),
 "RZ3":(RZ3,3)} 

def extrapolate(xs,ys,extrapfnc,pts=2,bound=None):
    xmod = xs[len(xs)-pts:] #cuts down zeta vals to those used in extrap
    ymod = ys[len(ys)-pts:] #"" for energies
    #print(xmod)
    #print(ymod)
    #print(extrapfnc)
    if extrapfnc == RZ2 or extrapfnc == RZ3:
        return extrapfnc(xs,ys,pts=pts)
    if bound:
        fit,fitcov = curve_fit(extrapfnc,xmod,ymod,bounds=bound)
    else:
        fit,fitcov = curve_fit(extrapfnc,xmod,ymod)
    return fit,fitcov,extrapfnc

def errest(xs,ys,pts=2):
    ps = []
    for entry in lookup:
        if lookup[entry][1] <= pts:
            if "RZ" in entry:
                myextrap = extrapolator(lookup[entry][0],fit="manual",pts=lookup[entry][1])
            else:
                myextrap = extrapolator(lookup[entry][0],pts=lookup[entry][1])
            myextrap.extrapolate(xs,ys)
            ps.append([entry,myextrap.p[0]])
    return np.array(ps)



class extrapolator:
    def __init__(self,fnc,pts=2,fit=extrapolate):
        self.fit = fit
        self.pts = pts
        if type(fnc) == str:
            if fnc in lookup:
                self.fnc = lookup[fnc][0]
            else:
                pass
        else:
            self.fnc = fnc
        
    def extrapolate(self,xs,ys):
        if self.fit == extrapolate:
            p,cov,fnc = extrapolate(xs,ys,self.fnc,pts=self.pts)
        elif self.fit == "m" or self.fit == "man" or self.fit == "manual":
            p,cov,fnc = self.fnc(xs,ys,pts=self.pts)
        else:
            p,cov,fnc = (None,None,None)
        self.p = p
        self.cov = cov
        assert self.fnc == fnc
        
