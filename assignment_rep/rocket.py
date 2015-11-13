# Script to solve rocket assignment for pure vertical flight of module 1.

from math import sin, cos, log, ceil
import numpy
from matplotlib import pyplot
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16


# model parameters:
g = 9.81        # gravity in m s^{-2}   
C_D = 0.15      # drag coefficient --- or D/L if C_L=1
ms = 50.        # rocket shell mass in kg
ro = 1.091      # average air density kg/m^3 (constant through flight)
r = 0.5         # radius of rocket in m
ve = 325.       # exhaust speed in m/s
A = numpy.pi*r**2 # diameter of rocket

### set initial conditions ###
v0 = 0.
mp0 = 100.     # mass of rocket propellant at time t=0
x0 = 0.        # horizotal position is arbitrary and constant
y0 = 0.        # initial altitude


def f(u,t):
    """Returns the right-hand side of the phugoid system of equations.
    
    Parameters
    ----------
    u : array of float
        array containing the solution at time n.
        
    Returns
    -------
    dudt : array of float
        array containing the RHS given u.
    """
    
    v = u[0]
    mp = u[1]
    y = u[2]
    
    if t<=5:
        return numpy.array([-g + 20.*ve/(ms+mp) - 0.5/(ms+mp)*ro*v*abs(v)*A*C_D,
                      -20.,
                      v])
    else:
        return numpy.array([-g - 0.5/ms*ro*v*abs(v)*A*C_D,
                      0.,
                      v])                  
    
                      
def euler_step(u, f, dt):
    """Returns the solution at the next time-step using Euler's method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equation.
    dt : float
        time-increment.
    
    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """
    
    return u + dt * f(u,t)
    
                          

T = 100                          # final time
dt = 0.1                         # time increment
N = int(T/dt) + 1                # number of time-steps
t = 0.
# initialize the array containing the solution for each time-step
u = numpy.empty((N, 3))
u[0] = numpy.array([v0, mp0, y0])# fill 1st element with initial values

# time loop - Euler method
for n in range(N-1):
        
        t += dt
        
        u[n+1] = euler_step(u[n], f, dt) 
        
        print "speed of rocket %s" % u[n+1][0]
        
        if u[n+1][2]<0:
            break

u = numpy.delete(u,numpy.arange(n+2,len(u)),axis=0)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            