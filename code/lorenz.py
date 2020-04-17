import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from rk4 import rk4

# Lorenz paramters and initial conditions
u0, v0, w0 = 0, 1, 1.05

#def lorenz(X, t, sigma, beta, rho):
def lorenz(X, t):
    sigma, beta, rho = 10, 2.667, 28
    dx = np.zeros(3);
    """The Lorenz equations."""
    u, v, w = X
    dx[0] = -sigma*(u - v)
    dx[1] = rho*u - v - u*w
    dx[2] = -beta*w + u*v
    
    #return up, vp, wp
    return dx

def get_lorenz(n=10000,timestep=0.01):
# Integrate the Lorenz equations on the time grid t
  #t = np.linspace(0, tmax, n+1)
  #t = np.arange(n)*timestep
  #f = odeint(lorenz, (u0, v0, w0), t, args=(sigma, beta, rho))
  #x, y, z = f.T
  
  [data, t] = rk4(lorenz, [u0, v0, w0], 0, timestep, n);  

  #array = np.concatenate((np.reshape(x,(-1,1)), np.reshape(y,(-1,1)),np.reshape(z,(-1,1))), axis=1);
  #np.savetxt('lorenz.txt', array);
  return [np.reshape(t,-1,1)[int(n/10):], data[int(n/10):,:]];

