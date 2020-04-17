import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Rossler paramters and initial conditions
x0, y0, z0 = 0, 0, 0
a, b, c = 0.13, 0.2, 6.5

def rossler(X, t, a, b, c):
    """The Lorenz equations."""
    x, y, z = X
    xp = -y - z
    yp = x + a*y
    zp = b + z*(x - c)
    return xp, yp, zp

def get_rossler(tmax=100, n=10000, start=[0,0,0]):
# Integrate the Lorenz equations on the time grid t
  
  x0 = start[0];
  y0 = start[1];
  z0 = start[2];

  t = np.linspace(0, tmax, n+1)
  f = odeint(rossler, (x0, y0, z0), t, args=(a, b, c))
  x, y, z = f.T

  array = np.concatenate((np.reshape(x,(-1,1)), np.reshape(y,(-1,1)),np.reshape(z,(-1,1))), axis=1);
  np.savetxt('rossler.txt', array);
  return [np.reshape(t,-1,1)[int(n/5):], array[int(n/5):]];
