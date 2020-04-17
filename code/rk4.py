#!/bin/python
import numpy as np

def rk4(derivative, x0, t0, h, n):
    
    x = np.zeros((n, len(x0)));
    print(x.shape)
    t = np.zeros(n);
    x[0,:] = x0;
    t[0] = t0;
    for i in range(1,n):
        k1 = derivative(x[i-1,:], t[i-1]);
        k2 = derivative(x[i-1,:] + h*k1/2, t[i-1] + h/2);        
        k3 = derivative(x[i-1,:] + h*k2/2, t[i-1] + h/2);        
        k4 = derivative(x[i-1,:] + h*k3, t[i-1] + h);
        x[i,:] = x[i-1,:] + h/6*(k1 + 2*k2 + 2*k3 + k4);        
        t[i]  = t[i-1] + h

    #return [x[int(n/10):,:],t[int(n/10):]]
    return [x, t]
