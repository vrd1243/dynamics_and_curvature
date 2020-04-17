import numpy as np
from lorenz import get_lorenz
from rossler import get_rossler
from pendulum import derivative, forced_derivative, changeA
from rk4 import rk4

def data_lorenz(noise_level,n=50000,tmax=500):

    [t, data] = get_lorenz(n=n, tmax=tmax)
    max_val = np.max(np.abs(data[:,0]))
    series = data[:,0] + noise_level*max_val*(np.random.random((data.shape[0])) - .5)

    return series

def data_rossler(noise_level):

    [t, data] = get_rossler(n=50000, tmax=10000);

    max_val = np.max(np.abs(data[:,0]))
    series = data[:,0] + noise_level*max_val*(np.random.random((data.shape[0])) - .5); 

    return series;

def data_pendulum(noise_level):
    
    changeA(0.886); 
    [data,t] = rk4(derivative, [3.14, 50], 0, 0.005, 21000); 
    data = data[1000:,:];
    max_val = np.max(np.abs(data[:,0]))
    series = data[:,0] + noise_level*max_val*(np.random.random((data.shape[0])) - .5); 

    return series;

def data_pendulum_omega(noise_level):
    
    changeA(0.886); 
    [data,t] = rk4(derivative, [3.14, 50], 0, 0.005, 21000); 
    data = data[1000:,:];
    max_val = np.max(np.abs(data[:,1]))
    series = data[:,1] + noise_level*max_val*(np.random.random((data.shape[0])) - .5); 

    return series;

def data_musical(noise_level): 
      
    data = np.loadtxt('c_maj_scale/02- C scale-consolidated.txt')
    data = data[40000:50000]
    
    max_val = np.max(data) - np.min(data)
    series = data + noise_level*max_val*(np.random.random((data.shape[0])) - .5); 
    
    return series;

def data_pendulum_portrait(noise_level):
   
    data = np.loadtxt('pendulum_portrait.txt'); 
    print(data[:,0].shape)
    idx = np.where(np.logical_and(data[:,0] < 10, data[:,0] > -10))
    data = data[idx, 0].reshape(-1,1)
    print(data[:,0].shape)
    max_val = np.max(np.abs(data[:,0]))
    series = data[:,0] + noise_level*max_val*(np.random.random((data.shape[0])) - .5); 

    return series;
