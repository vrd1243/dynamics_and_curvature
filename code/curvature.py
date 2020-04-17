#!/bin/python

import numpy as np
from numpy import linalg as LA
from pytisean import tiseanio

def menger_curvature(x,y,z):
    a = LA.norm(x-y);
    b = LA.norm(y-z);
    c = LA.norm(x-z);

    s = (a + b + c) / 2;

    if (a*b*c == 0):
        print("Hitting infinite curvature");
        return np.inf;
    
    if (s-a)*(s-b)*(s-c) < 0:
        print("Semi-perimeter is less than side", s, a, b, c);
        return 0;

    return 4*np.sqrt(s*(s-a)*(s-b)*(s-c)) / (a*b*c);

def get_menger(series, tau):
    
    menger_series = [];
    print(1, len(series) - tau - 1);
    for t in range(1, len(series) - tau - 1): 
        menger_series.append(menger_curvature(np.asarray([series[t-1], series[t+tau-1]]),
                                              np.asarray([series[t], series[t+tau]]),
                                              np.asarray([series[t+1], series[t+tau+1]])));
    
     
    return np.asarray(menger_series); 

def get_menger_perturb_avg(series, tau):
    
    menger_series = [];
    scale = 0.1;
    for t in range(1, len(series) - tau - 1): 
        
        menger_avg = []; 
        for i in range(100):
            first = np.asarray([series[t-1], series[t+tau-1]]) + scale*np.random.random((2,1))
            second = np.asarray([series[t], series[t+tau]]) + scale*np.random.random((2,1))
            third = np.asarray([series[t+1], series[t+tau+1]]) + scale*np.random.random((2,1))

            menger_avg.append(menger_curvature(first,
                                              second,
                                              third));
    
        menger_series.append(np.sum(menger_avg) / 100);
    
    return np.asarray(menger_series); 


def get_menger_centered_avg_embed(embed, tau_avg = 1, n_neighbors = 4):

    menger_series = []; 
    trange = range(tau_avg + n_neighbors, embed.shape[0] - tau_avg - n_neighbors)
    for t in trange: 
        menger_avg = []; 
        for k in range(-n_neighbors, n_neighbors+1):
            
            first = embed[t - tau_avg + k]; 
            second = embed[t + k];
            third = embed[t + tau_avg + k];
                                                 
            menger_avg.append(menger_curvature(first,
                                              second,
                                              third));
         
        menger_series.append(np.sum(menger_avg) / (2*n_neighbors + 1));
    
    return np.asarray(menger_series); 

def get_menger_centered_avg(series, tau, tau_avg = 1, n_neighbors = 4):

    menger_series = [];

    #print(tau_avg + n_neighbors, len(series) - tau - tau_avg - n_neighbors);
    
    for t in range(tau_avg + n_neighbors, len(series) - tau - tau_avg - n_neighbors): 
        menger_avg = []; 
        for k in range(-n_neighbors, n_neighbors+1):
            first = np.asarray([series[t-tau_avg+k], series[t+tau-tau_avg+k]]);
            second = np.asarray([series[t+k], series[t+tau+k]]);
            third = np.asarray([series[t+tau_avg+k], series[t+tau+tau_avg+k]]);
    
            menger_avg.append(menger_curvature(first,
                                              second,
                                              third));
         
        menger_series.append(np.sum(menger_avg) / (2*n_neighbors + 1));
    
    return np.asarray(menger_series); 

def get_menger_new_centered_avg(series, tau, tau_avg = 1, n_neighbors = 4):

    menger_series = [];

    print(tau_avg + n_neighbors, len(series) - tau - tau_avg - n_neighbors);
    for t in range(tau_avg + n_neighbors, len(series) - tau - tau_avg - n_neighbors): 
        menger_avg = []; 
        for k in range(-n_neighbors, n_neighbors+1):
           
           reach_avg = 0;

           for reach in range(1, tau_avg + 1):
                
              first = np.asarray([series[t-reach+k], series[t+tau-reach+k]]);
              second = np.asarray([series[t+k], series[t+reach+k]]);
              third = np.asarray([series[t+reach+k], series[t+tau+reach+k]]);
              reach_avg += menger_curvature(first, second, third)
           
           menger_avg.append(reach_avg / tau_avg);
         
        menger_series.append(np.sum(menger_avg) / (2*n_neighbors + 1));

    return np.asarray(menger_series); 
