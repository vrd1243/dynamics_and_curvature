#!/bin/python    
import numpy as np

beta = 0.25
m = 0.1
l = 0.1
g = 9.8
A = 0.75
alpha = np.sqrt(g/l)*0.75

def derivative(x, t):
    dx = np.zeros(2);
    dx[0] = x[1];
    dx[1] = A*np.cos(alpha*t) / (m*l) - beta/m*x[1] - g/l*np.sin(x[0])
    return dx

def changeA(a):
    global A;
    A = a;


def forced_derivative(x, t):
    
    dx = np.zeros(2)
    dx[0] = x[1]
    dx[1] = -0.3*np.sin(x[0]) + 2*np.sin((2*x[0] - t))

    return dx;
