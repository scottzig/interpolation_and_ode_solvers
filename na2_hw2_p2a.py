#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 2 Problem 2a

This code plots the interpolating polynomial for 
f(x) = sin(1/x) on [0.05,1] with a piecewise linear
function over an increasing number of intervals
"""
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.sin(1/x)

def p(x,x_1,x_2):
    y_1 = np.sin(1/x_1)
    y_2 = np.sin(1/x_2)
    y = y_1 + (y_2-y_1)/(x_2-x_1)*(x-x_1)
    return y

N = 129

err_vals = np.zeros(N-2)

for i in range(2,N):
    err_sub = 0
    for j in range(0,i-1):
        pt_1 = 0.05+(1-0.05)*j/(i-1)
        pt_2 = 0.05+(1-0.05)*(j+1)/(i-1)
        oo = np.linspace(pt_1,pt_2,100)
        err_new = np.max(np.abs(f(oo)-p(oo,pt_1,pt_2)))
        if err_new > err_sub:
            err_sub = err_new
    err_vals[i-2] = err_sub
    
plt.figure()
plt.title('Number of interpolation points = 129')
yy = np.linspace(0.05,1,100)
plt.plot(yy,f(yy))
N_2 = 129
for j in range(0,N_2-1):
    pt_1 = 0.05+(1-0.05)*j/(N_2-1)
    pt_2 = 0.05+(1-0.05)*(j+1)/(N_2-1)
    oo = np.linspace(pt_1,pt_2,100)
    plt.plot(oo,p(oo,pt_1,pt_2))
    print(pt_1)
    if j == N_2-2:
        print(pt_2)
    
plt.figure()
plt.plot(err_vals)
plt.title('Maximum Error as a function of number of evaluation points')   