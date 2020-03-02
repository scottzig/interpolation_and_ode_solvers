#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 2 Problem 2b

This code plots the interpolating polynomial for 
f(x) = sin(1/x) on [0.05,1] with a piecewise linear
function over an increasing number of intervals where we 
only replace intervals with bad approximations.
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

int_pts = (0.05,1)

for i in range(2,N):
    #txt = str(i)
    #plt.figure()
    #plt.title('Number of interpolation points = ' + txt)
    #yy = np.linspace(0.05,1,100)
    #plt.plot(yy,f(yy))
    err_sub = 0
    replace = 0
    if i == N-1:
        txt = str(i+1)
        plt.figure()
        plt.title('Number of interpolation points = ' + txt)
        yy = np.linspace(0.05,1,100)
        plt.plot(yy,f(yy))
    for j in range(0,i-1):
        pt_1 = int_pts[j]
        pt_2 = int_pts[j+1]
        oo = np.linspace(pt_1,pt_2,100)
        if i == N-1:
            plt.plot(oo,p(oo,pt_1,pt_2))
        err_new = np.max(np.abs(f(oo)-p(oo,pt_1,pt_2)))
        if err_new > err_sub:
            err_sub = err_new
            replace = j
    if i < N-1:
        err_vals[i-2] = err_sub
        r_1 = replace
        r_2 = replace + 1
        new_int_pt = (int_pts[r_1]+int_pts[r_2])/2
        int_pts = np.append(int_pts,new_int_pt)
        int_pts = np.sort(int_pts)

for i in range(0,128):
    print(int_pts[i])
    
plt.figure()
plt.plot(err_vals)
plt.title('Maximum Error as a function of number of evaluation points')