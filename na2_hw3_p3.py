#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 3 Problem 3

This code approximates the integral of sin(1/x) on 0.05 to 1
using a trapezoidal rule where the intervals are chosen using an
adaptive method.
"""
import numpy as np
import matplotlib.pyplot as plt

true_val = 0.50283962;

def f(x):
    return np.sin(1/x)

def nu(x_f,x_l):
    h = x_l-x_f
    mid = (x_l+x_f)/2
    der_2 = (2/(mid**3))*np.cos(1/mid)-(1/(mid**4))*(np.sin(1/mid))
    return 1/2*np.abs(der_2)*(h**3)
    
# Adaptive trapezoidal Rule 

plt.figure()  
K = 1
err = 1
eval_pts = (0.05,1)
while err > 10**(-6):
    sum = 0
    err_sub = 0
    replace = 0
    for j in range(0,K):
        x_1 = eval_pts[j]
        x_2 = eval_pts[j+1]
        sum = sum + (x_2-x_1)*(f(x_1)+f(x_2))/2
        err_new = nu(x_1,x_2)
        if err_new > err_sub:
            err_sub = err_new
            replace = j
    err = np.abs(true_val - sum)
    plt.plot(K+1,np.log10(err),'.')
    plt.title('Trapezoidal Rule with Adaptive Intervals')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Error (log scale)')
    r_1 = replace
    r_2 = replace + 1
    new_eval_pt = (eval_pts[r_1]+eval_pts[r_2])/2
    eval_pts = np.append(eval_pts,new_eval_pt)
    eval_pts = np.sort(eval_pts)
    K = K+1
    
txt = str(K+2)
print('Adaptive trapezoidal function evaluations required = ' + txt)