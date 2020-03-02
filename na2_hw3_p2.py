#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 3 Problem 1

This code approximates the integral of sin(1/x) on 0.05 to 1
using the Romberg algorithm
"""
import numpy as np
import matplotlib.pyplot as plt

true_val = 0.50283962;

def f(x):
    return np.sin(1/x)

# Romberg Algorithm
    
plt.figure()
err = 1
a = 0.05
b = 1
M = 1
# col = 3
while err > 10**(-6):
    h = b-a
    R = np.zeros([M,M])
    R_val = h*0.5*(f(a)+f(b))
    R[0,0] = R_val
    if M > 1:
        for n in range(1,M):
            h = h/2
            stop = 2**(n-1)+1
            sum = 0
            for i in range(1,stop):
                sum = sum + f(a+(2*i-1)*h)
            R[n,0] = 0.5*R[n-1,0] + h*sum
            for m in range(1,n+1):
                R[n,m] = R[n,m-1] + (R[n,m-1]-R[n-1,m-1])/(4**m-1)
    #if M > col:
    err = np.abs(R[M-1,M-1] - true_val)
    plt.plot(2**M+1,np.log10(err),'.')
    plt.title('Romberg Algorithm')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Error (log scale)')
    M = M+1
    
txt = str(2**M+2)
print('Romberg function evaluations required = ' + txt)