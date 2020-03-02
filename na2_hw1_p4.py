#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 1 Problem 4

This code plots the interpolating polynomial for 
f(x) = 1/(1+x^2) on [-5,5] using N equidistant points.
"""
import numpy as np
import matplotlib.pyplot as plt
import pylab

def f(x):
    return 1/(1+x**2)

def p(x,N):
    add = 0
    for i in range(1,N):
        prod = 1
        y = 1/(1+((-5*N+10*i-5)/(N-1))**2)
        for j in range(1,N):
            if j != i:
                prod = prod*(N*x-x+5*N-10*j+5)/(10*i-10*j)
        add = add + y*prod
    return add

oo = np.linspace(-5,5)

plt.figure()
plt.title('N=4')
y = f(oo)
pylab.plot(oo,y, label = 'f(x)')
plt.legend()
z = p(oo,4)
pylab.plot(oo,z, label='p_N(x)')
plt.legend()
z_2 = np.abs(f(oo)-p(oo,4))
pylab.plot(oo,z_2, label = '|f(x)-p_N(x)|')
plt.legend()

plt.figure()
plt.title('N=10')
y = f(oo)
pylab.plot(oo,y, label = 'f(x)')
plt.legend()
z = p(oo,10)
pylab.plot(oo,z, label='p_N(x)')
plt.legend()
z_2 = np.abs(f(oo)-p(oo,10))
pylab.plot(oo,z_2, label = '|f(x)-p_N(x)|')
plt.legend()

plt.figure()
plt.title('N=26')
y = f(oo)
pylab.plot(oo,y, label = 'f(x)')
plt.legend()
z = p(oo,26)
pylab.plot(oo,z, label='p_N(x)')
plt.legend()
z_2 = np.abs(f(oo)-p(oo,26))
pylab.plot(oo,z_2, label = '|f(x)-p_N(x)|')
plt.legend()

err_vals = np.zeros(30)

for i in range(1,30):
    max_err = np.amax(np.abs(f(oo)-p(oo,i)))
    err_vals[i] = max_err
    
plt.figure()
plt.plot(err_vals)
plt.title('Maximum error as a function of N')
plt.xlabel('N')
plt.ylabel('Maximum error')