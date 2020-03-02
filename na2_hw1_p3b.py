#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 1 Problem 3b

This code plots the interpolating polynomial for 
f(x) = 1/x on [1/2,3/2] using N points determined by
the roots of the Chebyshev polynomial.
"""
import numpy as np
import matplotlib.pyplot as plt
import pylab

def f(x):
    return 1/x

def p(x,N):
    add = 0
    for i in range(1,N):
        prod = 1
        y = 2/(2+np.cos((2*i-1)*np.pi/(2*N)))
        for j in range(1,N):
            if j != i:
                prod = prod*(2*x-2-np.cos((2*j-1)*np.pi/(2*N)))/(np.cos((2*i-1)*np.pi/(2*N))-np.cos((2*j-1)*np.pi/(2*N)))
        add = add + y*prod
    return add

oo = np.linspace(1/2,3/2)

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
plt.title('N=15')
y = f(oo)
pylab.plot(oo,y, label = 'f(x)')
plt.legend()
z = p(oo,15)
pylab.plot(oo,z, label='p_N(x)')
plt.legend()
z_2 = np.abs(f(oo)-p(oo,15))
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