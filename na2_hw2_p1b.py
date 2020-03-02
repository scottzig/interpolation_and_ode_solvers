#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 2 Problem 1a

This code plots the interpolating polynomial for 
f(x) = 0 on [-1,1] using N points chosen as
the roots of the Chebyshev polynomial with data
y_i corrupted by noise.
"""
import numpy as np
import matplotlib.pyplot as plt
import pylab

oo = np.linspace(-1,1,100)

def poly(x):
    return x-x

def p(x,N):
    y = np.zeros(N)
    for i in range(0,N-1):
        y[i] = np.random.uniform(-0.001,0.001)
    add = 0
    for i in range(1,N):
        prod = 1
        for j in range(1,N):
            if j != i:
                prod = prod*(x-np.cos((2*j-1)*np.pi/(2*N)))/(np.cos((2*i-1)*np.pi/(2*N))-np.cos((2*j-1)*np.pi/(2*N)))
        add = add + y[i]*prod
    return add

for i in range(3,31):
    plt.figure()
    num = str(i-1)
    plt.title(num + ' points')
    z_1 = poly(oo)
    pylab.plot(oo,z_1, label = 'p_{N-1}')
    plt.legend()
    z_2 = p(oo,i)
    pylab.plot(oo,z_2, label='tilde{p}_{N-1}(x)')
    plt.legend()
    z_3 = np.abs(z_2-z_1)
    pylab.plot(oo,z_3, label = 'Error')
    plt.legend()
    
err_vals = np.zeros(30)
    
for i in range(0,30):
    max_err = np.amax(np.abs(poly(oo)-p(oo,i)))
    err_vals[i] = max_err
    
plt.figure()
plt.plot(err_vals)
plt.title('Maximum error as a function of N')
plt.xlabel('N')
plt.ylabel('Maximum error')