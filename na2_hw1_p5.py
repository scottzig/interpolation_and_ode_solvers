#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 1 Problem 5

This code plots the piecewise linear interpolating polynomial for 
f(x) = 1/(1+x^2) on [-5,5] using N equidistant points.
"""
import numpy as np
import matplotlib.pyplot as plt
import pylab

def f(x):
    return 1/(1+x**2)

def x(i,N):
    return (-5*N+10*i-5)/(N-1)

def y(i,N):
    return 1/(1+((-5*N+10*i-5)/(N-1))**2)

def p(x_1,x_2,y_1,y_2,x):
    return y_1+((y_2-y_1)/(x_2-x_1))*(x-x_1) 

plt.figure()
plt.title('N=5')
N = 5
xx = np.linspace(-5,5)
plt.plot(xx,f(xx))
for i in range(1,N):
    oo = np.linspace(x(i,N),x(i+1,N))
    fun = p(x(i,N),x(i+1,N),y(i,N),y(i+1,N),oo)
    diff = np.abs(f(oo)-fun)
    plt.plot(oo,fun)
    plt.plot(oo,diff)
    
plt.figure()
plt.title('N=20')
N = 20
xx = np.linspace(-5,5)
plt.plot(xx,f(xx))
for i in range(1,N):
    oo = np.linspace(x(i,N),x(i+1,N))
    fun = p(x(i,N),x(i+1,N),y(i,N),y(i+1,N),oo)
    diff = np.abs(f(oo)-fun)
    plt.plot(oo,fun)
    plt.plot(oo,diff)
    
plt.figure()
plt.title('N=45')
N = 45
xx = np.linspace(-5,5)
plt.plot(xx,f(xx))
for i in range(1,N):
    oo = np.linspace(x(i,N),x(i+1,N))
    fun = p(x(i,N),x(i+1,N),y(i,N),y(i+1,N),oo)
    diff = np.abs(f(oo)-fun)
    plt.plot(oo,fun)
    plt.plot(oo,diff)

err_vals = np.zeros(1000)

for i in range(1,1000):
    err_vals_sub = 0
    for j in range(1,i):
        oo = np.linspace(x(j,N),x(j+1,N))
        fun = p(x(j,N),x(j+1,N),y(j,N),y(j+1,N),oo)
        pot_err = np.amax(np.abs(f(oo)-fun))
        if pot_err > err_vals_sub:
            err_vals_sub = pot_err
    err_vals[i] = err_vals_sub
    
plt.figure()
plt.plot(err_vals)
plt.title('Maximum error as a function of N')
plt.xlabel('N')
plt.ylabel('Maximum error')
