#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 1 Problem 5

This code plots the piecewise quadratic interpolating polynomial for 
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

def p(x_1,x_2,x_3,y_1,y_2,y_3,x):
    return y_1*((x-x_2)*(x-x_3))/((x_1-x_2)*(x_1-x_3)) \
            +y_2*((x-x_1)*(x-x_3))/((x_2-x_1)*(x_2-x_3)) \
            +y_3*((x-x_1)*(x-x_2))/((x_3-x_1)*(x_3-x_2)) 

plt.figure()
plt.title('N=5')
N = 5
xx = np.linspace(-5,5)
plt.plot(xx,f(xx))
for i in range(1,N):
    if i % 2 != 0:
        oo = np.linspace(x(i,N),x(i+2,N))
        fun = p(x(i,N),x(i+1,N),x(i+2,N),y(i,N),y(i+1,N),y(i+2,N),oo)
        diff = np.abs(f(oo)-fun)
        plt.plot(oo,fun)
        plt.plot(oo,diff)
    
plt.figure()
plt.title('N=21')
N = 21
xx = np.linspace(-5,5)
plt.plot(xx,f(xx))
for i in range(1,N):
    if i % 2 != 0:
        oo = np.linspace(x(i,N),x(i+2,N))
        fun = p(x(i,N),x(i+1,N),x(i+2,N),y(i,N),y(i+1,N),y(i+2,N),oo)
        diff = np.abs(f(oo)-fun)
        plt.plot(oo,fun)
        plt.plot(oo,diff)
    
plt.figure()
plt.title('N=45')
N = 45
xx = np.linspace(-5,5)
plt.plot(xx,f(xx))
for i in range(1,N):
    if i % 2 != 0:
        oo = np.linspace(x(i,N),x(i+2,N))
        fun = p(x(i,N),x(i+1,N),x(i+2,N),y(i,N),y(i+1,N),y(i+2,N),oo)
        diff = np.abs(f(oo)-fun)
        plt.plot(oo,fun)
        plt.plot(oo,diff)

err_vals = np.zeros(1000)

for i in range(1,1000):
    if i % 2 != 0:
        err_vals_sub = 0
        for j in range(1,i):
            if j % 2 != 0:
                oo = np.linspace(x(i,N),x(i+2,N))
                fun = p(x(i,N),x(i+1,N),x(i+2,N),y(i,N),y(i+1,N),y(i+2,N),oo)
                pot_err = np.amax(np.abs(f(oo)-fun))
                if pot_err > err_vals_sub:
                    err_vals_sub = pot_err
        err_vals[i] = err_vals_sub
    
plt.figure()
plt.plot(err_vals)
plt.title('Maximum error as a function of N')
plt.xlabel('N')
plt.ylabel('Maximum error')
