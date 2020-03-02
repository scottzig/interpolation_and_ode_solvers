#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 4 Problem 4

This code computes solutions to the differential equation
x' = (1/2)x w/ initial condition x(0)=1 using four different methods
with a variety of step sizes.
"""
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return (1/2)*x

def x_true(t):
    return np.exp(1/2*t)

# first order Taylor expansion method
delt_t = 2
tf_vec = np.zeros(7)
i = 0
while delt_t >= 1/32:
    t = 0
    x = 1
    while t < 4:
        x_new = (1+delt_t/2)*x
        x = x_new
        t = t+delt_t
    err = np.abs(x-x_true(4))
    tf_vec[i] = err
    plt.loglog(delt_t,err,'ko')
    plt.loglog(delt_t,delt_t,'ro')
    plt.title('First order Taylor in black, h^1 in red')
    plt.xlabel('delt_t')
    plt.ylabel('Error')
    delt_t = delt_t/2
    i = i+1
    
# second order Taylor expansion method
plt.figure()
delt_t = 2
ts_vec = np.zeros(7)
i = 0
while delt_t >= 1/32:
    t = 0
    x = 1
    while t < 4:
        x_new = x + delt_t*f(x) + delt_t**2*(1/4)*f(x)
        x = x_new
        t = t+delt_t
    err = np.abs(x-x_true(4))
    ts_vec[i] = err
    plt.loglog(delt_t,err,'ko')
    plt.loglog(delt_t,delt_t**2,'ro')
    plt.title('Second order Taylor in black, h^2 in red')
    plt.xlabel('delt_t')
    plt.ylabel('Error')
    delt_t = delt_t/2
    i = i+1
    
# implicit Euler method
plt.figure()
delt_t = 1
ie_vec = np.zeros(7)
i = 1
while delt_t >= 1/32:
    t = 0
    x = 1
    while t < 4:
        x_new = x/(1-(1/2)*delt_t)
        x = x_new
        t = t+delt_t
    err = np.abs(x-x_true(4))
    ie_vec[i] = err
    plt.loglog(delt_t,err,'ko')
    plt.loglog(delt_t,delt_t,'ro')
    plt.title('Implicit Euler in black, h^1 in red')
    plt.xlabel('delt_t')
    plt.ylabel('Error')
    delt_t = delt_t/2
    i = i+1
    
# crank-nicolson
plt.figure()
delt_t = 2
cn_vec = np.zeros(7)
i = 0
while delt_t >= 1/32:
    t = 0
    x = 1
    while t < 4:
        x_new = (x+(1/4)*delt_t*x)/(1-(1/4)*delt_t)
        x = x_new
        t = t+delt_t
    err = np.abs(x-x_true(4))
    cn_vec[i] = err
    plt.loglog(delt_t,err,'ko')
    plt.loglog(delt_t,delt_t**2,'ro')
    plt.title('Crank-Nicolson in black, h^2 in red')
    plt.xlabel('delt_t')
    plt.ylabel('Error')
    delt_t = delt_t/2
    i = i+1
    
plt.figure()
pt_vec = [2,1,1/2,1/4,1/8,1/16,1/32]
plt.semilogy(pt_vec,tf_vec,'bo',label = 'First order Taylor')
plt.semilogy(pt_vec,ts_vec,'co',label = 'Second order Taylor')
plt.semilogy(pt_vec,ie_vec,'yo',label = 'Implicit Euler')
plt.semilogy(pt_vec,cn_vec,'mo',label = 'Crank-Nicolson')  
plt.title('Comparison of accuracy of ODE solvers') 
plt.xlabel('delt_t')
plt.ylabel('Error')
plt.legend()