#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 4 Problem 3b

This code approximates the minimum of f(x) = xe^x-1 by using
Newton's method on the function g(x)=f'(x) and approximating
f'(x) and f''(x) using symmetric finite difference operators.
"""
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x*np.exp(x)-1

def g(x,h):
    return float(f(x+h)-f(x-h))/float(2*h)

def der_g(x,h):
    return float(f(x+h)-2*f(x)+f(x-h))/float(h**2)

x_0 = -0.5
err = 100
x = x_0
h_0 = 10**(0)
while err > 10**(-6):
    top = g(x,h_0)
    bottom = der_g(x,h_0)
    if bottom == 0.0:
        break
    x_new = x - top/bottom
    err = np.abs(x-x_new)
    x = x_new
txt_0 = str(x)
print('The minimum using h=1 is at point ' + txt_0)
    
x_0 = -0.5
err = 100
x = x_0
h_2 = 10**(-2)
while err > 10**(-6):
    top = g(x,h_2)
    bottom = der_g(x,h_2)
    if bottom == 0.0:
        break
    x_new = x - top/bottom
    err = np.abs(x-x_new)
    x = x_new
    
txt_2 = str(x)
print('The minimum using h=10^(-2) is at point ' + txt_2)
    
x_0 = -0.5
err = 100
x = x_0
h_4 = 10**(-4)
while err > 10**(-6):
    top = g(x,h_4)
    bottom = der_g(x,h_4)
    if bottom == 0.0:
        break
    x_new = x - top/bottom
    err = np.abs(x-x_new)
    x = x_new
    
txt_4 = str(x)
print('The minimum using h=10^(-4) is at point ' + txt_4)

x_0 = -0.5
err = 100
x = x_0
h_6 = 10**(-6)
while err > 10**(-6):
    top = g(x,h_6)
    bottom = der_g(x,h_6)
    if bottom == 0.0:
        break
    x_new = x - top/bottom
    err = np.abs(x-x_new)
    x = x_new
    
txt_6 = str(x)
print('The minimum using h=10^(-6) is at point ' + txt_6)

h = 1
while h > 10**(-6):
    x_0 = -0.5
    err = 100
    x = x_0
    while err > 10**(-6):
        top = g(x,h)
        bottom = der_g(x,h)
        if bottom == 0.0:
            break
        x_new = x - top/bottom
        err = np.abs(x-x_new)
        x = x_new
    err_1 = np.abs(x+1)
    plt.loglog(h,err_1,'.')
    plt.title('Root Approximation')
    plt.xlabel('h')
    plt.ylabel('Error')
    h = h/2