#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 4 Problem 2d

This code approximates the derivative of f(x)=e^x at x_0=0
using a three-point and five-point finite difference operator.
"""
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.exp(x)

# three point stencil
err = 100
h = 2
while h > 1/64:
    der = (-f(-h)+f(h))/(2*h)
    print(der)
    err = np.abs(der - 1)
    plt.loglog(h,err,'ko')
    plt.loglog(h,h**2,'ro')
    plt.title('Three point stencil in black, h^2 in red')
    plt.xlabel('h')
    plt.ylabel('Error in derivative approx.')
    h = h/2
    
# five point stencil
plt.figure()
err = 100
h = 2
while h > 1/64:
    der = (f(-2*h)-8*f(-h)+8*f(h)-f(2*h))/(12*h)
    print(der)
    err = np.abs(der - 1)
    plt.loglog(h,err,'ko')
    plt.loglog(h,h**4,'ro')
    plt.title('Five point stencil in black, h^4 in red')
    plt.xlabel('h')
    plt.ylabel('Error in derivative approx.')
    h = h/2
    