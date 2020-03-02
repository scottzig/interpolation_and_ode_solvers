#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 1 Problem 2

This code plots the function log(x) as well as
its second degree interpolating polynomial.
"""
import numpy as np
import matplotlib.pyplot as plt
import pylab

def f(x):
    return np.log(x)

def p(x):
    return (-3.33-5*np.log(2)+40*np.log(1.5)+41.67*np.log(1.6))*x**3 \
              +(17-20.5*np.log(2)-184*np.log(1.5)+187.5*np.log(1.6))*x**2 \
                +(-28.67+27.5*np.log(2)+272*np.log(1.5)-270.833*np.log(1.6))*x \
                  +(16-12*np.log(2)-128*np.log(1.5)+125*np.log(1.6))

oo = np.linspace(1,2)

plt.figure()
plt.title('f(x), p_2(x) and |f(x)-p_2(x)|')
y = f(oo)
pylab.plot(oo,y, label = 'f(x)')
plt.legend()
z = p(oo)
pylab.plot(oo,z, label='p_2(x)')
plt.legend()
z_2 = np.abs(f(oo)-p(oo))
pylab.plot(oo,z_2, label = '|f(x)-p_2(x)|')
plt.legend()