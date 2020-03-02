#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 4 Problem 3a

This code approximates the minimum of f(x) = xe^x-1 by using
Newton's method on the function g(x)=f'(x).
"""
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x*np.exp(x)-1

def g(x):
    return x*np.exp(x)+np.exp(x)

def der_g(x):
    return x*np.exp(x)+np.exp(x)+np.exp(x)

x_0 = -0.5
err = 100
x = x_0
while err > 10**(-6):
    x_new = x - g(x)/der_g(x)
    err = np.abs(x-x_new)
    x = x_new

txt = str(x)
print('The minimum is at point ' + txt)

