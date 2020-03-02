#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 2 Problem 3a

This code runs Newton's method to find the root of
f(x) = xe^x - 1 to six digits of accuracy starting at
x_0 = 0.5.s
"""
import numpy as np

def f(x):
    return x*np.exp(x)-1

def der_f(x):
    return x*np.exp(x)+np.exp(x)

x_0 = 0.5
err = 100
while err > .000001:
    x_new = x_0 - f(x_0)/der_f(x_0)
    err = np.abs(x_new-x_0)
    x_0 = x_new
    
print(x_new)

