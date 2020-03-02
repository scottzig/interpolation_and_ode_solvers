#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 4 Problem 1

This code approximates the integral of sin(1/x) on 0.05 to 1
using a Monte Carlo integration method.
"""
import numpy as np
import matplotlib.pyplot as plt

true_val = 0.50283962;

def f(x):
    return np.sin(1/x)

a = 0.05
b = 1

N = 2000 
for i in range(1,N):
    sum = 0
    for j in range(1,i):
        x = np.random.uniform(a,b)
        sum = sum+f(x)
    sum = (1/i)*sum
    err = np.abs(sum-true_val)
    plt.semilogy(i,err,'.')
    plt.title('Monte Carlo Integration')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Error')
    