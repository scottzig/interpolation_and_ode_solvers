#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 3 Problem 1

This code approximates the integral of sin(1/x) on 0.05 to 1
using 5 different methods: midpoint rule, trapezoidal rule,
Simpson rule, 2-point Gauss rule, and 3-point Gauss rule
"""
import numpy as np
import matplotlib.pyplot as plt

true_val = 0.50283962;

def f(x):
    return np.sin(1/x)

def gauss(x,x_1,x_2):
    return ((x_2-x_1)/2)*np.sin(1/(((x_2-x_1)*x+x_2+x_1)/2))

# Midpoint Rule
    
plt.figure()
K = 1
err = 1
while err > 10**(-6):
    sum = 0
    for j in range(0,K):
        x_1 = 0.05+(1-0.05)*j/K
        x_2 = 0.05+(1-0.05)*(j+1)/K
        x_cen = (x_1+x_2)/2
        sum = sum + (x_2-x_1)*f(x_cen)
    err = np.abs(true_val - sum)
    plt.plot(K,np.log10(err),'.')
    plt.title('Midpoint Rule')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Error (log scale)')
    K = K+1
    
txt = str(K+1)
print('Midpoint function evaluations required = ' + txt)
    
# Trapezoidal Rule

plt.figure()  
K = 1
err = 1
while err > 10**(-6):
    sum = 0
    for j in range(0,K):
        x_1 = 0.05+(1-0.05)*j/K
        x_2 = 0.05+(1-0.05)*(j+1)/K
        sum = sum + (x_2-x_1)*(f(x_1)+f(x_2))/2
    err = np.abs(true_val - sum)
    plt.plot(K+1,np.log10(err),'.')
    plt.title('Trapezoidal Rule')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Error (log scale)')
    K = K+1
    
txt = str(K+2)
print('Trapezoidal function evaluations required = ' + txt)
    
# Simpson's Rule

plt.figure()  
K = 1
err = 1
while err > 10**(-6):
    sum = 0
    for j in range(0,K):
        x_1 = 0.05+(1-0.05)*j/K
        x_2 = 0.05+(1-0.05)*(j+1)/K
        x_cen = (x_1+x_2)/2
        sum = sum + (x_2-x_1)*(f(x_1)+4*f(x_cen)+f(x_2))/6
    err = np.abs(true_val - sum)
    plt.plot(2*K+1,np.log10(err),'.')
    plt.title('Simpsons Rule')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Error (log scale)')
    K = K+1
            
txt = str(2*K+2)
print('Simpsons function evaluations required = ' + txt)

# 2-point Gauss rule

plt.figure()  
K = 1
err = 1
while err > 10**(-6):
    sum = 0
    for j in range(0,K):
        x_1 = 0.05+(1-0.05)*j/K
        x_2 = 0.05+(1-0.05)*(j+1)/K
        sum = sum + (gauss(-1/np.sqrt(3),x_1,x_2)+gauss(1/np.sqrt(3),x_1,x_2))
    err = np.abs(true_val - sum)
    plt.plot(2*K,np.log10(err),'.')
    plt.title('2 Point Gauss Rule')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Error (log scale)')
    K = K+1
            
txt = str(2*K+1)
print('2 Point Gauss Rule function evaluations required = ' + txt)

# 3-point Gauss rule

plt.figure()  
K = 1
err = 1
while err > 10**(-6):
    sum = 0
    for j in range(0,K):
        x_1 = 0.05+(1-0.05)*j/K
        x_2 = 0.05+(1-0.05)*(j+1)/K
        sum = sum + ((8/9)*gauss(0,x_1,x_2)+(5/9)*gauss(np.sqrt(3/5),x_1,x_2)+(5/9)*gauss(-np.sqrt(3/5),x_1,x_2))
    err = np.abs(true_val - sum)
    plt.plot(3*K,np.log10(err),'.')
    plt.title('3 Point Gauss Rule')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Error (log scale)')
    K = K+1
            
txt = str(3*K+1)
print('3 Point Gauss Rule function evaluations required = ' + txt) 