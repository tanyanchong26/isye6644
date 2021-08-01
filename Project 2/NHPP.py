# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 10:56:48 2021

@author: tanya
"""

from math import sin, log
from numpy.random import uniform
import matplotlib.pyplot as plt

def rate_func(t):
    return 1+sin(t/5)

T = [0]
lambda_vec = [rate_func(T[-1])]
lambda_max = 2


while T[-1]<40:
    U = uniform()
    V = uniform()
    t = T[-1]
    
    while V > rate_func(t)/lambda_max:
        U = uniform()
        V = uniform()
        t += -1/lambda_max*log(U)
    
    lambda_ = rate_func(t)
    T.append(t)
    lambda_vec.append(lambda_)

plt.scatter(T, lambda_vec)
plt.show()


rate_func(t)/2