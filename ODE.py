# code = utf-8
# This program aims to using the fourth-order Runge-Kutta method to solve ordinary differential equations
import numpy as np
import math
# from scipy.integrate import odeint

# Input the fuction


def f(x, y):
    return 8-3*y

# Set iterative method


def yn(x, y, h):
    K1 = h*f(x, y)
    print(K1)
    K2 = h*f(x+h/2, y+K1/2)
    print(K2)
    K3 = h*f(x+h/2, y+K2/2)
    print(K3)
    K4 = h*f(x+h, y+K3)
    print(K4)
    return y+(K1+2*K2+2*K3+K4)/6


if __name__ == "__main__":

    # Set the step size
    h = 0.2
    # Input the initial number(x,y)
    num = np.array([0.0, 2.0])
    # Set the error limit
    error = 0.00005
    # Set the end x
    x_end = 0.4
    # Solve
    while num[0] < 0.4:
        num[1] = yn(num[0], num[1], h)
        num[0] += h
        print(num)
