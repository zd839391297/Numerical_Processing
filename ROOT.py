# This program aims to Find the root of the linear equation

import numpy as np
import sympy as sp
# Dichotomy


def Dichotomy(a, b, delta, p):
    i = np.array([a, b])
    while abs(i[0]-i[1]) > delta:
        temp = (i[0]+i[1])/2
        if p(temp)*p(i[0]) > 0:
            i[0] = temp
        else:
            i[1] = temp
        # print(i)
    print("the result of Dichotomy is", i[0])

# Iterative method


def Iterative(a, b, delta, p, method):
    i = np.array([a, b])
    while abs(i[1]-i[0]) > delta:
        i[0] = method(i[1])
        # print(i)
        i[0], i[1] = i[1], i[0]
    print("the result of Iterative is", i[1])

# Wegstein


def Wegstein(a, b, delta, p, method):
    i = np.array([a, b, 0])

    def method_we(x1, x2):
        return (x2*method(x1)-x1*method(x2))/((x2-x1)-(method(x2)-method(x1)))

    while abs(i[0]-i[1]) > delta:
        i[2] = method_we(i[0], i[1])
        # print(i)
        i[0], i[1], i[2] = i[1], i[2], i[0]
    print("the result of Wegstein method is", i[1])

# Newton


def Newton(a, b, delta, p):
    i = np.array([a, b])

    def method_no(x):
        m = sp.symbols("m")
        f = p(m)
        return x-p(x)/(f.diff().evalf(subs={m: x}))

    while abs(i[0]-i[1]) > delta:
        i[0] = method_no(i[1])
        # print(i)
        i[0], i[1] = i[1], i[0]
    print("the result of Newton method is", i[1])

# Chord cut


def Chord(a, b, delta, p):
    i = np.array([a, b, 0])

    def method_ch(x1, x2):
        return x2-p(x2)*(x2-x1)/(p(x2)-p(x1))

    while abs(i[0]-i[1]) > delta:
        i[2] = method_ch(i[0], i[1])
        # print(i)
        i[0], i[1], i[2] = i[1], i[2], i[0]
    print("the result of Chord cut method is", i[1])

# sympy


def sympy_m(p):
    x = sp.symbols("x")
    print("the result of sympy is", sp.solve(p(x)))


if __name__ == "__main__":
    # import the equation
    def p(x):
        return (x**3-x**2+0.2149*x-0.01438)

    # input the iterative method

    def method(x):
        return (x**2-0.2149*x+0.01438)**(1/3.0)

    # input the Error limit
    delta = 0.00005

    # search initial numbers
    a = 0.0
    b = 0.1
    step = 0.1
    while p(a)*p(b) > 0:
        a += 0.2
        a, b = b, a
        # print(p(b))
    Dichotomy(a, b, delta, p)
    Iterative(a, b, delta, p, method)
    Wegstein(a, b, delta, p, method)
    Newton(a, b, delta, p)
    Chord(a, b, delta, p)
    sympy_m(p)
