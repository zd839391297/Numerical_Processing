# This is an example for Lagrangian interpolation
import numpy as np
from sympy import *


point = np.loadtxt("point_file.txt")
x = Symbol("x")
# x = 0.872
# print(point)
point_shape = point.shape
if point_shape[1] > 2 or point_shape[1] < 2:
    print("Sorry, updating")
else:
    point_num = point_shape[0]
    # print(point_num)

    def LI(x, j):
        li = 1
        i = 1
        while i <= point_num:
            if i != j:
                # print("i=", i)
                li *= (x-point[i-1][0])/(point[j-1][0]-point[i-1][0])
                # print("li=", li)
                i = i+1
            else:
                i = i+1
        return li

    def LN(x):
        ln = 0
        j = 1
        while j <= point_num:
            # print("j=", j)
            ln += LI(x, j)*point[j-1][1]
            # print("ln=", ln)
            j = j+1
        return expand(ln)
    print(LN(x))