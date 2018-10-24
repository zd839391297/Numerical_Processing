# This program aims to do the numberical integration progress
# Input: a n*2 data frame :
# The first line is the independent variable(x),
#   the second line is the dependent variable(y)
# The method is the complex simpson
# However, while the number of independent variables is even,
#   this program will find a suitable position to integration one segment
#       by the trapezoidal formula
# The independent variable number must be evenly arranged


# pragram here
import numpy as np

point = np.loadtxt("point_ni.txt")
# print (point.shape[0])
point_num = point.shape[1]
x_in = point[0][1]-point[0][0]


def complex_simpson(x, x_in):
    x_num = x.shape[0]
    s_pr = x[0]+x[x_num-1]
    for i in range(1, x_num-1):
        if i % 2 != 0:
            s_pr += 4*x[i]
        else:
            s_pr += 2*x[i]
    s = x_in*s_pr/6
    return s


def sort(a):
    for i in range(0, a.shape[1]):
        for j in range(0, a.shape[1]):
            if a[1][i] < a[1][j]:
                a[:, [i, j]] = a[:, [j, i]]


def trapezoidal(fa, fb, x_in):
    s = (fa+fb)*x_in/2
    return s


if point.shape[1] % 2 == 0:
    # print(1)
    point_pg = np.zeros(shape=[2, point_num-1])
    for i in range(0, point_num-1):
        point_pg[0][i] = i
    for i in range(0, point_num-1):
        point_pg[1][i] = abs(point[1][i+1]-point[1][i])
    # position = np.where(point_pg == point_pg.min)
    sort(point_pg)
    j = 0
    while point_pg[0][j] % 2 != 0:
        j += 1
    i = int(point_pg[0][j])
    point_l = point[:, 0:i+1]
    point_r = point[:, (i+2):point_num]
    result = 0
    if point_l.size != 0:
        result += complex_simpson(point_l[1], x_in)
    # print(complex_simpson(point_l[1], x_in))
    if point_r.size != 0:
        result += complex_simpson(point_r[1], x_in)
    result += trapezoidal(point[1][i], point[1][i+1], x_in)
else:
    result = complex_simpson(point[1], x_in)
print("The numberical integration result is: ", result)
# print((0+4*0.5687+2*0.7909+4*0.5743+2*0.1350-4*0.1852-2*0.1802+4*0.0811+0.2917)/6)
# print((0.2917+0.3031)/2)
