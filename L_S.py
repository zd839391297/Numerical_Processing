# This is a simple example of using least squares to fit data
import numpy as np
import math

point = np.loadtxt("point_ls.txt")
point_shape = point.shape
if point_shape[1] > 2 or point_shape[1] < 2:
    print("Sorry,updating")
else:
    point_num = point_shape[0]
    x_sum = sum(point[:, 0])
    y_sum = sum(point[:, 1])
    i = 0
    xy_num = np.zeros(shape=(8, 1))
    while i < point_num:
        xy_num[i][0] = point[i][0]*point[i][1]
        i = i+1
    xy_sum = sum(xy_num)
    i = 0
    xx_num = np.zeros(shape=(8, 1))
    while i < point_num:
        xx_num[i][0] = point[i][0]*point[i][0]
        i = i+1
    xx_sum = sum(xx_num)
    b = (xy_sum-(x_sum*y_sum)/point_num)/(xx_sum-(x_sum*x_sum)/point_num)
    x_mean = np.mean(point[:, 0])
    y_mean = np.mean(point[:, 1])
    a = y_mean-b*x_mean
    print("y=", b, "x+", a)
    b_s = math.e**b
    print(b_s)
