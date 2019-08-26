import numpy as np
import sympy as sp

g1 = 3400  #处理湿物料量
#g1=input()
x_ = 0  #平衡含水量
rous = 1730  #颗粒密度
xc = 0.013  #临界含水量
roub = 800  #堆积密度
dm = 0.14  #颗粒平均直径
cs = 1.47  #干物料比热容
theta1 = 30  #进口温度
t0 = 24  #进预热器温度
#t0=input()
t1 = 80  #进干燥器温度
h0 = 0.024  #初始湿度
#h0=imput()
x1 = 0.04  #初始干基含水量
x2 = 0.0004  #结束干基含水量
yita = 0.3  #干燥装置热损失

omega1 = x1 / (1 + x1)  #初始湿基含水量
omega2 = x2 / (1 + x2)  #结束湿基含水量
print("初始湿基含水量ω1：", omega1)
print("结束湿基含水量ω2：", omega2)

g = g1 * (1 - omega1)  #固体物料质量流量
print("固体物料质量流量G：", g)

w = g * (x1 - x2)  #水分的蒸发量
print("水分的蒸发量W：", w)

h1 = h0  #湿空气离开预热器湿度
#h2 = 1  #湿空气离开干燥器湿度
#L = w / (h2 - h1)  #绝干空气流量
#print("绝干空气流量L：", L)

deltat = 35  #空气出口温度比出口处湿球温度高出的温度
#deltat=input()
tw1 = 37.5  #进口湿球温度，查湿度图得到
#tw1=input()
tw2 = tw1  #出口湿球温度
t2 = tw2 + deltat  #空气出口温度
print("空气出口温度t2：", t2)
rtw2 = 2412.9  #在tw2温度下水的汽化热，aspen算的
#rtw2=input()

x = sp.symbols("x")
a = sp.solve([(t2 - x) / (t2 - tw2) - (rtw2 * (x2 - x_) - cs * (t2 - tw2) *
                                       ((x2 - x_) /
                                        (xc - x_))**(rtw2 * (xc - x_) /
                                                     (cs * (t2 - tw2)))) /
              (rtw2 * (xc - x_) - cs * (t2 - tw2))], [x])
theta2 = a[x]  #物料离开干燥器的温度
print("物料离开干燥器的温度θ2：", theta2)

qd = 0  #干燥器补充热量
q1 = w * (2490 + 1.88 * t2) / 3600
print("Q1：", q1, "kW")
cm = cs + 4.187 * x2
q2 = g * cm * (theta2 - theta1) / 3600
print("Q2：", q2, "kW")
ql = yita * (q1 + q2)
print("Ql：", ql, "kW")
x = sp.symbols("x")
q3=x * (1.01 + 1.88 * h0) * (t2 - t0) / 3600
qp=x * (1.01 + 1.88 * h0) * (t1 - t0) / 3600
a = sp.solve([qp - q1 - q2 - q3 - ql], [x])
L = a[x]  #绝干空气流量
print("绝干空气流量L：", L)

x = sp.symbols("x")
temp1 = w / (x - h1)
a = sp.solve([L-temp1], [x])
h2 = a[x]  #湿空气离开干燥器湿度
print("湿空气离开干燥器湿度H2：", h2)

