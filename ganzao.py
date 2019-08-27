import math
import sympy as sp
from scipy.optimize import root, fsolve

g1 = 3400  #处理湿物料量
#g1=input()
x_ = 0  #平衡含水量
rous = 1730  #颗粒密度
xc = 0.013  #临界含水量
roub = 800  #堆积密度
dm = 0.00014  #颗粒平均直径
cs = 1.47  #干物料比热容
theta1 = 30  #物料进口温度
t0 = 24  #空气进预热器温度
#t0=input()
t1 = 80  #空气进干燥器温度
h0 = 0.024  #初始湿度
#h0=imput()
x1 = 0.04  #初始干基含水量
x2 = 0.0004  #结束干基含水量
yita = 0.15  #干燥装置热损失

omega1 = x1 / (1 + x1)  #初始湿基含水量
omega2 = x2 / (1 + x2)  #结束湿基含水量
print("初始湿基含水量ω1：", omega1)
print("结束湿基含水量ω2：", omega2)

G = g1 * (1 - omega1)  #固体物料质量流量
print("固体物料质量流量G：", G)

w = G * (x1 - x2)  #水分的蒸发量
print("水分的蒸发量W：", w)

h1 = h0  #湿空气离开预热器湿度

deltat = 20  #空气出口温度比出口处湿球温度高出的温度
#deltat=input()
tw1 = 34.2  #进口湿球温度，aspen计算得到
#tw1=input()
tw2 = tw1  #出口湿球温度
t2 = tw2 + deltat  #空气出口温度
print("空气出口温度t2：", t2)
rtw2 = 2420.48  #在tw2温度下水的汽化热，aspen算的
#rtw2=input()

x = sp.symbols("x")
temp1 = (rtw2 * (x2 - x_) - cs * (t2 - tw2) *
         ((x2 - x_) / (xc - x_))**(rtw2 * (xc - x_) /
                                   (cs *
                                    (t2 - tw2)))) / (rtw2 * (xc - x_) - cs *
                                                     (t2 - tw2))
a = sp.solve([(t2 - x) / (t2 - tw2) - temp1], [x])
theta2 = a[x]  #物料离开干燥器的温度
print("物料离开干燥器的温度θ2：", theta2)

qd = 0  #干燥器补充热量
q1 = w * (2490 + 1.88 * t2) / 3600
print("Q1：", q1, "kW")
cm = cs + 4.187 * x2
q2 = G * cm * (theta2 - theta1) / 3600
print("Q2：", q2, "kW")
ql = yita * (q1 + q2)
print("Ql：", ql, "kW")
x = sp.symbols("x")
q3 = x * (1.01 + 1.88 * h0) * (t2 - t0) / 3600
qp = x * (1.01 + 1.88 * h0) * (t1 - t0) / 3600
a = sp.solve([qp - q1 - q2 - q3 - ql], [x])
L = a[x]  #绝干空气流量
print("绝干空气流量L：", L)

x = sp.symbols("x")
temp1 = w / (x - h1)
a = sp.solve([L - temp1], [x])
h2 = a[x]  #湿空气离开干燥器湿度
print("湿空气离开干燥器湿度H2：", h2)

qp = L * (1.01 + 1.88 * h0) * (t1 - t0) / 3600
print("预热器的热负荷Qp：", qp)

ts = 125  #加热蒸汽温度
gamma = 2190.21  #加热蒸汽冷凝潜热
yita2 = 0.15  #预热器热损失
wh = qp * 3600 / (gamma * (1 - yita2))  #蒸汽消耗量
print("蒸汽消耗量Wh：", wh)

yitaz = q1 / qp  #干燥系统热效率
print("干燥系统热效率η：", yitaz)

rou = 0.972166281370025  #80度空气的密度
miu = 0.0000209013074229616  #80度空气的黏度
lamda = 0.0295522535656798  #80度空气的导热系数
g = 9.81  #重力加速度
ar = (dm**3) * (rous - rou) * rou * g / (miu**2)  #阿基米德准数
print("阿基米德准数Ar：", ar)
emf = 0.4  #球形颗粒床层临界流化点
lymf = 2 / 1000000  #临界流化李森科准数,由emf与Ar查图得到
umf = (lymf * miu * rous * g / (rou**2))**(1 / 3)  #临界流化速度
print("临界流化速度umf：", umf, "m/s")
e = 1  #球形颗粒带出时床层孔隙率
ly = 0.55  #带出李森科指数，由e与Ar查图得到
ut = (ly * miu * rous * g / (rou**2))**(1 / 3)  #带出速度
print("带出速度ut：", ut, "m/s")
u = 0.7 * ut  #操作流化速度
print("操作流化速度u：", u, "m/s")

z0 = 0.15  #静止床层厚度
L_ = rou * u  #干空气的质量流速
print("干空气的质量流速L-：", L_)
a = 6 * roub / (rous * dm)  #静止床层的比表面积
print("静止床层的比表面积a：", a)
re = dm * u * rou / miu  #雷诺数
print("雷诺数Re：", re)
alpha = 4 * (10**-3) * lamda * (re**1.5) / dm  #流化床层的对流传热系数
print("流化床层的对流传热系数α：", alpha)
alphaa = alpha * a  #流化床层的体积传热系数
print("流化床层的体积传热系数αa：", alphaa)
c = 0.11  #αa校正因子，由dm查图得到
alphaa_ = c * alphaa  #校正后流化床层的体积传热系数
print("校正后流化床层的体积传热系数αa'：", alphaa_)
x = sp.symbols("x")
a = sp.solve([
    alphaa_ * z0 - (1.01 + 1.88 * h0) * L_ /
    (((1.01 + 1.88 * h0) * L_ * x * (t1 - tw1)) /
     (G * (x1 - x2) * rtw2 / 3600) - 1)
], [x])
print((1.01 + 1.88 * h0) * L_)  #分子系数
print(((1.01 + 1.88 * h0) * L_ * (t1 - tw1)) /
      (G * (x1 - x2) * rtw2 / 3600))  #分母系数
a1 = a[x]  #干燥第一阶段所需底面积
print("干燥第一阶段所需底面积A1：", a1)

cm2 = cs + 4.187 * x2  #干燥产品的比热容
print("干燥产品的比热容cm2：", cm2)
x = sp.symbols("x")
a = sp.solve([
    alphaa_ * z0 - (1.01 + 1.88 * h0) * L_ /
    (((1.01 + 1.88 * h0) * L_ * x) / (G * cm2 * math.log(
        (t1 - theta1) / (t1 - theta2)) / 3600) - 1)
], [x])
print((1.01 + 1.88 * h0) * L_ / (G * cm2 * math.log(
    (t1 - theta1) / (t1 - theta2)) / 3600))  #分母系数
a2 = a[x]  #物料升温阶段所需底面积
print("物料升温阶段所需底面积A2：", a2)

A0 = a1 + a2  #流化床层总底面积
print("流化床层总底面积A：", A0)

length = 2.6  #长，估计值，理论小于2.5-2.7
width = 2.5  #宽，估计值，理论小于2，超过需要设置特殊的物料散布装置
A = length * width  #流化床层实际总底面积
print("流化床层实际总底面积A：", A)

G2 = G * (1 + x2)
tao = z0 * A * roub * 60 / G2  #物料在床层的停留时间
print("物料在床层的停留时间t：", tao, 'min')

e = ((18 * re + 0.36 * re * re) / ar)**0.21  #床层孔隙率
print("床层孔隙率e：", e)
z = z0 * (1 - emf) / (1 - e)  #浓相区高度
print("浓相区高度Z1：", z)

de = 2 * (width * length / 4) / (width + length / 4)  #当量直径
print("当量直径De：", de)
k = 2.3  #z2/de的值，由u、De查图得到
z2 = k * de  #分离段高度
print("分离段高度Z2：", z2)

e0 = 1 - roub / rous  #静态床层孔隙率
deltap = 0.15  #分布板压降占床层压降的比重
deltapd = deltap * z0 * (1 - e0) * (rous - rou) * g  #气体通过分布板的压力降
print("气体通过分布板的压力降ΔPd：", deltapd)
kezan = 2  #分布板的阻力系数
u0 = (2 * deltapd / (kezan * rou))**0.5  #气体通过筛板的速度
print("气体通过筛板的速度u0：", u0)
vs = L * (0.772 + 1.244 * h0) * (t1 + 273) / (273 * 3600)  #干燥介质热空气的体积流量
print("干燥介质热空气的体积流量Vs：", vs)
d0 = 0.0015  #筛孔直径，m
n0 = vs / (math.pi * d0 * d0 * u0 / 4)  #总筛孔数
print("总筛孔数n0：", n0)
phi = math.pi * d0 * d0 * n0 / (4 * A)  #实际开孔率
print("实际开孔率φ：", phi)
t = 0.952 * d0 / (math.sqrt(phi))  #孔心距
print("孔心距t：", t, 'm')

ret = dm * ut * rou / miu  #带出雷诺数
print("带出雷诺数Ret：", ret)
x = sp.symbols("x")
a = sp.solve([(x - 1) / (u - umf) - 25 / (ret**0.44)], [x])
ev = a[x]  #床层膨胀率
print("床层膨胀率Ev：", ev)


def f(x):
    x0 = float(x[0])
    return [
        2.14 * (z0 - x0 / ev) / (((1 / ev)**(1 / 3)) *
                                 (G / (3600 * width * roub))) - 18 +
        1.52 * math.log(ret / (5 * x0))
    ]


h = fsolve(f, [0.1])  #溢流堰高度
print("溢流堰高度h：", h[0])

v1=L*(0.772+1.244*h0)*(273+t0)/273#送风机风量
print("送风机风量v1：", v1,"m3/h")
v2=L*(0.772+1.244*h0)*(273+t2)/273#送风机风量
print("送风机风量v2：", v2,"m3/h")