import math
import sympy as sp
from scipy.optimize import root, fsolve
from PIL import Image

print("请输入设计任务：")
g1 = float(input("输入处理湿物料量（kg/h）<ep：2400>:"))  #处理湿物料量
h0 = float(input("输入空气初始湿度<ep:0.024>:"))  #空气初始湿度
t0 = float(input("输入空气进预热器温度（℃）<ep:24>:"))  #空气进预热器温度
x_ = 0  #平衡含水量
rous = 1730  #颗粒密度
xc = 0.013  #临界含水量
roub = 800  #堆积密度
dm = 0.00014  #颗粒平均直径
cs = 1.47  #干物料比热容
theta1 = 30  #物料进口温度
t1 = 80  #空气进干燥器温度
x1 = 0.04  #初始干基含水量
x2 = 0.0004  #结束干基含水量
yita = 0.15  #干燥装置热损失
z0 = 0.15  #静止床层厚度

print("2.1 物料衡算")
omega1 = x1 / (1 + x1)  #初始湿基含水量
omega2 = x2 / (1 + x2)  #结束湿基含水量
print("初始湿基含水量ω1：", omega1)
print("结束湿基含水量ω2：", omega2)
G = g1 * (1 - omega1)  #固体物料质量流量
print("固体物料质量流量G：", G, "kg/h")
w = G * (x1 - x2)  #水分的蒸发量
print("水分的蒸发量W：", w, "kg/h")
h1 = h0  #湿空气离开预热器湿度

print("2.2 空气和物料出口温度的确定")
deltat = float(
    input("输入空气出口温度比出口处湿球温度高出的温度,在20-50℃之间<ep:20>:"))  #空气出口温度比出口处湿球温度高出的温度
print("根据右图读出进口湿球温度，不会读的填35:")
img = Image.open('h_i.jpg')
img.show()
tw1 = float(input("输入进口湿球温度<ep:34.2>:"))  #进口湿球温度,aspen算的
tw2 = tw1  #出口湿球温度
t2 = tw2 + deltat  #空气出口温度
print("空气出口温度t2：", t2)
rtw2 = float(input("输入在tw2温度下水的汽化热<ep:2420.48>:"))  #在tw2温度下水的汽化热，aspen算的
x = sp.symbols("x")
temp1 = (rtw2 * (x2 - x_) - cs * (t2 - tw2) *
         ((x2 - x_) / (xc - x_))**(rtw2 * (xc - x_) /
                                   (cs *
                                    (t2 - tw2)))) / (rtw2 * (xc - x_) - cs *
                                                     (t2 - tw2))
a = sp.solve([(t2 - x) / (t2 - tw2) - temp1], [x])
theta2 = a[x]  #物料离开干燥器的温度
print("物料离开干燥器的温度θ2：", theta2)

print("3.3 干燥器的热量衡算")
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

print("3.4 预热器的热负荷和加热蒸汽消耗量")
qp = L * (1.01 + 1.88 * h0) * (t1 - t0) / 3600
print("预热器的热负荷Qp：", qp)
ts = float(input("输入加热蒸汽温度<ep:125>:"))  #加热蒸汽温度
gamma = float(input("输入加热蒸汽冷凝潜热<ep:2190.21>:"))  #加热蒸汽冷凝潜热
yita2 = float(input("输入预热器热损失比<ep:0.15>:"))  #预热器热损失
wh = qp * 3600 / (gamma * (1 - yita2))  #蒸汽消耗量
print("蒸汽消耗量Wh：", wh)
yitaz = q1 / qp  #干燥系统热效率
print("干燥系统热效率η：", yitaz)
print("请检查热效率是否大于30%，若不是，输入0退出，重新取空气出口温度比出口处湿球温度高出的温度或者重新读温焓图，若是，输入1继续")
temp = int(input("请输入0/1："))
if (temp == 0):
    exit()

print("3.1 流化速度的确定")
rou = 0.972166281370025  #80度空气的密度
miu = 0.0000209013074229616  #80度空气的黏度
lamda = 0.0295522535656798  #80度空气的导热系数
g = 9.81  #重力加速度
ar = (dm**3) * (rous - rou) * rou * g / (miu**2)  #阿基米德准数
print("80℃下空气的密度为0.9722kg/m3,黏度为2.0901e-5Pa*s,导热系数为0.02955W/m*℃")
print("重力加速度为9.81，临界流化点取0.4")
print("阿基米德准数Ar：", ar)
emf = 0.4  #球形颗粒床层临界流化点
print("根据右图读出临界流化李森科准数，不会读的填0.000002:")
img = Image.open('ly_ar.png')
img.show()
lymf = float(input("输入临界流化李森科准数<ep:0.000002>:"))  #临界流化李森科准数,由emf与Ar查图得到
umf = (lymf * miu * rous * g / (rou**2))**(1 / 3)  #临界流化速度
print("临界流化速度umf：", umf, "m/s")
e = 1  #球形颗粒带出时床层孔隙率
ly = 0.55  #带出李森科指数，由e与Ar查图得到
ut = (ly * miu * rous * g / (rou**2))**(1 / 3)  #带出速度
print("带出速度ut：", ut, "m/s")
temp=float(input("输入操作流化点，在0.55-0.75之间<ep:0.7>:"))#操作流化系数，0.55-0.75之间
u = temp * ut  #操作流化速度
print("操作流化速度u：", u, "m/s")

print("3.2 流化床层底面积的计算")
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
print("根据右图读出αa校正因子，不会读的填0.11")
img = Image.open('c.png')
img.show()
c = float(input("输入αa校正因子<ep:0.112>:"))  #αa校正因子，由dm查图得到
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

print("3.3 干燥器的宽度与长度")
print("根据总底面积估计长与宽，长应小于2.5-2.7,宽应小于2，若面积太小，可酌情增长宽")
length = float(input("输入长<ep:2.6>:"))  #长，估计值，理论小于2.5-2.7
width = float(input("输入宽<ep:2.5>:"))  #宽，估计值，理论小于2，超过需要设置特殊的物料散布装置
A = length * width  #流化床层实际总底面积
print("流化床层实际总底面积A：", A)

print("3.4 停留时间")
G2 = G * (1 + x2)
tao = z0 * A * roub * 60 / G2  #物料在床层的停留时间
print("物料在床层的停留时间t：", tao, 'min')

print("3.5 干燥器高度")
e = ((18 * re + 0.36 * re * re) / ar)**0.21  #床层孔隙率
print("床层孔隙率e：", e)
z = z0 * (1 - emf) / (1 - e)  #浓相区高度
print("浓相区高度Z1：", z)

de = 2 * (width * length / 4) / (width + length / 4)  #当量直径
print("当量直径De：", de)
print("根据右图读出Z2/De的值，不会读的填2.5")
img = Image.open('z2_de.png')
img.show()
k = float(input("输入Z2/De的值<ep:2.3>:"))  #z2/de的值，由u、De查图得到
z2 = k * de  #分离段高度
print("分离段高度Z2：", z2)
zz=z+z2#总高度
print("总高度Z：", zz)

print("3.5 干燥器结构设计")
print("3.5.1 布气装置--分布板")
print("分布板压降为床层压降的15%")
e0 = 1 - roub / rous  #静态床层孔隙率
deltap = 0.15  #分布板压降占床层压降的比重
deltapd = deltap * z0 * (1 - e0) * (rous - rou) * g  #气体通过分布板的压力降
print("气体通过分布板的压力降ΔPd：", deltapd)
kezan = float(input("输入分布板阻力系数，在1.1-2.5之间<ep:2>:"))  #分布板的阻力系数
u0 = (2 * deltapd / (kezan * rou))**0.5  #气体通过筛板的速度
print("气体通过筛板的速度u0：", u0)
vs = L * (0.772 + 1.244 * h0) * (t1 + 273) / (273 * 3600)  #干燥介质热空气的体积流量
print("干燥介质热空气的体积流量Vs：", vs)
d0 = float(input("输入筛孔直径<ep:0.0015>:"))  #筛孔直径，m
n0 = vs / (math.pi * d0 * d0 * u0 / 4)  #总筛孔数
print("总筛孔数n0：", n0)
phi = math.pi * d0 * d0 * n0 / (4 * A)  #实际开孔率
print("实际开孔率φ：", phi)
t = 0.952 * d0 / (math.sqrt(phi))  #孔心距
print("孔心距t：", t, 'm')
print("请检查开孔率是否在3-13%，若不是，输入0退出，重新取筛孔直径或阻力系数，若是，输入1继续")
temp = int(input("请输入0/1："))
if (temp == 0):
    exit()

print("3.5.3 物料出口溢流堰高")
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

print("4.1 送风机与排风机")
v1 = L * (0.772 + 1.244 * h0) * (273 + t0) / 273  #送风机风量
print("送风机风量v1：", v1, "m3/h")
v2 = L * (0.772 + 1.244 * h0) * (273 + t2) / 273  #送风机风量
print("送风机风量v2：", v2, "m3/h")
