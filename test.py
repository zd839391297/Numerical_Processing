# This is a test document
import NI
import LS
import LI
import ROOT


data_txt = "point_li.txt"
LI.li(data_txt)
data_txt = "point_ls.txt"
LS.ls(data_txt)
data_txt = 'point_ni.txt'
NI.ni(data_txt)
from numpy import poly1d
p = poly1d([1, -1, 0.2149, -0.01438])


def method(x):
    p = poly1d([1, -0.2149, 0.01438])
    return p(x)**(1/3.0)


ROOT.Dichotomy(0.7, 0.8, 0.00005, p)
ROOT.Iterative(0.7, 0.8, 0.00005, p, method)
ROOT.Wegstein(0.7, 0.8, 0.00005, p, method)
ROOT.Newton(0.7, 0.8, 0.00005, p)
ROOT.Chord(0.7, 0.8, 0.00005, p)
ROOT.sympy_m(p)


def q(x):
    return x**2-1


ROOT.Dichotomy(0.0, 2.0, 0.0005, q)
ROOT.Newton(0.0, 2.0, 0.0005, q)
