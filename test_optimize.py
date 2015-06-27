#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pylab as pl
from scipy.linalg import solve
from scipy.optimize import minimize, rosen, rosen_der

def tset(x):
    return (x[0] - 1)**2. + (x[1] - 2.5)**2.

def const1(x):
    return x[0] -2.*x[1]+2.

def const2(x):
    return -x[0] -2. * x[1] +6.

def const3(x):
    return -x[0] + 2. * x[0] +2.


def test2(x):
    return x[0]**2. +x[1]**2.

def const4(x):
    return -x[0]+x[1]-1

if __name__ == '__main__':
    #fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2
    #fun = tset(x)
    #cons = ({'type': 'ineq', 'fun': const1},
     #   {'type': 'ineq', 'fun': const2},
      #  {'type': 'ineq', 'fun': const3})
    cons = ({'type': 'eq', 'fun': const4})

    bnds = ((0, None), (0, None))
    import datetime
    print datetime.datetime.today()
    res = minimize(test2, (20, -20), method='SLSQP', bounds=bnds,constraints=cons)
    print res
