#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os

import weight
from fpa import post_process_operation
from fpa.FPA import Wing


def main(wing, velocity, temperature, aoaarray):
    """
    csvfile, number of cell, design cruise speed, ambient temperature
    :param velocity:
    :param temperature:
    :param aoaarray:
    :return:
    """

    wing.set_temperature(temperature)
    wing.set_velocity(velocity)
    wing.calc_reynolds()
    wing.set_aerodynamcis_data()
    """機体の特性を出す"""
    CLarray = []
    CDarray = []
    for i in aoaarray:
        print "calculating aerodynamics at alpha = {} [deg]".format(str(i))
        wing.calc_lift_slope_and_zero_lift_array()
        wing.calc_CL_Cdi_CD(i)
        CLarray.append(wing.CL)
        CDarray.append(wing.CD)
    return np.array(CLarray), np.array(CDarray)

    # wing.CLarray = CLarray
    # wing.CDarray = CDarray

    # """maxL/Dの線を引くためのリスト生成"""
    # # j = 0
    # maxslope = None
    # for i in np.array(CLarray)/np.array(CDarray):
    #     if i == max(np.array(CLarray)/np.array(CDarray)):
    #         maxslope = i
    #     # j += 1
    # wing.maxslope = maxslope
    #
    # xlist = list(CDarray)
    # ylist = list(np.array(wing.CDarray)*wing.maxslope)
    # xlist.insert(0, 0)
    # ylist.insert(0, 0)
    # wing.xmaxLDline = xlist
    # wing.ymaxLDline = ylist


if __name__ == '__main__':
    number_cell = 40
    velocity = 8.5
    temperature = 30
    aspect_ratio = 21.0
    surface_area = 24.0
    aoa_array = np.arange(-5, 15, 0.5)

    testWing = Wing('testplane.csv', number_cell, aspect_ratio, surface_area, optflag=0)
    cl, cd = main(testWing, velocity, temperature, aoa_array)
    data_for_record = np.array([aoa_array, cl, cd]).transpose()
    if not os.path.isdir("./results/{}".format(testWing.dirname)):
        os.mkdir("./results/{}".format(testWing.dirname))
    np.savetxt("./results/{}/aerodynamics_data.csv".format(testWing.dirname), data_for_record,
               delimiter=",", header="alpha, CL, CD")

    # testWing.calc_variedaoa(velocity, temperature, range(0, 10))
    # testWing.calc_planform()
    # ww = weight.calc_weight(testWing.span, "FX76-MP140")
    # post_process_operation.draw_spandirdata(testWing.yy,
    #                                         testWing.dL,
    #                                         testWing.clDist,
    #                                         testWing.circDist,
    #                                         testWing.ellipse,
    #                                         testWing.inducedAoa,
    #                                         testWing.planx,
    #                                         testWing.plany,
    #                                         testWing.dirname)
    #
