#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from FPA import Wing, Body, Tail
import weight
import io_fpa
import post_process_operation

if __name__ == '__main__':
    number_cell = 40
    velocity = 8.5
    temperature = 30
    # array_angle_of_attack = range(-3, 7)
    # aircraft_cg = 0.30
    aspect_ratio = 21.0
    surface_area = 22.0

    testWing = Wing('testplane.csv', number_cell, aspect_ratio, surface_area, optflag=0)
    testWing.temperature = temperature
    testWing.velocity = velocity
    # testWing.wingshape()
    testWing.calc_variedaoa(velocity, temperature, [0, 1, 2])
    ww = weight.calc_weight(testWing.span, "FX76-MP140")
    # testWing.calc_withconstWeight(ww, velocity, temperature)
    testWing.calc_planform()
    post_process_operation.draw_spandirdata(testWing.yy,
                                            testWing.dL,
                                            testWing.clDist,
                                            testWing.circDist,
                                            testWing.ellipse,
                                            testWing.inducedAoa,
                                            testWing.planx,
                                            testWing.plany,
                                            testWing.dirname)


    # zz = []
    # for j in [21, 23, 25, 27, 29, 31]: #aspect ratio
    #     for k in [20, 22, 24, 26]: #surface area
    #
    #         testWing = Wing('testplane.csv', number_cell, k, j, optical_wing_flag)
    #         testWing.temperature = temperature
    #         testWing.velocity = velocity
    #
    #         import random
    #         res = []
    #         for i in range(100000): #400000
    #             x1 = random.randrange(80, 101)
    #             x2 = random.randrange(70, 101)
    #             x3 = random.randrange(50, 101)
    #             x4 = random.randrange(45, 101)
    #             x5 = random.randrange(40, 101)
    #             x6 = random.randrange(40, 101)
    #             #print i
    #             x = [x1, x2, x3, x4, x5, x6]
    #             if x1 > x2 > x3 > x4 > x5 > x6:
    #                 res.append(testWing.opt_circ(x))
    #         res = np.array(res)
    #         resed = res[res[:, 8].argsort()]
    #         np.savetxt(testWing.dirname + "/" + "optresult.csv", np.array(resed), delimiter=',')
    #         print resed[len(resed)-1][8]
    #         print resed[len(resed)-2][8]
    #         print resed[len(resed)-3][8]
    #         testWing.opt_circ(resed[len(resed)-1][2:8])
    #         ww = weight.calc_weight(testWing.span, "FX76-MP140")
    #
    #         testWing.calc_withconstWeight(ww, velocity, temperature)
    #         testWing.calc_planform()
    #         post_process_operation.draw_spandirdata(testWing.yy,
    #                                                 testWing.dL,
    #                                                 testWing.clDist,
    #                                                 testWing.circDist,
    #                                                 testWing.ellipse,
    #                                                 testWing.inducedAoa,
    #                                                 testWing.planx,
    #                                                 testWing.plany,
    #                                                 testWing.dirname)
    #
    #
    #
    #         testBody = Body(velocity, temperature)
    #         testBody.fairdragcalc(1.136) #1.26がF-TECの値
    #         testBody.framedragcalc(0.0078)
    #         testTail = Tail(velocity, temperature)
    #         testTail.calc_htaildrag(2.0)
    #         testTail.calc_vtaildrag(1.5)
    #         power = io.gen_result(testWing, testBody, testTail)
    #         zz.append([j, k, power])
    #
    #         plt.clf()
    #         plt.figure(figsize=(8, 8))
    #         plt.plot(testWing.CD, testWing.CL, 'o')
    #
    #         testWing.calc_variedaoa(velocity, temperature, array_angle_of_attack)
    #         plt.plot(testWing.CDarray, testWing.CLarray)
    #         plt.plot(testWing.xmaxLDline, testWing.ymaxLDline)
    #         plt.xlim(xmin=0)
    #         plt.ylim(ymin=0)
    #         plt.xlabel("CD")
    #         plt.ylabel("CL")
    #         plt.legend(("Design Cruise Point", "Polar Curve", "maxL/D line"))
    #         plt.savefig(testWing.dirname + "/" + "PolarCurbe.png")
    #
    #         plt.clf()
    #         plt.plot(array_angle_of_attack, testWing.CDarray)
    #         plt.savefig(testWing.dirname + "/" + "alpha-CD.png")
    #
    #         plt.clf()
    #         plt.plot(array_angle_of_attack, testWing.CLarray)
    #         plt.savefig(testWing.dirname + "/" + "alpha-CL.png")
    #
    # print zz
    # zz = np.array(zz).transpose()
    # CS = plt.contourf(zz[0], zz[1], zz[2])
    # plt.clabel(CS, inline=1, fontsize=10)
    # plt.title('Required power by AR and S')
    # plt.savefig("contour.png")
