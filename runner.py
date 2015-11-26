#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from FPA import Wing, Body, Tail
import weight
import iomod
import postmod

if __name__ == '__main__':
    # 4分割で航空力学の基礎第2版 p.147の計算になる
    # wing('kokurikgaku_p.147.csv',4)
    ##    testWing = wing('testplane.csv',40)
    ##    testWing.calc_reynolds(9, 30)
    ncell = 40
    velocity = 8.5
    temperature = 30
    aoaarray = range(-3, 7)
    aircraftcg = 0.30
    optwing = 1
    #directshape = [[0.0, 100.0], [40.0, 100.0], [70.0, 75.0], [100.0, 50.0]]

    zz = []
    for j in [21, 23, 25, 27, 29, 31]: #aspect ratio
        for k in [20, 22, 24, 26]: #surface area
        
            testWing = Wing('testplane.csv', ncell, k, j, optwing)
            testWing.temperature = temperature
            testWing.velocity = velocity
        #x = [70.,50.]
        #a1 = testWing.opt_circ(x)
        #print a1
        #
        #res = minimize(testWing.opt_circ,(90.,80.,70.,60.,50.,40.,30.),method = "SLSQP")
        #res = minimize(testWing.opt_circ,(86.,89.2,88.7,80.1,58.4,38.6,28.3),method = "L-BFGS-B")
        #res = basinhopping(testWing.opt_circ,(86.,89.2,88.7.,80.1.,58.4.,38.6.,28.3.),niter=1,stepsize=0.1 )

            import random
            res = []
            for i in range(100000): #400000
                x1 = random.randrange(80, 101)
                x2 = random.randrange(70, 101)
                x3 = random.randrange(50, 101)
                x4 = random.randrange(45, 101)
                x5 = random.randrange(40, 101)
                x6 = random.randrange(40, 101)
                #print i
                x = [x1, x2, x3, x4, x5, x6]
                if x1 > x2 and x2 > x3 and x3 > x4 and x4 > x5 and x5 > x6:
                    #print "hit", i
                    res.append(testWing.opt_circ(x))
            #np.savetxt('test.csv', np.array(res), delimiter=',')
            res = np.array(res)
            resed = res[res[:,8].argsort()]
            np.savetxt(testWing.dirname + "/" + "optresult.csv", np.array(resed), delimiter=',')
            print resed[len(resed)-1][8]
            print resed[len(resed)-2][8]
            print resed[len(resed)-3][8]
            testWing.opt_circ(resed[len(resed)-1][2:8])
            ww = weight.calc_weight(testWing.span,"FX76-MP140")

            testWing.calc_withconstWeight(ww, velocity, temperature) #96.007
            #Constrained Optimization can be calculated by scipy.minimize.
            #func = self.ellipse
            #
            #
            testWing.calc_planform()
            postmod.draw_spandirdata(testWing.yy,
                                     testWing.dL,
                                     testWing.clDist,
                                     testWing.circDist,
                                     testWing.ellipse,
                                     testWing.inducedAoa,
                                     testWing.planx,
                                     testWing.plany,
                                     testWing.dirname)

            #print "total weight  ",ww


            testBody = Body(velocity, temperature)
            testBody.fairdragcalc(1.136) #1.26がF-TECの値
            testBody.framedragcalc(0.0078)
            testTail = Tail(velocity, temperature)
            testTail.calc_htaildrag(2.0)
            testTail.calc_vtaildrag(1.5)
            power = iomod.gen_result(testWing,testBody,testTail)
            zz.append([j, k, power])

            plt.clf()
            plt.figure(figsize=(8, 8))
            plt.plot(testWing.CD, testWing.CL, 'o')

            testWing.calc_variedaoa(velocity,temperature,aoaarray)
            plt.plot(testWing.CDarray, testWing.CLarray)
            plt.plot(testWing.xmaxLDline, testWing.ymaxLDline)
            plt.xlim(xmin=0)
            plt.ylim(ymin=0)
            plt.xlabel("CD")
            plt.ylabel("CL")
            plt.legend(("Design Cruise Point", "Polar Curve", "maxL/D line"))
            plt.savefig(testWing.dirname + "/" + "PolarCurbe.png")

            plt.clf()
            plt.plot(aoaarray, testWing.CDarray)
            plt.savefig(testWing.dirname + "/" + "alpha-CD.png")

            plt.clf()
            plt.plot(aoaarray, testWing.CLarray)
            plt.savefig(testWing.dirname + "/" + "alpha-CL.png")



    print zz
    zz = np.array(zz).transpose()
    CS = plt.contourf(zz[0], zz[1], zz[2])
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('Required power by AR and S')
    plt.savefig("contour.png")
