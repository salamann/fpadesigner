#!/usr/bin/env python
# -*- coding: utf-8 -*-

from FPA import Wing
import weight
import post_process_operation

if __name__ == '__main__':
    number_cell = 40
    velocity = 8.5
    temperature = 30
    aspect_ratio = 21.0
    surface_area = 22.0

    testWing = Wing('testplane.csv', number_cell, aspect_ratio, surface_area, optflag=0)
    testWing.temperature = temperature
    testWing.velocity = velocity
    testWing.calc_variedaoa(velocity, temperature, [0, 1, 2])
    testWing.calc_planform()
    ww = weight.calc_weight(testWing.span, "FX76-MP140")
    post_process_operation.draw_spandirdata(testWing.yy,
                                            testWing.dL,
                                            testWing.clDist,
                                            testWing.circDist,
                                            testWing.ellipse,
                                            testWing.inducedAoa,
                                            testWing.planx,
                                            testWing.plany,
                                            testWing.dirname)

