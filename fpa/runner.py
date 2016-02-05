#!/usr/bin/env python
# -*- coding: utf-8 -*-

import weight
from fpa import post_process_operation
from fpa.FPA import Wing

if __name__ == '__main__':
    number_cell = 40
    velocity = 8.5
    temperature = 30
    aspect_ratio = 21.0
    surface_area = 24.0

    testWing = Wing('testplane.csv', number_cell, aspect_ratio, surface_area, optflag=0)
    testWing.calc_variedaoa(velocity, temperature, range(0, 10))
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

