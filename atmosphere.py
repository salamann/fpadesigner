#!/usr/bin/env python
# -*- coding: utf-8 -*-

class Air(object):
    def __init__(self, temperature):
        self.kinetic_viscosty = 1.34 * 10 ** - 5. + 9.31477 * 10 ** - 8. * temperature
        self.airdensity = 1.28912 - 0.004122391 * temperature

    def get_dynamic_pressure(self, airspeed):
        self.dynamic_pressure = 0.5 * self.airdensity * airspeed ** 2.0