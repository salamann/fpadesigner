#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        Post implementation of 
# Purpose:
#
# Author:      shunyo
#
# Created:     31/03/2013
# Copyright:   (c) shunyo 2013
# Licence:     MIT
#-------------------------------------------------------------------------------

import numpy as np
import pylab as pl

def draw_spandirdata(yy,dL,clDist,circDist,ellipse,inducedAoa,planx,plany,dirname):
    pl.figure(figsize=(12,10))
    pl.axis("equal")
    pl.subplot(511)
    pl.plot(yy,dL,label = "dL")
    pl.legend()
    pl.subplot(512)
    pl.plot(yy,clDist,label = "dCL")
    pl.legend()
    pl.subplot(513)
    pl.plot(yy,circDist,label = "Gamma")
    pl.plot(yy,ellipse,label = "Ideal")
    pl.legend()
    pl.subplot(514)
    pl.plot(yy,inducedAoa,label = "Alpha_i")

    pl.legend()
    pl.subplot(515)
    pl.plot(planx,plany,label = "Planform")
    pl.legend()
    pl.xlabel("y [m]")
    pl.legend()
    pl.savefig(dirname + "/" +"span")

if __name__ == '__main__':
    pass
