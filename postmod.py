#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt


def draw_spandirdata(yy, dL, clDist, circDist, ellipse, inducedAoa, planx, plany, dirname):
    plt.figure(figsize=(12, 10))
    plt.axis("equal")
    plt.subplot(511)
    plt.plot(yy, dL, label="dL")
    plt.legend()
    plt.subplot(512)
    plt.plot(yy, clDist, label="dCL")
    plt.legend()
    plt.subplot(513)
    plt.plot(yy, circDist, label="Gamma")
    plt.plot(yy, ellipse, label="Ideal")
    plt.legend()
    plt.subplot(514)
    plt.plot(yy, inducedAoa, label="Alpha_i")

    plt.legend()
    plt.subplot(515)
    plt.plot(planx, plany, label="Planform")
    plt.legend()
    plt.xlabel("y [m]")
    plt.legend()
    plt.savefig(dirname + "/" + "span")

if __name__ == '__main__':
    pass
