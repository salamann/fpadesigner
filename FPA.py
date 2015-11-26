#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.interpolate import interp2d, interp1d
from scipy.linalg import solve
from scipy import optimize

import io_fpa
from aerodynamics_2d import read_xflr5_data


# Load wing configuration file and create "wing" instance
#
# -----------------------------------------------
# sourceFile  - source csv file created by following structure
# analyzeStep - analyze slice steps スパン方向の分割数
# -----------------------------------------------


def calc_thickness_of_wing(XFOILdirectory, chordArray2):
    """
    calculation wing thickness list
    """
    # open airfoil data
    data = io_fpa.open2read(u"{}\\foil.dat".format(XFOILdirectory))

    # make airfoil list
    xlist = [float(i.split()[0]) for i in data[1:]]
    ylist = [float(i.split()[1]) for i in data[1:]]

    # divide upper and lower
    zeropoint = None
    for i in range(len(xlist)):
        if xlist[i] == ylist[i]:
            zeropoint = i

    upperx = np.array(xlist[:zeropoint+1])[::-1]
    uppery = np.array(ylist[:zeropoint+1])[::-1]

    lowerx = np.array(xlist[zeropoint:])
    lowery = np.array(ylist[zeropoint:])

    # interpolate uppwer and lower file in order to be able to different yposition of both upper and lower
    linear_interp_upper = interp1d(upperx, uppery)
    linear_interp_lower = interp1d(lowerx, lowery)

    xx = np.linspace(0., 1., 100)
    newylower = linear_interp_lower(xx)
    newyupper = linear_interp_upper(xx)

    thickness = newyupper - newylower
    maxthickness = max(thickness)

    # make thickness list of span direction
    thickness = [i * maxthickness for i in chordArray2]
    return thickness

    # plt.plot(self.yy, self.thickness)
    # plt.savefig(self.dirname + "/" + "thickness")


def calc_chord(lambda1, lambda2, y1, y2, Cr, yy):
    """
    コード長を求めるメソッド
    :param lambda1: 一つ前のテーパ比
    :param lambda2: 1つ後のテーパ比
    :param y1: 一つ前のスパン
    :param y2: 1つ後のスパン
    :param Cr: ルートコード長
    :param yy: 求めるスパン長
    :return:
    """
    return lambda1 * Cr - (yy - y1) * (lambda1 - lambda2) / (y2 - y1) * Cr


class Wing(object):
    def __init__(self, sourceFile, halfStep, surface, aspect, optflag=0):
        self.aerodynamics_2d_data = None
        self.lift_slope_array = None
        self.zero_lift_angle_array = None
        self.local_angle = None
        self.xcp_array = None
        self.chord_cp_array = None
        self.cd0_array = None
        self.span = None
        self.cr = None
        self.span_lambda = None
        self.yy = None
        self.chord_array = None
        self.temperature = None
        self.velocity = None
        self.Re = None
        self.airDensity = None

        self._sourceFile = sourceFile
        self.halfStep = halfStep
        data = io_fpa.readcsv(self._sourceFile)
        #self.surface = data[5][0]
        #self.aspect = data[5][1]
        self.surface = surface
        self.aspect = aspect
        self.dihedral = data[5][2]
        self.XFOILdirectory = str(data[5][3])
        # if optflag == 1:
        #    pass#self.shapeData = directshape
        # elif optflag == 0:
        #    self.shapeData = data[9:]
        self.dirname = str(self.XFOILdirectory) + "S" + str(self.surface) + "AR" + str(self.aspect)
        if not os.path.isdir(self.dirname):
            os.mkdir(self.dirname)

        if optflag == 1:
            pass
        elif optflag == 0:
            self.shapeData = data[9:]
            self.wing_shape()

    # TODO: グラフを描くための機能が混じっている。リファクタリングが必要。
    def wing_shape(self):
        """
        calc the chord array according to slices
        """
        # スパン長の計算
        span = np.sqrt(self.aspect * self.surface)
        self.span = span

        # ルートコード長crの計算
        bunbo = 0.
        for i in range(0, len(self.shapeData)-1, 1):
            bunbo += (self.shapeData[i][1]/100. + self.shapeData[i+1][1]/100.) \
                     * (self.shapeData[i+1][0]/100. - self.shapeData[i][0]/100.)
        bunbo = bunbo * span / 2.
        cr = self.surface / bunbo
        self.cr = cr

        # 平面形のラムダとスパン位置のリスト作成
        spanlambda = [self.shapeData[i][1]/100. for i in range(len(self.shapeData))]
        spanratio = [self.shapeData[i][0]/100. for i in range(len(self.shapeData))]
        self.span_lambda = spanlambda

        # スパン位置
        y = [(self.span / 2.0) * np.cos((i + 1) * (np.pi / 2.0) / self.halfStep) for i in range(self.halfStep)]

        # コード長の配列 cosθ基準
        chord_array = []
        for yy in y:
            for i in range(len(spanlambda)-1):
                if spanratio[i] * self.span/2.0 <= yy < spanratio[i + 1]*self.span/2.0:
                    chord_array.append(calc_chord(spanlambda[i],
                                                  spanlambda[i + 1],
                                                  spanratio[i]*self.span/2.0,
                                                  spanratio[i + 1]*self.span/2.0,
                                                  cr,
                                                  yy))
        # スパン位置の配列
        self.yy = y
        # コード長の配列
        self.chord_array = np.array(chord_array)

        figx = cr * (1.0 - np.array(spanlambda))
        figy = self.span / 2.0 * np.array(spanratio)
        figy = list(figy)
        figy1 = figy[:]
        figy.reverse()
        figy2 = figy
        figyy = figy1 + figy2 + [figy2[len(figy2)-1]]
        figx = list(figx)
        for i in range(len(figy2)):
            figx.append(cr)
        figx.append(0)
        figxx = figx

        self.figxx = figxx
        self.figyy = figyy

        calc_thickness_of_wing(self.XFOILdirectory, self.chord_array)

    def set_temperature(self, temperature):
        self.temperature = temperature

    def set_velocity(self, velocity):
        self.velocity = velocity

    def calc_reynolds(self):
        kine_vis = 1.34 * 10 ** - 5. + 9.31477 * 10 ** - 8. * self.temperature
        self.Re = np.array(self.chord_array) * self.velocity / kine_vis
        self.airDensity = 1.28912 - 0.004122391 * self.temperature
        return self.Re

    def calc_lift_slope_and_zero_lift_array(self, start_angle=-3.0, end_angle=4.0):
        """
        揚力傾斜をRe数ごと(スパン方向ごと)に計算
        This method is to calculate lift slope
        for each span position using reynolds number.
        :param start_angle: float
        :param end_angle: float
        """
        cl = self.aerodynamics_2d_data["CL"]
        alpha = np.radians(self.aerodynamics_2d_data["alpha"])
        func = interp2d(alpha[0], np.arange(0, 1.01, 0.1), cl, kind='linear')
        data = func(alpha[0], self.Re / 10 ** 6.)
        """
        Calculate the range that is used for getting regression curve for lift slope
        """
        start_index = np.where(alpha[0] == np.radians([start_angle]))[0]
        end_index = np.where(alpha[0] == np.radians([end_angle]))[0]
        lift_slope_array = []
        zero_lift_angle_array = []
        for datum in data:
            slope, intercept, r_value, p_value, std_err = stats.linregress(alpha[0][start_index:end_index],
                                                                           datum[start_index:end_index])
            lift_slope_array.append(slope)

            zero_lift_angle_array.append(-intercept / slope)
        self.lift_slope_array = lift_slope_array
        self.zero_lift_angle_array = zero_lift_angle_array

    def set_local_angle(self):
        self.local_angle = [self.angle - i for i in np.degrees(self.inducedAoa)]

    def calc_cd0_array(self):
        """
        スパン方向のCD0の計算
        :return: None
        """
        local_angles = self.local_angle #吹き下ろしを弾いた後の迎角配列

        cd = self.aerodynamics_2d_data["CD"]
        alpha = np.radians(self.aerodynamics_2d_data["alpha"])
        func = interp2d(alpha[0], np.arange(0, 1.01, 0.1), cd, kind='linear')
        # data = func(alpha[0], self.Re / 10 ** 6.)
        cd0_array = []
        for local_angle, reynolds_number in zip(local_angles, self.Re):
            cd0_array.append(func(np.radians(local_angle), reynolds_number)[0])
        self.cd0_array = cd0_array

    def calc_xcp_array(self):
        """
        #Calculate center of pressure
        :return: None
        """
        local_angles = self.local_angle #吹き下ろしを弾いた後の迎角配列
        xcp = self.aerodynamics_2d_data["XCp"]
        alpha = np.radians(self.aerodynamics_2d_data["alpha"])
        func = interp2d(alpha[0], np.arange(0, 1.01, 0.1), xcp, kind='linear')
        # data = func(alpha[0], self.Re / 10 ** 6.)
        xcp_array = []
        for localangle, reynolds_number in zip(local_angles, self.Re):
            xcp_array.append(func(np.radians(localangle), reynolds_number)[0])
        self.xcp_array = xcp_array
        self.chord_cp_array = xcp_array * self.chord_array

    def calc_down_wash(self, thetas, An):
        """
        #吹き下ろしの計算
        #航空力学の基礎第2版 p.141 式3.97より
        """
        dwArray = []
        for theta in thetas:
            dw = sum([(2*i+1)*An[i] * np.sin((2*i+1)*theta) / np.sin(theta) for i in range(len(An))]) * self.velocity
            dwArray.append(dw)
        self.dwArray = dwArray

    def calc_induced_angle_of_attack(self):
        """
        calculation induced alpha [radians]
        """
        inducedAoa = np.array(self.dwArray) / self.velocity
        self.inducedAoa = inducedAoa

    def calc_CL_Cdi_CD(self, angle, oddOReven=1):
        """
        calculate induced drag.
        :param angle:
        :param oddOReven:
        :return:
        """
        self.angle = angle
        self.oddOReven = oddOReven
        # 航空力学の基礎第2版 p.142 式3.99の下の式より
        # slopeはRe数ごとに与えることにした

        # μの計算。
        # self.liftSlopeArray[i] -> 5.5 で航空力学の基礎第2版 p.147の計算になる
        self.calc_lift_slope_and_zero_lift_array()
        u = [self.lift_slope_array[i] * self.chord_array[i] / 4.0 / self.span for i in range(len(self.chord_array))]

        # 絶対迎角
        # self.calc_zero_lift_angle_array()
        # θの設定、θ[0]はpi/n, θ[1]は2pi/n, θ[2]は3pi/n・・・θ[ラスト]はpi/2
        thetas = [float(i)/self.halfStep * np.pi / 2.0 for i in range(1, self.halfStep+1)]

        # 行列を作る計算
        left_matrix = []
        right_matrix = []
        i = 0
        for theta in thetas:
            tmp = []
            for n in range(1, self.halfStep * 2, oddOReven+1):
                tmp.append((n * u[i] + np.sin(theta)) * np.sin(n*theta))
            right_hand_term = u[i] * np.sin(theta)
            left_matrix.append(np.array(tmp) / right_hand_term)

            # これを入れると航空力学の基礎p.147の連立方程式が再現できる
            # print tmp / righthand

            # absoluteAlpha = 1.0で航空力学の基礎第2版 p.147の計算になる
            absoluteAlpha = np.radians(self.angle) - self.zero_lift_angle_array[i]
            # absoluteAlpha = 1.0
            right_matrix.append(absoluteAlpha)

            i += 1

        # 連立方程式を行列を使って解く
        An = solve(np.array(left_matrix), np.array(right_matrix))
        # これをいれると航空力学の基礎p.147のAnの値がでる

        # calc CL
        self.CL = np.pi * self.aspect * An[0]

        # calc Cdi
        sigma = 0
        # 航空力学の基礎第2版 p.147 式3.122の下の式より3からはじめる

        j = 3
        for i in range(self.oddOReven, self.halfStep):
            sigma += j * An[i] ** 2.0 / An[0] ** 2.0
            j += 2

        self.Cdi = (1.0 + sigma) * self.CL ** 2.0 / np.pi / self.aspect

        # down washと誘導迎角の計算
        self.calc_down_wash(thetas, An)
        self.calc_induced_angle_of_attack()
        # new angleは吹き下ろしを考慮した迎角
        # newangle = self.angle - np.degrees(self.inducedAoa)

        self.set_local_angle()
        self.calc_cd0_array()
        # CD0の計算
        # calc CD2
        deltaD = []
        for i in range(len(self.chord_array)):
            if i == 0:
                deltaD.append((self.cr * self.span_lambda[len(self.span_lambda) - 1] + self.chord_array[0]) *
                              (self.span / 2.0 - self.yy[i]) * 0.5 * self.cd0_array[0])
            else:
                deltaD.append((self.chord_array[i - 1] + self.chord_array[i]) *
                              (self.yy[i - 1] - self.yy[i]) * 0.5 *
                              (self.cd0_array[i - 1] + self.cd0_array[i]) / 2.0)
        D2 = 0.5 * self.airDensity * self.velocity ** 2.0 * sum(deltaD)
        dCD = D2 / (0.5 * self.airDensity * self.velocity ** 2.0 * self.surface / 2.0)

        # 風圧中心の計算
        self.calc_xcp_array()

        # スパン方向の揚力係数・循環・揚力計算
        circDist = []
        clDist = []
        j = 0
        for theta in thetas:
            summation = sum([An[i] * np.sin((2*i + 1) * theta) for i in range(0, self.halfStep)])

            # 循環の計算
            circ = 2.0 * self.span * self.velocity * summation
            circDist.append(circ)

            # 局所揚力係数の計算
            cllocal = 4.0 * self.span / self.chord_array[j] * summation
            clDist.append(cllocal)
            j += 1

        self.circDist = circDist
        self.clDist = clDist

        # 循環分布からΓ-yの楕円等価面積を求める
        # 循環の積分
        # numpy.trapzは台形近似で区分求積する。
        # Ellipse distribution will be calculated by circulation distribution.
        # numerical integration will be ensured by numpy.trapz using trapezoital approzimation.
        sum_gamma1 = -np.trapz(self.circDist,self.yy)
        minor_axis = 8./np.pi/self.span * sum_gamma1
        self.ellipse = minor_axis*(1.-np.array(self.yy)**2./(self.span/2.)**2.)**0.5
        self.eval_func = -np.trapz((self.circDist - self.ellipse)**2., self.yy)
        # print self.eval_func

        dL = 0.5 * self.airDensity * self.velocity ** 2.0 * self.chord_array * clDist

        self.dL = dL

        # 吹き下ろしの計算
        # print "---dw---",self.dwArray
        # 誘導抵抗以外の抵抗：CD0
        self.CD0 = dCD
        self.CD = dCD + self.Cdi

        # 揚力の計算：単位はNで出る

        # def calc_L(self):
        self.L = 0.5 * self.airDensity * self.velocity ** 2.0 * self.surface * self.CL * \
                 np.cos(np.radians(self.dihedral))
        self.L = self.L / 9.80665

        # 抵抗の計算：単位はNで出る
        # def calc_D(self):
        self.D = 0.5 * self.airDensity * self.velocity ** 2.0 * self.surface * self.CD

        # ワット数の計算：単位はｗ
        # def calc_W(self):
        self.W = self.D * self.velocity

    #TODO: これはよく分からない。たぶん循環分布を楕円型にして最適化を行うものと思うが,空力データそのものとは関係ない
    def optimize_circulation(self, x):
        """
        Calculating integration of circulation
        :param x:
        :return:
        """
        self.shapeData = [[0.0, 100.0],
                          [40.0, 100],
                          [50., x[0]],
                          [70, x[1]],
                          [82, x[2]],
                          [90, x[3]],
                          [95, x[4]],
                          [100.0, x[5]]]
        self.wing_shape()
        self.calc_reynolds()
        self.calc_lift_slope_and_zero_lift_array()
        self.calc_CL_Cdi_CD(3.)
        print "calculating...", round(self.eval_func*1000, 1),\
            round(self.W, 1),\
            round(x[0], 1),\
            round(x[1], 1),\
            round(x[2], 1),\
            round(x[3], 1),\
            round(x[4], 1),\
            round(x[5], 1),\
            round(self.L*9.80665/self.D, 2), \
            round(self.chord_array[0], 2)
        return round(self.eval_func*1000, 1),\
               round(self.W, 1),\
               round(x[0], 1),\
               round(x[1], 1),\
               round(x[2], 1),\
               round(x[3], 1),\
               round(x[4], 1),\
               round(x[5], 1),\
               round(self.L*9.80665/self.D, 2), \
               round(self.chord_array[0], 2)
        #
        # Use module shown below if you do optimization
        #
        # return self.eval_func*10
        # return 1./round(self.L*9.80665/self.D,2)

    # TODO: solve_CLは，与えられた重量に対して釣り合うCLを探すプログラム。空力データとは関係ないので移動させる
    def solve_CL(self, angle):
        """
        This method is to calculate lift coefficient,
        where wing area, airspeed, and weight are given.
        """
        self.calc_reynolds()
        self.calc_lift_slope_and_zero_lift_array()
        # self.calc_zero_lift_angle_array()
        self.calc_CL_Cdi_CD(angle)

        Cw = self.weight * 9.81 / (0.5 * self.airDensity * self.velocity ** 2.0 * self.surface)
        print "solving CL of constant weight..."
        return Cw - self.CL * np.cos(np.radians(self.dihedral))

    # TODO: weightは，与えられた重量に対して釣り合うCLを探すプログラム。空力データとは関係ないので移動させる
    def calc_weight(self, weight):
        self.weight = weight

    # TODO: calc_withconstWeightは，与えられた重量に対して釣り合うCLを探すプログラム。空力データとは関係ないので移動させる
    def calc_withconstWeight(self, objweight, velocity, temperature):
        self.calc_reynolds()
        self.calc_weight(objweight)

        optimize.brenth(self.solve_CL, -5, 10)

    ##        L = self.L/9.81

    def set_aerodynamcis_data(self):
        self.aerodynamics_2d_data = read_xflr5_data(self.XFOILdirectory)

    # TODO: グラフを描くための機能と迎角を振るための機能が混じっている
    def calc_variedaoa(self, velocity, temperature, aoaarray):
        """
        csvfile, number of cell, design cruise speed, ambient temperature
        :param velocity:
        :param temperature:
        :param aoaarray:
        :return:
        """

        #testWing = wing(wingcsv,ncell)
        self.set_temperature(temperature)
        self.set_velocity(velocity)
        self.calc_reynolds()
        self.set_aerodynamcis_data()
        """機体の特性を出す"""
        CLarray = []
        CDarray = []
        for i in aoaarray:
            print "alpha = {} [deg]".format(str(i))
            self.calc_lift_slope_and_zero_lift_array()
            self.calc_CL_Cdi_CD(i)
            CLarray.append(self.CL)
            CDarray.append(self.CD)

        self.CLarray = CLarray
        self.CDarray = CDarray

        """maxL/Dの線を引くためのリスト生成"""
        j = 0
        maxslope = None
        for i in np.array(CLarray)/np.array(CDarray):
            if i == max(np.array(CLarray)/np.array(CDarray)):
                maxslope = i
            j += 1
        self.maxslope = maxslope

        xlist = list(CDarray)
        ylist = list(np.array(self.CDarray)*self.maxslope)
        xlist.insert(0, 0)
        ylist.insert(0, 0)
        self.xmaxLDline = xlist
        self.ymaxLDline = ylist

    def calc_planform(self):
        """
        calculation of planform data for drawing planform
        :return:
        """
        y1 = [self.chord_array[len(self.xcp_array) - 1] * (1.0 - self.xcp_array[len(self.xcp_array) - 1]) - self.chord_array[len(self.xcp_array) - i] * (1.0 - self.xcp_array[len(self.xcp_array) - i]) for i in range(1, len(self.xcp_array) + 1)]
        x1 = [self.yy[len(self.xcp_array) - i] for i in range(1, len(self.xcp_array) + 1)]
        y2 = [self.chord_array[len(self.xcp_array) - 1] * self.xcp_array[len(self.xcp_array) - 1] + self.chord_array[len(self.xcp_array) - i] * self.xcp_array[len(self.xcp_array) - i] for i in range(1, len(self.xcp_array) + 1)]

        x2 = x1 + x1[::-1]
        y2 = y1 + y2[::-1]
        plt.figure(figsize=(12, 4))
        self.planx = x2
        self.plany = y2
        plt.plot(x2, y2)
        xcp0 = self.chord_array[len(self.xcp_array) - 1] * (1.0 - self.xcp_array[len(self.xcp_array) - 1])
        plt.plot([0, self.span / 2.], [xcp0, xcp0])
        plt.axis("equal")
        plt.xlabel("y [m]")
        plt.ylabel("x [m]")
        plt.legend(("planform", "pressure center"))
        plt.savefig(self.dirname + "/" + "testwing")
#        pl.clf()


if __name__ == '__main__':
    pass
