#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve
from scipy import stats
from scipy.interpolate import interp1d, interp2d

import io_fpa
from aerodynamics_2d import read_xflr5_data
import post_process_operation
import weight
# Load wing configuration file and create "wing" instance
#
# -----------------------------------------------
# sourceFile  - source csv file created by following structure
# analyzeStep - analyze slice steps スパン方向の分割数
# -----------------------------------------------
#
# from runner import velocity, temperature, aoaarray, zz, testWing, res, x1, x2, x3, x4, x5, x6, x


class Wing(object):
    def __init__(self, sourceFile, halfStep, surface, aspect, optflag=0):
        self._sourceFile = sourceFile
        self.halfStep = halfStep
        data = io_fpa.readcsv(self._sourceFile)
        #self.surface = data[5][0]
        #self.aspect = data[5][1]
        self.surface = surface
        self.aspect = aspect
        self.dihedral = data[5][2]
        self.XFOILdirectory = str(data[5][3])
##        if optflag == 1:
##            pass#self.shapeData = directshape
##        elif optflag == 0:
##            self.shapeData = data[9:]
        import os
        self.dirname = str(self.XFOILdirectory) + "S" + str(self.surface) + "AR" + str(self.aspect)
        try:
            os.mkdir(self.dirname)
        except:
            pass

        if optflag == 1:
            pass
        elif optflag == 0:
            self.shapeData = data[9:]
            self.wingshape()
        self.aerodynamics_2d_data = None
        self.lift_slope_array = None
        self.zero_lift_angle_array = None
        self.local_angle = None

    def wingshape(self):
        """
        calc the chord array according to slices
        """
        #スパン長の計算
        span = np.sqrt(self.aspect * self.surface)
        self.span = span

        #ルートコード長crの計算
        bunbo = 0.
        for i in range(0, len(self.shapeData)-1, 1):
            bunbo += (self.shapeData[i][1]/100. + self.shapeData[i+1][1]/100.) * (self.shapeData[i+1][0]/100. - self.shapeData[i][0]/100.)
        bunbo = bunbo * span / 2.
        cr = self.surface / bunbo
        self.cr = cr


        #平面形のラムダとスパン位置のリスト作成
        spanlambda = [self.shapeData[i][1]/100. for i in range(len(self.shapeData))]
        spanratio = [self.shapeData[i][0]/100. for i in range(len(self.shapeData))]
        self.spanlambda = spanlambda

        #スパン位置
        y = [(self.span / 2.0) *np.cos((i + 1) * (np.pi / 2.0) / (self.halfStep)) for i in range(self.halfStep)]

        #コード長の配列 cosθ基準
        chordArray2 = []
        for yy in y:
            for i in range(len(spanlambda)-1):
                if spanratio[i]*self.span/2.0 <= yy and yy < spanratio[i + 1]*self.span/2.0:
                    chordArray2.append(self.calc_chord(spanlambda[i],
                                                       spanlambda[i + 1],
                                                       spanratio[i]*self.span/2.0,
                                                       spanratio[i + 1]*self.span/2.0,
                                                       cr, yy))
        #スパン位置の配列
        self.yy = y
        #コード長の配列
        self.chordArray2 = np.array(chordArray2)

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

        self.calc_wingthickness()
    """
    コード長を求めるメソッド
    lambda1:一つ前のテーパ比,lambda2:1つ後のテーパ比,y1:一つ前のスパン,y2:1つ後のスパン
    Cr:ルートコード長、yy:求めるスパン長
    """
    def calc_chord(self, lambda1, lambda2, y1, y2, Cr, yy):
        return lambda1 * Cr - (yy - y1) * (lambda1 - lambda2) / (y2 - y1) * Cr

    def calc_wingthickness(self):
        """calculation wing thickness list"""
        from scipy.interpolate import interp1d

        #open airfoil data
        data = io_fpa.open2read(self.XFOILdirectory + "/" + "foil.dat")

        #make airfoil list
        xlist = [float(i.split()[0]) for i in data[1:]]
        ylist = [float(i.split()[1]) for i in data[1:]]

        #divide upper and lower
        for i in range(len(xlist)):
            if xlist[i] == ylist[i]:
                zeropoint = i
        upperx = np.array(xlist[:zeropoint+1])[::-1]
        uppery = np.array(ylist[:zeropoint+1])[::-1]

        lowerx = np.array(xlist[zeropoint:])
        lowery = np.array(ylist[zeropoint:])

        #interpolate uppwer and lower file in order to be able to different yposition of both upper and lower
        linear_interp_upper = interp1d(upperx, uppery)
        linear_interp_lower = interp1d(lowerx, lowery)
##
        xx = np.linspace(0., 1., 100)
        newylower = linear_interp_lower(xx)
        newyupper = linear_interp_upper(xx)

        thickness = newyupper - newylower
        maxthickness = max(thickness)
        #thickness36 = thickness[36]
        #print thickness36

        #make thickness list of span direction
        self.thickness = [i * maxthickness for i in self.chordArray2]
        #self.thick36 = [i * thickness36 for i in self.chordArray2]
        #print self.thick36

        plt.plot(self.yy, self.thickness)
        plt.savefig(self.dirname + "/" + "thickness")

    # Analyzes the wing object
    # -----------------------------------------------
    # velocity   - airspeed
    # temperature - temperature
    # angle      - angle of flight
    # -----------------------------------------------
    #
    #
    def calc_reynolds(self, velocity=0, temperature=0):
        kine_vis = 1.34 * 10 ** - 5. + 9.31477 * 10 ** - 8. * self.temperature
        self.Re = np.array(self.chordArray2) * self.velocity / kine_vis
        self.airDensity = 1.28912 - 0.004122391 * self.temperature
        return self.Re

    def calc_lift_slope_array(self):
        """
        揚力傾斜をRe数ごと(スパン方向ごと)に計算
        This method is to calculate lift slope
        for each span position using reynolds number.
        """
        cl = self.aerodynamics_2d_data["CL"]
        alpha = np.radians(self.aerodynamics_2d_data["alpha"])
        func = interp2d(alpha[0], np.arange(0, 1.01, 0.1), cl, kind='linear')
        data = func(alpha[0], self.Re / 10 ** 6.)
        """
        Calculate the range that is used for getting regression curve for lift slope
        """
        start_index = np.where(alpha[0] == np.radians([-3.0]))[0]
        end_index = np.where(alpha[0] == np.radians([4.0]))[0]
        lift_slope_array = []
        zero_lift_angle_array = []
        for datum in data:
            slope, intercept, r_value, p_value, std_err = stats.linregress(alpha[0][start_index:end_index],
                                                                           datum[start_index:end_index])
            lift_slope_array.append(slope)

            zero_lift_angle_array.append(-intercept / slope)
        self.lift_slope_array = lift_slope_array
        self.zero_lift_angle_array = zero_lift_angle_array


        # datas = [self.calc_interpolate(i) for i in self.Re / 10 ** 6.]

        # lift_slope_array = []
        # for data in datas:
        #     start_angle = -3 #deg
        #     end_angle = 4 #deg
        #
        #     slope, intercept, r_value, p_value, std_err = stats.linregress(np.radians(data.transpose()[0][start_angle+10:end_angle+10]), data.transpose()[1][start_angle+10:end_angle+10])
        #     lift_slope_array.append(slope)
        # self.lift_slope_array = lift_slope_array

    # def calc_zeroLiftAngle(self, data):
    #     """
    #     This method is to calculate zero lift angle of attack.
    #     :param data:
    #     :return:
    #     """
    #     from scipy import stats
    #
    #     #揚力傾斜を出すのに使う迎角範囲
    #     start_angle = 0 #deg
    #     end_angle = 5 #deg
    #
    #     slope, intercept, r_value, p_value, std_err = stats.linregress(np.radians(data.transpose()[0][start_angle+10:end_angle+10]),data.transpose()[1][start_angle+10:end_angle+10])
    #     self.zeroLiftAngle = - data[10][1] / slope
    #
    #     return self.zeroLiftAngle

    #ゼロ揚力角の配列の計算
    # def calc_zero_lift_angle_array(self):
    #     cl = self.aerodynamics_2d_data["CL"]
    #     alpha = np.radians(self.aerodynamics_2d_data["alpha"])
    #     func = interp2d(alpha[0], np.arange(0, 1.01, 0.1), cl, kind='linear')
    #     data = func(alpha[0], self.Re / 10 ** 6.)
    #     """
    #     Calculate the range that is used for getting regression curve for lift slope.
    #     """
    #     start_index = np.where(alpha[0] == np.radians([-3.0]))[0]
    #     end_index = np.where(alpha[0] == np.radians([4.0]))[0]
    #     zeroliftangleArray = []
    #     for datum in data:
    #         slope, intercept, r_value, p_value, std_err = stats.linregress(alpha[0][start_index:end_index],
    #                                                                        datum[start_index:end_index])
    #         zeroliftangleArray.append(-intercept / slope)
    #     self.zeroliftangleArray = zeroliftangleArray

    # def calc_interpolate(self, reynolds):
    #     rey1 = round(reynolds, 1)
    #     rey2 = round(reynolds, 2)
    #
    #     #レイノルズ数の小数点第2桁が0のとき,つまり内挿しなくていいとき 0.80とか
    #     if rey1 == rey2:
    #         txtFile = str(rey1)
    #         if len(txtFile)==1:
    #             txtFile+='.'
    #         for i in range(4-len(txtFile)):
    #             txtFile = txtFile+'0'
    #
    #         self.XFOILfile = '/'+txtFile+'.txt'
    #         data = io_fpa.read_data(self.XFOILdirectory + self.XFOILfile)
    #     #内挿すべきとき たとえば0.84
    #     else:
    #         txtFile = str(rey1-0.1)
    #         if len(txtFile)==1:
    #             txtFile+='.'
    #         for i in range(4-len(txtFile)):
    #             txtFile = txtFile+'0'
    #         self.XFOILfile = '/'+txtFile+'.txt'
    #         data1 = io_fpa.read_data(self.XFOILdirectory + self.XFOILfile)
    #
    #         txtFile = str(rey1)
    #         if len(txtFile)==1:
    #             txtFile+='.'
    #         for i in range(4-len(txtFile)):
    #             txtFile = txtFile+'0'
    #         self.XFOILfile = '/'+txtFile+'.txt'
    #         data2 = io_fpa.read_data(self.XFOILdirectory + self.XFOILfile)
    #
    #         data = (data2-data1)*(rey2-rey1)/(float(str(rey1 + 0.1))-rey1)+data1
    #
    #     return data

    def set_local_angle(self):
        self.local_angle = [self.angle - i for i in np.degrees(self.inducedAoa)]

    #スパン方向のCD0の計算
    def calc_CD0Array(self):
        local_angles = self.local_angle #吹き下ろしを弾いた後の迎角配列

        cd = self.aerodynamics_2d_data["CD"]
        alpha = np.radians(self.aerodynamics_2d_data["alpha"])
        func = interp2d(alpha[0], np.arange(0, 1.01, 0.1), cd, kind='linear')
        # data = func(alpha[0], self.Re / 10 ** 6.)
        cd0Array = []
        for localangle, reynolds_number in zip(local_angles, self.Re):
            cd0Array.append(func(localangle, reynolds_number)[0])
        self.cd0Array = cd0Array

        # start_index = np.where(alpha[0] == np.radians([-3.0]))[0]
        # end_index = np.where(alpha[0] == np.radians([4.0]))[0]
        # for datum in data:
        #     slope, intercept, r_value, p_value, std_err = stats.linregress(alpha[0][start_index:end_index],
        #                                                                    datum[start_index:end_index])
        # self.lift_slope_array = lift_slope_array
        # self.zero_lift_angle_array = zero_lift_angle_array
        #
        #
        #
        # datas = [self.calc_interpolate(i) for i in self.Re / 10 ** 6.]
        # cd0Array = []
        # for i in range(len(self.Re)):
        #     #datas[i].transpose() #[0]と[2]でinterpolateする。
        #     linear_interp = interp1d(datas[i].transpose()[0], datas[i].transpose()[2])
        #     cd0Array.append(float(linear_interp(localangle[i])))
        # self.cd0Array = cd0Array

    #Calculate center of pressure
    def calc_xcp(self):
        local_angles = self.local_angle #吹き下ろしを弾いた後の迎角配列
        xcp = self.aerodynamics_2d_data["XCp"]
        alpha = np.radians(self.aerodynamics_2d_data["alpha"])
        func = interp2d(alpha[0], np.arange(0, 1.01, 0.1), xcp, kind='linear')
        # data = func(alpha[0], self.Re / 10 ** 6.)
        xcpArray = []
        for localangle, reynolds_number in zip(local_angles, self.Re):
            xcpArray.append(func(localangle, reynolds_number)[0])
        self.xcpArray = xcpArray
        self.chordcpArray = xcpArray * self.chordArray2


        # from scipy.interpolate import interp1d
        # datas = [self.calc_interpolate(i) for i in self.Re / 10 ** 6.]
        #
        # xcpArray = []
        # for i in range(len(self.Re)):
        #     linear_interp = interp1d(datas[i].transpose()[0], datas[i].transpose()[9])
        #     xcpArray.append(float(linear_interp(local_angles[i])))
        #
        # self.xcpArray = xcpArray
        # self.chordcpArray = xcpArray * self.chordArray2

    #吹き下ろしの計算
    #航空力学の基礎第2版 p.141 式3.97より
    def calc_downwash(self, thetas, An):
        dwArray = []
        for theta in thetas:
            dw = sum([(2*i+1)*An[i] * np.sin((2*i+1)*theta) / np.sin(theta) for i in range(len(An))]) * self.velocity
            dwArray.append(dw)
        self.dwArray = dwArray

    """
    calculation induced alpha [radians]
    """
    def calc_inducedAoa(self):
        inducedAoa = np.array(self.dwArray) / self.velocity
        self.inducedAoa = inducedAoa

    """
    calculation induced drag
    """
    def calc_CL_Cdi_CD(self, angle, oddOReven=1):
        self.angle = angle
        self.oddOReven = oddOReven
        #航空力学の基礎第2版 p.142 式3.99の下の式より
        #slopeはRe数ごとに与えることにした

        # μの計算。
        # self.liftSlopeArray[i] -> 5.5 で航空力学の基礎第2版 p.147の計算になる
        self.calc_lift_slope_array()
        u = [self.lift_slope_array[i] * self.chordArray2[i] / 4.0 / self.span for i in range(len(self.chordArray2))]

        # 絶対迎角
        # self.calc_zero_lift_angle_array()
        #θの設定、θ[0]はpi/n, θ[1]は2pi/n, θ[2]は3pi/n・・・θ[ラスト]はpi/2
        thetas = [float(i)/(self.halfStep) * np.pi / 2.0 for i in range(1, self.halfStep+1)]

        #行列を作る計算
        lmatrix = []
        rmatrix = []
        i = 0
        for theta in thetas:
            tmp = []
            for n in range(1, self.halfStep * 2, oddOReven+1):

                tmp.append((n * u[i] + np.sin(theta)) * np.sin(n*theta))
            righthand = u[i] * np.sin(theta)
            lmatrix.append(tmp / righthand)

            #これを入れると航空力学の基礎p.147の連立方程式が再現できる
            #print tmp / righthand

            # absoluteAlpha = 1.0で航空力学の基礎第2版 p.147の計算になる
            absoluteAlpha = np.radians(self.angle) - self.zero_lift_angle_array[i]
            #absoluteAlpha = 1.0
            rmatrix.append(absoluteAlpha)

            i += 1

        #連立方程式を行列を使って解く
        An = solve(np.array(lmatrix), np.array(rmatrix))
        #これをいれると航空力学の基礎p.147のAnの値がでる

        # calc CL
        self.CL = np.pi*self.aspect*An[0]

        # calc Cdi
        sigma = 0
        #航空力学の基礎第2版 p.147 式3.122の下の式より3からはじめる

        j = 3
        for i in range(self.oddOReven, self.halfStep):
            sigma += j * An[i] ** 2.0 / An[0] ** 2.0
            j += 2

        self.Cdi = (1.0 + sigma)*self.CL**2.0/np.pi/self.aspect

        #donw washと誘導迎角の計算
        self.calc_downwash(thetas,An)
        self.calc_inducedAoa()
        #new angleは吹き下ろしを考慮した迎角
##        newangle = self.angle - np.degrees(self.inducedAoa)

        self.set_local_angle()
        self.calc_CD0Array()
        #CD0の計算
        #calc CD2
        deltaD = []
        for i in range(len(self.chordArray2)):
            if i == 0:
                deltaD.append((self.cr * self.spanlambda[len(self.spanlambda)-1] + self.chordArray2[0]) * (self.span/2.0 - self.yy[i]) * 0.5 * self.cd0Array[0])
            else:
                deltaD.append((self.chordArray2[i-1] + self.chordArray2[i]) * (self.yy[i-1] - self.yy[i]) * 0.5 * (self.cd0Array[i-1] + self.cd0Array[i])/2.0)
        D2 = 0.5 * self.airDensity * self.velocity ** 2.0 * sum(deltaD)
        dCD = D2 / (0.5 * self.airDensity * self.velocity ** 2.0 * self.surface / 2.0)

        #風圧中心の計算
        self.calc_xcp()

        #スパン方向の揚力係数・循環・揚力計算
        circDist = []
        clDist = []
        j = 0
        for theta in thetas:
            summation = sum([An[i] * np.sin((2*i + 1) * theta) for i in range(0,self.halfStep)])

            #循環の計算
            circ = 2.0 * self.span * self.velocity * summation
            circDist.append(circ)

            #局所揚力係数の計算
            cllocal = 4.0 * self.span / self.chordArray2[j] * summation
            clDist.append(cllocal)
            j += 1

        self.circDist = circDist
        self.clDist = clDist

        #循環分布からΓ-yの楕円等価面積を求める
		#循環の積分
        #numpy.trapzは台形近似で区分求積する。
        #Ellipse distribution will be calculated by circulation distribution.
        #numerical integration will be ensured by numpy.trapz using trapezoital approzimation.
        sum_gamma1 = -np.trapz(self.circDist,self.yy)
        minor_axis = 8./np.pi/self.span * sum_gamma1
        self.ellipse = minor_axis*(1.-np.array(self.yy)**2./(self.span/2.)**2.)**0.5
        self.eval_func = -np.trapz((self.circDist - self.ellipse)**2.,self.yy)
        #print self.eval_func



        dL = 0.5 * self.airDensity * self.velocity **2.0 * self.chordArray2 * clDist

        self.dL = dL

        #吹き下ろしの計算
        #print "---dw---",self.dwArray
        #誘導抵抗以外の抵抗：CD0
        self.CD0 = dCD
        self.CD = dCD + self.Cdi



    #揚力の計算：単位はNで出る

##    def calc_L(self):
        self.L = 0.5*self.airDensity*self.velocity**2.0*self.surface*self.CL*np.cos(np.radians(self.dihedral))
        self.L = self.L / 9.80665

    #抵抗の計算：単位はNで出る
##    def calc_D(self):
        self.D = 0.5*self.airDensity*self.velocity**2.0*self.surface*self.CD

    #ワット数の計算：単位はｗ
##    def calc_W(self):
        self.W = self.D * self.velocity

##
##    def calc_all(self, angle, Reynolds):
##        self.calc_liftSlope(Reynolds)
##        self.calc_zeroLiftAngle()
##        self.calc_CL_Cdi_CD(angle)


    #Calculating integration of circulation
    def opt_circ(self, x):
        #self.shapeData = [[0.0, 100.0], [40.0, 100], [50., x[0]],[60, x[1]],[70, x[2]],[80, x[3]],[90, x[4]],[95, x[5]], [100.0, x[6]]]
        self.shapeData = [[0.0, 100.0],
                          [40.0, 100],
                          [50., x[0]],
                          [70, x[1]],
                          [82, x[2]],
                          [90, x[3]],
                          [95, x[4]],
                          [100.0, x[5]]]
        self.wingshape()
        #print self.chordArray2
        self.calc_reynolds(self.velocity, self.temperature)
        self.calc_lift_slope_array()
        # self.calc_zero_lift_angle_array()
        self.calc_CL_Cdi_CD(3.)
        #print "calculating...",round(self.eval_func*1000,1),round(self.W,1),round(x[0],1),round(x[1],1),round(x[2],1),round(x[3],1),round(x[4],1),round(x[5],1),round(x[6],1),round(self.L*9.80665/self.D,2)
        print "calculating...", round(self.eval_func*1000,1),round(self.W,1),round(x[0],1),round(x[1],1),round(x[2],1),round(x[3],1),round(x[4],1),round(x[5],1),round(self.L*9.80665/self.D,2), round(self.chordArray2[0],2)
        return round(self.eval_func*1000,1),round(self.W,1),round(x[0],1),round(x[1],1),round(x[2],1),round(x[3],1),round(x[4],1),round(x[5],1),round(self.L*9.80665/self.D,2), round(self.chordArray2[0],2)
        #
        #Use module shonw below if you do optimization
        #
        #return self.eval_func*10
        #return 1./round(self.L*9.80665/self.D,2)

    def calcTrimdrag(self):
        airplaneCG=0.3

    def solve_CL(self, angle):
        """
        This method is to calculate lift coefficient,
        where wing area, airspeed, and weight are given.
        """
        self.calc_reynolds(self.velocity, self.temperature)
        self.calc_lift_slope_array()
        # self.calc_zero_lift_angle_array()
        self.calc_CL_Cdi_CD(angle)

        Cw = self.weight*9.81 / (0.5 * self.airDensity * self.velocity ** 2.0 *self.surface)
        print "solving CL of constant weight..."
        return Cw - self.CL * np.cos(np.radians(self.dihedral))

    def calc_weight(self, weight):
        self.weight = weight

    def calc_withconstWeight(self, objweight, velocity, temperature):
        self.calc_reynolds(velocity, temperature)
        self.calc_weight(objweight)
        from scipy import optimize

        optimize.brenth(self.solve_CL, -5, 10)

    ##        L = self.L/9.81

    def set_aerodynamcis_data(self):
        self.aerodynamics_2d_data = read_xflr5_data(self.XFOILdirectory)

    def calc_variedaoa(self, velocity, temperature, aoaarray):
        """
        csvfile, number of cell, design cruise speed, ambient temperature
        """
        #testWing = wing(wingcsv,ncell)
        self.calc_reynolds(velocity, temperature)
        self.set_aerodynamcis_data()
        """機体の特性を出す"""
        CLarray = []
        CDarray = []
        for i in aoaarray:
            print "alpha = " + str(i) +" deg"
            self.calc_lift_slope_array()
            # self.calc_zero_lift_angle_array()
            self.calc_CL_Cdi_CD(i)
            CLarray.append(self.CL)
            CDarray.append(self.CD)

        self.CLarray = CLarray
        self.CDarray = CDarray

        """maxL/Dの線を引くためのリスト生成"""
        j = 0
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

    """calculation of planform data for drawing planform"""

    def calc_planform(self):
        y1 = [self.chordArray2[len(self.xcpArray)-1] * (1.0 - self.xcpArray[len(self.xcpArray)-1]) - self.chordArray2[len(self.xcpArray)-i] * (1.0 - self.xcpArray[len(self.xcpArray)-i]) for i in range(1,len(self.xcpArray)+1)]
        x1 = [self.yy[len(self.xcpArray)-i] for i in range(1, len(self.xcpArray)+1)]
        y2 = [self.chordArray2[len(self.xcpArray)-1] * self.xcpArray[len(self.xcpArray)-1] + self.chordArray2[len(self.xcpArray)-i] * self.xcpArray[len(self.xcpArray)-i] for i in range(1,len(self.xcpArray)+1)]

        x2 = x1 + x1[::-1]
        y2 = y1 + y2[::-1]
        plt.figure(figsize=(12, 4))
        self.planx = x2
        self.plany = y2
        plt.plot(x2, y2)
        xcp0 = self.chordArray2[len(self.xcpArray)-1] * (1.0 - self.xcpArray[len(self.xcpArray)-1])
        plt.plot([0, self.span / 2.], [xcp0, xcp0])
        plt.axis("equal")
        plt.xlabel("y [m]")
        plt.ylabel("x [m]")
        plt.legend(("planform", "pressure center"))
        plt.savefig(self.dirname + "/" + "testwing")
#        pl.clf()


class Body(object):
    def __init__(self, velocity, temperature):
        self.temperature = temperature
        self.velocity = velocity
        kine_vis = 1.34 * 10 ** - 5. + 9.31477 * 10 ** - 8. * self.temperature
        self.airDensity = 1.28912 - 0.004122391 * self.temperature
        self.dynpres = 0.5 * self.airDensity * self.velocity ** 2.0

    def fairdragcalc(self, fairArea):
        """CD=0.10 : F-TEC technical report"""
        """CD=0.20 : T-MIT windtunnel data"""
        CDfair = 0.15
        self.fairringDrag = CDfair * fairArea * self.dynpres
        self.fairringPower = self.fairringDrag * self.velocity

    def framedragcalc(self, framearea):
        CDframe = 0.1
        self.frameDrag = CDframe * framearea * self.dynpres
        self.framePower = self.frameDrag * self.velocity


class Tail(object):
    def __init__(self, velocity, temperature):
        self.temperature = temperature
        self.velocity = velocity
        kine_vis = 1.34 * 10 ** - 5. + 9.31477 * 10 ** - 8. * self.temperature
        self.airDensity = 1.28912 - 0.004122391 * self.temperature
        self.dynpres = 0.5 * self.airDensity * self.velocity ** 2.0

    def calc_htaildrag(self, htailArea):
        CDhtail = 0.008
        self.htailDrag = CDhtail * htailArea * self.dynpres
        self.htailPower = self.htailDrag * self.velocity

    def calc_vtaildrag(self, vtailArea):
        CDvtail = 0.008
        self.vtailDrag = CDvtail * vtailArea * self.dynpres
        self.vtailPower = self.vtailDrag * self.velocity

if __name__ == '__main__':
    pass