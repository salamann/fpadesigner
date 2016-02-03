#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os


def readcsv(filename):
    """
    Reads csv file
    All numerical data will be converted to float, otherwise string

    filename - filepath and filename
    """
    import csv
    data = []
    csvfile = file(filename)
    reader = csv.reader(csvfile)
    for row in reader:
        tmp = []
        for col in row:
            try:
                tmp.append(float(col))
            except:
                tmp.append(col)
        data.append(tmp)
    return data


def open2read(filename):
    """
    :param filename:
    :return:
    """
    with open(os.path.relpath(filename), "r") as data1:
        lines = data1.readlines()
    # data1.close()
    return lines


def read_data(filename):
    """
    Read aerodynamics data generated by xlfr5
    :param filename:
    :return:
    """
    #Xfoil data: 12, XFLR5 data: 11
    start_row = 11
    data = open2read(filename)
    data = [i.rsplit() for i in data[start_row:]]
    data = np.array(data, dtype=float)
    return data


def gen_result(testWing, testBody, testTail):
    """

    """
    br ="""
"""
    f = open(testWing.dirname + "/" +'text.txt', 'w') # 書き込みモードで開く
    str1 = ""
    str1 += "CL = %0.3f" % testWing.CL

    str1 += br
    str1 += "CD0 = %0.4f" % testWing.CD0
    str1 += br
    str1 += "CDi = %0.4f" % testWing.Cdi
    str1 += br
    str1 += "CD = %0.4f" % testWing.CD
    str1 += br

    str1 += "Lift = %0.3f[kgf]" % testWing.L
    str1 += br
    str1 += "WingDrag = %0.3f[N]" % testWing.D
    str1 += br
    str1 += "WingPower = %0.3f[W]" % testWing.W
    str1 += br

    str1 += "Fair Power = %0.3f[W]" % testBody.fairringPower
    str1 += br
    str1 += "Frame Power = %0.3f[W]" % testBody.framePower
    str1 += br

    str1 += "htail Power = %0.3f[W]" % testTail.htailPower
    str1 += br
    str1 += "vtail Power = %0.3f[W]" % testTail.vtailPower
    str1 += br

    RequiredPower = testBody.framePower \
                    + testBody.fairringPower \
                    + testWing.W \
                    + testTail.htailPower \
                    + testTail.vtailPower
    str1 += "Required Power = %0.3f[W]" % RequiredPower
    str1 += br

    f.write(str1) # 引数の文字列をファイルに書き込む
    f.close() # ファイルを閉じる
    return RequiredPower


if __name__ == '__main__':
    pass
