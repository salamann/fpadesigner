#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def calc_weight(span, airfoil):
    if airfoil == "DAE31":
        ff = 1.04
    elif airfoil == "DAE21" or airfoil == "FX76-MP140":
        ff = 1.02
    elif airfoil == "FX76-MP160":
        ff = 1.00
    else:
        raise ValueError("There is no wing foil you selected.")

    weight_wing = ff * 0.0101 * span ** 2. + 0.1137 *span + 1.1324

    #---------------
    #weight_pilot include the following items.
    # human being
    # clothes
    # shoes
    # helmets
    # eyeware
    #---------------
    weight_pilot = 60.

    #---------------
    #weight_others includes the following components and items.
    # propellers
    # frame
    # chains or shafts
    # crank
    # pedal
    # tail pipe
    # electrical equipment
    # battely
    # pilot seat
    # water bottle
    # camera
    # fairing
    # tires
    #---------------
    weight_others = 16.#15.699

    weight_wing2nd = wing2(span)
    weight_htail = htail(3)
    weight_vtail = vtail(2)
    summation = sum([weight_wing, weight_pilot, weight_others, weight_wing2nd, weight_htail, weight_vtail])

    return summation


def wing2(span):
    max_code = 1.2        #矩形部での翼弦長[m]
    min_code = 0.6        #翼端での翼弦長[m]
    spar_pos = 0.3        #桁位置
    rectangle_span = 10   #矩形部長さ[m]
    taper_span = 10       #テーパ部長さ[m]

    #spanが変わったときに重量感度を出すため
    rectangle_span = span /2.
    taper_span = span/2.

    airfoil_area = 0.098201235     #翼型の面積
    front_edge_length = 0.83781283 #前縁部長さ[m]
    rib_cap_length = 0.866192514   #リブキャップ長さ[m]
    spar_cap_area = 0.009646018    #桁キャップ断面積[m2]
    rear_edge_area = 0.001086871   #後縁材断面積[m2]
    rear_cap_area = 0.001778896    #後縁キャップ断面積[m2]

    rho_ib = 28         #スタイロ密度(IB)
    rho_balsa = 140     #バルサ密度
    rho_rear_edge = 140 #後縁材密度

    rib_area_coeff = 0.97115806 #リブ面積係数(翼型の面積のうちリブ材の占める割合)

    rib_thick = 0.006        #リブ厚み[m]
    plank_thick = 0.003      #プランク厚み[m]
    balsa_thick = 0.001      #バルサ厚み[m]
    stringer_area = 0.000008 #ストリンガー断面積[m2]
    stringer_num = 8         #ストリンガー本数[本]

    limit_rect_rib2rib = 0.30 #矩形部のリブ間隔をこれ以下にする[m]
    limit_tape_rib2rib = 0.30 #テーパ部リブ間隔をこれ以下にする[m]

    bond = 0.35 #接着剤重量[kg/m]

    rib_area =  airfoil_area * max_code*max_code #矩形部の翼型面積
    rect_rib_num = rectangle_span//limit_rect_rib2rib +2  #リブ枚数(矩形部)[枚]
    tape_rib_num = taper_span//limit_tape_rib2rib       #リブ枚数(テーパ部)[枚/片翼](矩形含まず）
    rect_rib2rib = rectangle_span / (rect_rib_num-1) #矩形部のリブ間隔[m]
    tape_rib2rib = taper_span / tape_rib_num     #テーパ部リブ間隔[m]

    max_code2 = max_code * max_code

    rib_weight = rib_area * rib_area_coeff * rib_thick * rho_ib #リブ本体
    rib_cap_weight = rib_cap_length*max_code * rib_thick * balsa_thick * rho_balsa #リブキャップ
    spar_cap_weight = spar_cap_area * max_code2 * balsa_thick * rho_balsa ##桁キャップ
    rear_cap_weight = rear_cap_area * max_code2 * balsa_thick * rho_balsa ##後縁キャップ

    def sum_taper_scale(num):
        length_scale = 0
        area_scale = 0
        for i in range(1, int(num)):
            taper_code = max_code-(max_code-min_code)*tape_rib2rib*i/taper_span
            length_scale += taper_code
            area_scale += (taper_code/max_code)*(taper_code/max_code)
        return [length_scale, area_scale, taper_code]

    taper_scale = sum_taper_scale(tape_rib_num)

    sum_rib  = rib_weight     *(rect_rib_num+taper_scale[1]*2) #リブ本体
    sum_rib += rib_cap_weight *(rect_rib_num+taper_scale[0]*2) #リブキャップ
    sum_rib += spar_cap_weight*(rect_rib_num+taper_scale[1]*2) #桁キャップ
    sum_rib += rear_cap_weight*(rect_rib_num+taper_scale[1]*2) #後縁キャップ

    sum_plank  = front_edge_length*max_code *rectangle_span*plank_thick*rho_ib
    sum_plank += front_edge_length*(max_code+min_code)/2*taper_span*plank_thick*rho_ib*2

    #後縁総重量計算
    sum_rear_edge  = rear_edge_area*max_code2*rectangle_span*rho_rear_edge
##    print sum_rear_edge
    buf = (front_edge_length*(max_code-min_code)*(1-spar_pos)/taper_span)
    sum_rear_edge += rear_edge_area*max_code2*taper_span*(1+buf*buf/2)*rho_rear_edge*2
##    print sum_rear_edge

    sum_box = (taper_scale[2]*2)*(1-spar_pos)/2*tape_rib2rib*plank_thick*rho_ib * 4

    sum_stringer  = rectangle_span*stringer_area*rho_balsa *stringer_num
    buf = (front_edge_length*(max_code-min_code)*spar_pos/taper_span)
    sum_stringer += taper_span*(1+buf*buf/2)*stringer_area*rho_balsa *stringer_num*2

    sum_bond = bond*(rectangle_span+taper_span*2)

    sum_film = 0.1

    DAE21= 0.12412  #直径 #%36%かな？T-MITの古いプログラムから引用
    DAE31= 0.10889 #[m]

    chokei = DAE21*np.pi
    bond = 0.200 #接着剤重量[kg/m2]
    rectbondweight = chokei * rib_thick * rect_rib_num
    tapebondweight = chokei * rib_thick * tape_rib_num

    #ハーフスパンなので2倍、安全率で2倍
    bondweight = (rectbondweight + tapebondweight)*2. * 2.
    sum_bond = bondweight


    weightTE = span * 0.008 * 0.050 * 1/2. *140
    sum_rear_edge = weightTE
    #print "Rib      %0.3f [kg]"% sum_rib
    #print "Plank    %0.3f [kg]" % sum_plank
    #print "TE       %0.3f [kg]" % sum_rear_edge
    #print "Box      %0.3f [kg]" % sum_box
    #print "Stringer %0.3f [kg]" % sum_stringer
    #print "Bonds    %0.3f [kg]" % sum_bond
    summation =float(sum_rib + sum_plank + sum_rear_edge + sum_box + sum_stringer + sum_bond)
    #print "Weight   %0.3f [kg]" % summation
    summation = 1.3 * summation
    #print "Weight(1.3)%0.3f [kg]" % summation

    return summation


def htail(span):
    # エレベーター重量推算

    #メインの変数
    max_code = 0.6        #矩形部での翼弦長[m]
    min_code = 0.4        #翼端での翼弦長[m]
    spar_pos = 0.2849     #桁位置
    rectangle_span = 0.7  #矩形部長さ[m]
    taper_span = 1.25     #テーパ部長さ[m]

    rectangle_span = span * 0.4
    taper_span = span * 0.6

    #定数

    ##最大翼弦を1とした時の各定数(DAE31を想定)
    airfoil_area = 0.057901602     #翼型の面積
    front_edge_length = 0.678751688#前縁部長さ[m]
    rib_cap_length = 0.985979359   #リブキャップ長さ[m]
    spar_cap_area = 0.000665142    #桁キャップ断面積[m2]
    rear_edge_area = 0.000786169   #後縁材断面積[m2]
    rear_cap_area = 0.001426112    #後縁キャップ断面積[m2]

    ##密度[kg/m3]
    rho_ib = 28         #スタイロ密度(IB)
    rho_balsa = 140     #バルサ密度
    rho_rear_edge = 140 #後縁材密度

    ##他計算用定数
    rib_area_coeff = 0.888734796 #リブ面積係数(翼型の面積のうちリブ材の占める割合)

    rib_thick = 0.006        #リブ厚み[m]
    plank_thick = 0.003      #プランク厚み[m]
    balsa_thick = 0.001      #バルサ厚み[m]
    stringer_area = 0.000009 #ストリンガー断面積[m2]
    stringer_num = 3         #ストリンガー本数[本]

    limit_rect_rib2rib = 0.26 #矩形部のリブ間隔をこれ以下にする[m]
    limit_tape_rib2rib = 0.25 #テーパ部リブ間隔をこれ以下にする[m]

    bond = 0.2#接着剤重量[kg/m]

    #計算用変数
    rib_area =  airfoil_area * max_code*max_code #矩形部の翼型面積
    rect_rib_num = rectangle_span//limit_rect_rib2rib +2  #リブ枚数(矩形部)[枚]
    tape_rib_num = taper_span//limit_tape_rib2rib       #リブ枚数(テーパ部)[枚/片翼](矩形含まず）
    rect_rib2rib = rectangle_span / (rect_rib_num-1) #矩形部のリブ間隔[m]
    tape_rib2rib = taper_span / tape_rib_num     #テーパ部リブ間隔[m]

    max_code2 = max_code * max_code

    #リブまわり重量計算
    ##各1枚(本)あたりの重量
    rib_weight = rib_area * rib_area_coeff * rib_thick * rho_ib #リブ本体
    rib_cap_weight = rib_cap_length*max_code * rib_thick * balsa_thick * rho_balsa #リブキャップ
    spar_cap_weight = spar_cap_area*max_code2 * balsa_thick * rho_balsa #桁キャップ
    rear_cap_weight = rear_cap_area*max_code2 * balsa_thick * rho_balsa #後縁キャップ

    ##テーパ部の積算係数
    def sum_taper_scale(num):
        length_scale = 0
        area_scale = 0
        for i in range(1, int(num)):
            taper_code = max_code-(max_code-min_code)*tape_rib2rib*i/taper_span
            length_scale += taper_code
            area_scale += (taper_code/max_code)*(taper_code/max_code)
        return [length_scale, area_scale, taper_code]

    taper_scale = sum_taper_scale(tape_rib_num)
    #print rect_rib_num,tape_rib_num

    ##リブまわりの総重量
    sum_rib  = rib_weight     *(rect_rib_num+taper_scale[1]*2) #リブ本体
    sum_rib += rib_cap_weight *(rect_rib_num+taper_scale[0]*2) #リブキャップ
    sum_rib += spar_cap_weight*(rect_rib_num+taper_scale[1]*2) #桁キャップ
    sum_rib += rear_cap_weight*(rect_rib_num+taper_scale[1]*2) #後縁キャップ
    #print sum_rib

    #プランク総重量計算
    sum_plank  = front_edge_length*max_code *rectangle_span*plank_thick*rho_ib
    sum_plank += front_edge_length*(max_code+min_code)/2*taper_span*plank_thick*rho_ib*2
    #print sum_plank

    #後縁総重量計算
    sum_rear_edge  = rear_edge_area*max_code2*rectangle_span*rho_rear_edge
    buf = (front_edge_length*(max_code-min_code)*(1-spar_pos)/taper_span)
    sum_rear_edge += rear_edge_area*max_code2*taper_span*(1+buf*buf/2)*rho_rear_edge*2
    #print sum_rear_edge

    #ボックス部総重量計算
    sum_box = (taper_scale[2]*2)*(1-spar_pos)/2*tape_rib2rib*plank_thick*rho_ib * 4
    #print sum_box

    #ストリンガー総重量計算
    sum_stringer  = rectangle_span*stringer_area*rho_balsa *stringer_num
    buf = (front_edge_length*(max_code-min_code)*spar_pos/taper_span)
    sum_stringer += taper_span*(1+buf*buf/2)*stringer_area*rho_balsa *stringer_num*2
    #print sum_stringer

    #接着剤総重量計算
    sum_bond = bond*(rectangle_span+taper_span*2)
    #print sum_bond

    #桁重量
    pipeweight = span * 0.1


    chokei = 0.03*np.pi
    bond = 0.200 #接着剤重量[kg/m2]
    rectbondweight = chokei * rib_thick * rect_rib_num
    tapebondweight = chokei * rib_thick * tape_rib_num

    #ハーフスパンなので2倍、安全率で2倍
    bondweight = (rectbondweight + tapebondweight)*2. * 2.
    sum_bond = bondweight

    #総計
    #print sum_rib + sum_plank + sum_rear_edge + sum_box + sum_stringer + sum_bond + pipeweight

    #print "----htail----"
    #print "Rib      %0.3f [kg]"% sum_rib
    #print "Plank    %0.3f [kg]" % sum_plank
    #print "TE       %0.3f [kg]" % sum_rear_edge
    #print "Box      %0.3f [kg]" % sum_box
    #print "Stringer %0.3f [kg]" % sum_stringer
    #print "Bonds    %0.3f [kg]" % sum_bond
    #print "pipe     %0.3f [kg]" % pipeweight
    summation =float(sum_rib + sum_plank + sum_rear_edge + sum_box + sum_stringer + sum_bond + pipeweight)
    #print "Weight   %0.3f [kg]" % summation
    summation = 1.3 * summation
    #print "Weight(1.3)%0.3f [kg]" % summation

    return summation


def vtail(span):
    # ラダー重量推算

    #メインの変数
    max_code = 0.82      #矩形部での翼弦長[m]
    min_code_u = 0.55    #翼端での翼弦長(上)[m]
    min_code_l = 0.69     #翼端での翼弦長(下)[m]
    spar_pos = 0.2849    #桁位置
    rectangle_span = 0.5 #矩形部長さ[m]
    taper_span_u = 1.71  #テーパ部長さ(上)[m]
    taper_span_l = 0.41  #テーパ部長さ(下)[m]

    rectangle_span = span * 0.2
    taper_span_u = span * 0.65
    taper_span_l = span * 0.15

    #定数

    ##最大翼弦を1とした時の各定数(sd8020を想定)
    airfoil_area = 0.067941686     #翼型の面積
    front_edge_length = 0.822595956#前縁部長さ[m]
    rib_cap_length = 0.901990298   #リブキャップ長さ[m]
    spar_cap_area = 0.001301126    #桁キャップ断面積[m2]
    rear_edge_area = 0.001122057   #後縁材断面積[m2]
    rear_cap_area = 0.001282682    #後縁キャップ断面積[m2]

    ##密度[kg/m3]
    rho_ib = 28         #スタイロ密度(IB)
    rho_balsa = 140     #バルサ密度
    rho_rear_edge = 140 #後縁材密度

    ##他計算用定数
    rib_area_coeff = 0.929466342#リブ面積係数(翼型の面積のうちリブ材の占める割合)

    rib_thick = 0.006        #リブ厚み[m]
    plank_thick = 0.003      #プランク厚み[m]
    balsa_thick = 0.001      #バルサ厚み[m]
    stringer_area = 0.000009 #ストリンガー断面積[m2]
    stringer_num = 3         #ストリンガー本数[本]

    limit_rect_rib2rib = 0.25 #矩形部のリブ間隔をこれ以下にする[m]
    limit_tape_rib2rib = 0.25 #テーパ部リブ間隔をこれ以下にする[m]

    bond = 0.2#接着剤重量[kg/m]

    #計算用変数
    rib_area =  airfoil_area * max_code*max_code #矩形部の翼型面積
    rect_rib_num = rectangle_span//limit_rect_rib2rib +2  #リブ枚数(矩形部)[枚]
    tape_rib_num_u = taper_span_u//limit_tape_rib2rib       #リブ枚数(テーパ部上)[枚/片翼](矩形含まず）
    tape_rib_num_l = taper_span_l//limit_tape_rib2rib       #リブ枚数(テーパ部下)[枚/片翼](矩形含まず）
    rect_rib2rib = rectangle_span / (rect_rib_num-1) #矩形部のリブ間隔[m]
    tape_rib2rib_u = taper_span_u / tape_rib_num_u     #テーパ部リブ間隔(上)[m]
    tape_rib2rib_l = taper_span_l / tape_rib_num_l     #テーパ部リブ間隔(下)[m]

    max_code2 = max_code * max_code

    #リブまわり重量計算
    ##各1枚(本)あたりの重量
    rib_weight = rib_area * rib_area_coeff * rib_thick * rho_ib #リブ本体
    rib_cap_weight = rib_cap_length*max_code * rib_thick * balsa_thick * rho_balsa #リブキャップ
    spar_cap_weight = spar_cap_area*max_code2 * balsa_thick * rho_balsa #桁キャップ
    rear_cap_weight = rear_cap_area*max_code2 * balsa_thick * rho_balsa #後縁キャップ

    ##テーパ部の積算係数
    def sum_taper_scale(num, min_code, tape_rib2rib, taper_span):
        length_scale = 0
        area_scale = 0
        taper_code = 0
        for i in range(1, int(num)):
            taper_code = max_code-(max_code-min_code)*tape_rib2rib*i/taper_span
            length_scale += taper_code
            area_scale += (taper_code/max_code)*(taper_code/max_code)
        return [length_scale, area_scale, taper_code]

    taper_scale_u = sum_taper_scale(tape_rib_num_u, min_code_u, tape_rib2rib_u, taper_span_u)
    taper_scale_l = sum_taper_scale(tape_rib_num_l, min_code_l, tape_rib2rib_l, taper_span_l)

    #print rect_rib_num,tape_rib_num

    ##リブまわりの総重量
    sum_rib  = rib_weight     *(rect_rib_num+taper_scale_u[1]+taper_scale_l[1]) #リブ本体
    sum_rib += rib_cap_weight *(rect_rib_num+taper_scale_u[0]+taper_scale_l[1]) #リブキャップ
    sum_rib += spar_cap_weight*(rect_rib_num+taper_scale_u[1]+taper_scale_l[1]) #桁キャップ
    sum_rib += rear_cap_weight*(rect_rib_num+taper_scale_u[1]+taper_scale_l[1]) #後縁キャップ
    #print sum_rib

    #プランク総重量計算
    sum_plank  = front_edge_length*max_code *rectangle_span*plank_thick*rho_ib
    sum_plank += front_edge_length*(max_code+min_code_u)/2*taper_span_u*plank_thick*rho_ib
    sum_plank += front_edge_length*(max_code+min_code_l)/2*taper_span_l*plank_thick*rho_ib
    #print sum_plank

    #後縁総重量計算
    sum_rear_edge  = rear_edge_area*max_code2*rectangle_span*rho_rear_edge
    buf = (front_edge_length*(max_code-min_code_u)*(1-spar_pos)/taper_span_u)
    sum_rear_edge += rear_edge_area*max_code2*taper_span_u*(1+buf*buf/2)*rho_rear_edge
    buf = (front_edge_length*(max_code-min_code_l)*(1-spar_pos)/taper_span_l)
    sum_rear_edge += rear_edge_area*max_code2*taper_span_l*(1+buf*buf/2)*rho_rear_edge
    #print sum_rear_edge

    #ボックス部総重量計算
    sum_box  = (taper_scale_u[2]*2)*(1-spar_pos)/2*tape_rib2rib_u*plank_thick*rho_ib * 2
    sum_box += (taper_scale_l[2]*2)*(1-spar_pos)/2*tape_rib2rib_l*plank_thick*rho_ib * 2
    #print sum_box

    #ストリンガー総重量計算
    sum_stringer  = rectangle_span*stringer_area*rho_balsa *stringer_num
    buf = (front_edge_length*(max_code-min_code_u)*spar_pos/taper_span_u)
    sum_stringer += taper_span_u*(1+buf*buf/2)*stringer_area*rho_balsa *stringer_num
    buf = (front_edge_length*(max_code-min_code_l)*spar_pos/taper_span_l)
    sum_stringer += taper_span_l*(1+buf*buf/2)*stringer_area*rho_balsa *stringer_num
    #print sum_stringer

    #接着剤総重量計算
    sum_bond = bond*(rectangle_span+taper_span_u+taper_span_l)
    #print sum_bond


    chokei = 0.03*np.pi
    bond = 0.200 #接着剤重量[kg/m2]
    rectbondweight = chokei * rib_thick * rect_rib_num
    tapebondweight = chokei * rib_thick * tape_rib_num_l
    tapebondweight += chokei * rib_thick * tape_rib_num_u

    #ハーフスパンなので2倍、安全率で2倍
    bondweight = (rectbondweight + tapebondweight)*2. * 2.
    sum_bond = bondweight

    #桁重量
    pipeweight = span * 0.1

    #総計

    #print sum_rib + sum_plank + sum_rear_edge + sum_box + sum_stringer + sum_bond + pipeweight
    """
    print "----vtail----"
    print "Rib      %0.3f [kg]"% sum_rib
    print "Plank    %0.3f [kg]" % sum_plank
    print "TE       %0.3f [kg]" % sum_rear_edge
    print "Box      %0.3f [kg]" % sum_box
    print "Stringer %0.3f [kg]" % sum_stringer
    print "Bonds    %0.3f [kg]" % sum_bond
    print "pipe     %0.3f [kg]" % pipeweight
    """
    summation =float(sum_rib + sum_plank + sum_rear_edge + sum_box + sum_stringer + sum_bond + pipeweight)
    #print "Weight   %0.3f [kg]" % summation
    summation = 1.3 * summation
    #print "Weight(1.3)%0.3f [kg]" % summation

    return summation

if __name__ == '__main__':
    #スパンを入れる
    wing2(33.54)
    htail(3)
    vtail(2)
