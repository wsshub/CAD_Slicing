import math
import numpy as np
import Euler
import sympy
import ezdxf
from cmath import pi
import ezdxf
import sympy
from numpy import arctan


def printARCPoint(e, angle, xList, yList):
    # print(angle)
    global c1, c2
    angle = math.pi*(angle/180)  # 角度转弧度
    c1 = e.center[0] + e.radius*math.cos(angle)
    c2 = e.center[1] + e.radius*math.sin(angle)
    # print(c1,c2)
    xList.append(c1)
    yList.append(c2)


def printELLIPSEPoint(e, c, d, t, angle, xList, yList):
    # print(angle)
    global c1, c2
    angle = math.pi*(angle/180)  # 角度转弧度
    c1 = e.center[0] + c*math.cos(t)*math.cos(angle) - \
        d*math.sin(t)*math.sin(angle)
    c2 = e.center[1] + c*math.sin(t)*math.cos(angle) + \
        d*math.cos(t)*math.sin(angle)
    xList.append(c1)
    yList.append(c2)
    #print(c1, c2)


def getXY(filename: str, euler: Euler.Euler):
    dxf = ezdxf.readfile(filename)
    xList = []
    yList = []
    #print('x   y')
    i = 0
    a = 1  # 分割度数
    global c1, c2
    entities = euler.readDxf(filename)
    length = len(entities)
    c1 = euler.startPoint[0]
    c2 = euler.startPoint[1]
   # print(f'共有{length}条线')
    for e in entities:
        i = i+1
       # print(f'第{i}条线段')
        #print(e.dxftype, e.layer)
        if e.dxftype == 'LINE':
            # print(e.start[0], e.start[1])  # 线段起点
            xList.append(e.start[0])
            yList.append(e.start[1])
            c1 = e.end[0]
            c2 = e.end[1]
            if i == length:
               # print(e.end[0], e.end[1])  # 线段终点
                xList.append(e.end[0])
                yList.append(e.end[1])
        elif e.dxftype == 'ARC':
            #print(e.center, e.radius, e.start_angle, e.end_angle)
            angle1 = math.pi*(e.start_angle/180)
            # print(c1,c2,angle1)
            if abs(c1 - (e.center[0] + e.radius*math.cos(angle1))) <= 1 and abs(c2 - (e.center[1] + e.radius*math.sin(angle1))) <= 1:
                a = abs(a)
                if e.start_angle > e.end_angle:
                    for angle in np.arange(int(e.start_angle), 360, a):  # 圆弧按间隔1°分割成若干点
                        printARCPoint(e, angle, xList, yList)
                    for angle in np.arange(0, int(e.end_angle)+a, a):  # 圆弧按间隔1°分割成若干点
                        printARCPoint(e, angle, xList, yList)
                else:
                    # 圆弧按间隔1°分割成若干点
                    for angle in np.arange(int(e.start_angle), int(e.end_angle)+a, a):
                        printARCPoint(e, angle, xList, yList)
            else:
                a = -abs(a)
                if e.start_angle > e.end_angle:
                    for angle in np.arange(int(e.end_angle), 0, a):  # 圆弧按间隔1°分割成若干点
                        printARCPoint(e, angle, xList, yList)
                    for angle in np.arange(360, int(e.start_angle)+a, a):  # 圆弧按间隔1°分割成若干点
                        printARCPoint(e, angle, xList, yList)
                else:
                    # 圆弧按间隔1°分割成若干点
                    for angle in np.arange(int(e.end_angle), int(e.start_angle)+a, a):
                        printARCPoint(e, angle, xList, yList)
        elif e.dxftype == 'ELLIPSE':
            xa = e.major_axis[0]
            ya = e.major_axis[1]
            xb = e.minor_axis[0]
            yb = e.minor_axis[1]
            sn = (xa*yb-ya*xb)
            #print(xa, ya)
            c = (xa**2+ya**2)**0.5
            d = c*e.ratio
            # print('a=',c)
            # print('b=',d)
            t = arctan(ya/xa)

            # print('t=',t/pi*180)
            spx = round(e.start[0], 4)
            spy = round(e.start[1], 4)
            epx = round(e.end[0], 4)
            epy = round(e.end[1], 4)
            # print(spx,spy,type(spx))
            # print(epx,epy)
            x0, y0 = round(e.center[0], 4), round(e.center[1], 4)
            eangle1 = sympy.Symbol('eangle1')

            re = sympy.solve((x0 + c*sympy.cos(t)*sympy.cos(eangle1) -
                             d*sympy.sin(t)*sympy.sin(eangle1))-spx, eangle1)
            # print(re)

            eangle1 = round(re[0], 4)
            # print(eangle1)
            if abs((x0 + c*sympy.sin(t)*sympy.cos(eangle1)+d*sympy.cos(t)*sympy.sin(eangle1))-spy) < 0.001:
                # print((e.center[1] + c*sympy.sin(t)*sympy.cos(eangle1)+d*sympy.cos(t)*sympy.sin(eangle1))-spy)
                pass
            else:
                eangle1 = round(re[1], 4)
                # print((e.center[1] + c*sympy.sin(t)*sympy.cos(eangle1)+d*sympy.cos(t)*sympy.sin(eangle1))-spy)
                # print(eangle1)
            sa = eangle1/pi*180
            # print(sa)
            eangleend = sympy.Symbol('eangleend')
            re = sympy.solve((y0 + c*sympy.sin(t)*sympy.cos(eangleend) +
                             d*sympy.cos(t)*sympy.sin(eangleend))-epy, eangleend)
            # print(re)
            eangleend = round(re[0], 4)
            if abs((x0 + c*sympy.cos(t)*sympy.cos(eangleend)-d*sympy.sin(t)*sympy.sin(eangleend))-epx) < 0.001:
                # print((e.center[0] + c*sympy.cos(t)*sympy.cos(eangleend)-d*sympy.sin(t)*sympy.sin(eangleend))-epx)
                pass
            else:
                eangleend = round(re[1], 4)
                # print((e.center[1] + c*sympy.sin(t)*sympy.cos(eangleend)+d*sympy.cos(t)*sympy.sin(eangleend))-epy)
            ea = eangleend/pi*180

            if ea < sa:
                ea = ea+360
            if sn < 0:
                sa += 360
                if abs(c1 - spx) <= 0.1 and abs(c2 - spy) <= 0.1:
                    a = -abs(a)
                    for angle in range(int(sa), int(ea)+a, a):
                        printELLIPSEPoint(e, c, d, t, angle, xList, yList)
                else:
                    a = abs(a)
                    sa, ea = ea, sa
                    for angle in range(int(sa), int(ea)+a, a):
                        printELLIPSEPoint(e, c, d, t, angle, xList, yList)
            else:
                if abs(c1 - spx) <= 0.1 and abs(c2 - spy) <= 0.1:
                    a = abs(a)
                    for angle in range(int(sa), int(ea)+a, a):
                        printELLIPSEPoint(e, c, d, t, angle, xList, yList)
                else:
                    sa, ea = ea, sa
                    a = -abs(a)
                    for angle in range(int(sa), int(ea)+a, a):
                        printELLIPSEPoint(e, c, d, t, angle, xList, yList)
    # print('*****************')
    # for i in range(len(xList)):
    #     print(xList[i], yList[i])
    return xList, yList


def writeTxt(x, y, z, e):
    line = f'G1 X{x} Y{y} Z{z} E{int(e*1000000)/1000000} F{100}'
    print(line)


def start(args):
    euler = Euler.Euler()
    euler.delta = args['delta']
    euler.delta2 = args['delta2']
    xList, yList = getXY(args['path'], euler)
    print('开始生成gcode...')
    # file = 'E:\弯曲0.25-1.xlsx'  # 文件路径
    k = args['k']  # 0.2   # 倍率
    step = args['step']  # 0.25  # z 每次增加量
    zMax = args['zMax']  # 1    # z 最大值
    xBack = args['xBack']  # 40        # x 退刀
    yBack = args['xBack']  # 100        # y 退刀
    zBack = args['xBack']  # 140
    #excel = pd.read_excel(file, header=None)
    # print(excel)
    #height, width = excel.shape
    height = len(xList)
    x, y, z = 0, 0, 0
    e = 0
    i = 0
    down = True
    turn = False
    first = True
    gcode = f'D:\\{args["name"]}.gcode'
    txt = open(gcode, 'w')
    txt.write('G28 X0 Y0\n\
    G28 Z0\n\
    G1 X10 Y60 F400\n\
    G1 Z3 F400\n\
    G92 E0\n\
    G1 E20\n\
    G92 E0\n\
    G92 X0 Y0\n\
    G1 Z0 E0 F100\n')
    while(z <= zMax):
        if(i == height):
            z = z+step
            down = False
            turn = True
            i = i-1
            continue
        elif(i == -1):
            z = z+step
            down = True
            turn = True
            i = i+1
            continue
        #list = excel.loc[i].values
        # x2, y2 = str(xList[i]),str(yList[i])#list[1], list[2]
        x2 = xList[i]  # float(x2[1:len(x2)])
        y2 = yList[i]  # float(y2[1: len(y2)])
        #z2 = float(z2[1: len(z2)])
        # print(f'x :{x}  y: {y}  z: {z}')
        # print(f'x2:{x2} y2:{y2} z2:{z2}')
        if(first):
            x = x2
            y = y2
            first = False
        if((i == 0 or i == height-1) and turn):
            e = (abs((x2-x)**2+(y2-y)**2+(step)**2))**0.5*k+e
            turn = False
        else:
            e = (abs((x2-x)**2+(y2-y)**2))**0.5*k+e
        x, y = x2, y2
        #writeTxt(x, y, z, e)
        txt.write(f'G1 X{x} Y{y} Z{z} E{int(e*1000000)/1000000} F{100}\n')
        if(down):
            i = i+1
        else:
            i = i-1
    txt.write(f'G1 X{xBack} Y{yBack} Z{zBack} E{int(e*1000000)/1000000} F100\n\
    M104 T0 S0\n\
    M140 S0\n\
    ')
    txt.close()
    print('生成完毕，保存到', gcode)
