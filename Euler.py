import math
import ezdxf
import ezdxf.entities
import matplotlib.pyplot as plt
import numpy as np
import cvxpy as cp
import sympy


class Line:
    def __init__(self, e: ezdxf.entities.Line = None, start=None, end=None) -> None:
        self.dxftype = 'LINE'
        if e is not None:
            self.start = e.dxf.start
            self.end = e.dxf.end
        else:
            self.start = tuple(start)
            self.end = tuple(end)


class Arc:
    def __init__(self, e: ezdxf.entities.Arc) -> None:
        self.dxftype = 'ARC'
        self.center = e.dxf.center
        self.start_angle = e.dxf.start_angle
        self.end_angle = e.dxf.end_angle
        self.radius = e.dxf.radius
        self.start = e.start_point
        #[e.center[0]+e.radius*math.cos(math.pi*(e.start_angle/180)), e.center[1]+e.radius*math.sin(math.pi*(e.start_angle/180)), 0]
        self.end = e.end_point
        #[e.center[0]+e.radius*math.cos(math.pi*(e.end_angle/180)), e.center[1]+e.radius*math.sin(math.pi*(e.end_angle/180)), 0]


class Ellipse:
    def __init__(self, e: ezdxf.entities.Ellipse) -> None:
        self.dxftype = 'ELLIPSE'
        self.center = e.dxf.center
        self.major_axis = e.dxf.major_axis
        self.ratio = e.dxf.ratio
        self.start_param = e.dxf.start_param
        self.end_param = e.dxf.end_param
        self.start = e.start_point
        self.end = e.end_point
        self.minor_axis = e.minor_axis


class Euler:
    def __init__(self) -> None:
        self.pointList = []
        self.pointCount = 0
        self.pointDegree = []
        self.edges = []
        self.graph = None
        self.graphWeight = None
        self.edgeDouble = []
        self.floydPath = None
        self.xMax = 0
        self.yMax = 0
        self.inf = 10000
        self.delta = 1  # 重边偏移量
        self.delta2 = 1  # 立马返回时的延长量
        self.startPoint = None

    def getPointIndex(self, point):
        for i in range(len(self.pointList)):
            if self.isSamePoint(self.pointList[i], point):
                self.pointDegree[i] += 1
                return i
        i = len(self.pointList)
        self.pointList.append(point)
        self.pointDegree.append(1)
        return i

    def isSamePoint(self, p1, p2):
        for i in range(3):
            if abs(p1[i]-p2[i]) > 0.1:
                return False
        return True

    def optGraph(self, graph):
        graph2 = [[None]*self.pointCount for i in range(self.pointCount)]
        for i in range(self.pointCount):
            for j in range(self.pointCount):
                graph2[i][j] = graph[i][j]
        return graph2

    def weightGraph(self):
        graphWeight = [[-1]*self.pointCount for i in range(self.pointCount)]
        for i in range(self.pointCount):
            for j in range(self.pointCount):
                weight = -1
                e = self.graph[i][j]
                if e == None:
                    continue
                if e.dxftype == 'LINE':
                    weight = math.sqrt((e.start[0]-e.end[0]) **
                                       2+(e.start[1]-e.end[1])**2)
                else:
                    weight = self.inf
                graphWeight[i][j] = weight
        return graphWeight

    def floyd(self):
        g = self.graphWeight
        l = len(g)
        self.floydPath = [[i for i in range(self.pointCount)]
                          for i in range(self.pointCount)]
        for k in range(l):
            for i in range(l):
                for j in range(l):
                    if i == j:
                        continue
                    if g[i][k] != -1 and g[k][j] != -1:
                        if g[i][j] == -1 or g[i][j] > g[i][k]+g[k][j]:
                            g[i][j] = g[i][k]+g[k][j]
                            self.floydPath[i][j] = self.floydPath[i][k]

    def getCost(self, oddList: list):
        l = len(oddList)
        oddGraph = np.zeros((l, l))
        for i in range(l):
            for j in range(l):
                if i == j:
                    oddGraph[i][j] = self.inf
                else:
                    oddGraph[i][j] = self.graphWeight[oddList[i]][oddList[j]]
        return self.cvx2(oddGraph)

    def cvx2(self, c):
        """ 返回最优值 """
        n = len(c)
        x = cp.Variable((n, n), integer=True)
        obj = cp.Minimize(cp.sum(cp.multiply(c, x)))
        con = [0 <= x, x <= 1, cp.sum(x, axis=0, keepdims=True) == 1,
               cp.sum(x, axis=1, keepdims=True) == 1,
               x == cp.transpose(x)
               ]
        prob = cp.Problem(obj, con)
        prob.solve(solver='GLPK_MI')
        #print("最优值为", prob.value)
        return prob.value

    def cvx(self, c) -> set:
        """ 返回最优解 """
        n = len(c)
        x = cp.Variable((n, n), integer=True)
        obj = cp.Minimize(cp.sum(cp.multiply(c, x)))
        con = [0 <= x, x <= 1, cp.sum(x, axis=0, keepdims=True) == 1,
               cp.sum(x, axis=1, keepdims=True) == 1,
               x == cp.transpose(x)
               ]
        prob = cp.Problem(obj, con)
        prob.solve(solver='GLPK_MI')
        #print("最优值为", prob.value)
        # print("最优解为：", x.value)
        sol = x.value
        edges = set()
        i = 0
        for row in sol:
            j = 0
            for col in row:
                if col == 1:
                    #print(i, j)
                    edges.add((i, j))
                j += 1
            i += 1
        return edges

    def pointClosestToOrigin(self):
        """ 离原点最近的点的索引 """
        minDis = self.inf
        index = 0
        for i in range(self.pointCount):
            p = self.pointList[i]
            dis = math.sqrt(p[0]**2+p[1]**2)
            if dis < minDis:
                minDis = dis
                index = i
        return index

    def getStartEndIndex(self):
        """ 所有最优起止点坐标 """
        # 找出所有度为1的点
        singlePointList = []
        oddDegreeCount = 0
        for i in range(len(self.pointDegree)):
            if self.pointDegree[i] == 1:
                singlePointList.append(i)
            if self.pointDegree[i] % 2 == 1:
                oddDegreeCount += 1
        print('度为奇数的点有', oddDegreeCount, '个')
        # 枚举所有可能的起止点组合
        n = len(singlePointList)
        print('度为1的点有', n, '个')
        if n == 0 or n == 1:
            print('度为1的点少于两个, 选择离原点最近的点作为起点')
            i = self.pointClosestToOrigin()
            return [(i, i)]
        seList = []
        for i in range(n):
            for j in range(i+1, n):
                seList.append((singlePointList[i], singlePointList[j]))
        # 计算每种起止点组合的加边代价
        costList = []
        for se in seList:
            pd = list(self.pointDegree)
            pd[se[0]] = 2
            pd[se[1]] = 2
            oddListTemp = []
            for i in range(self.pointCount):
                if pd[i] % 2 == 1:
                    oddListTemp.append(i)
            costList.append(round(self.getCost(oddListTemp))/2)
        minCost = min(costList)
        print('加边最小代价：', minCost)
        seListBest = []
        for i in range(len(seList)):
            # print('起始点：', seList[i])
            # print('cost: ', costList[i])
            if costList[i] == minCost:
                seListBest.append(seList[i])
        return seListBest

    def addEdge(self, oddList):
        l = len(oddList)
        oddGraph = np.zeros((l, l))  # [[0]*l for i in range(l)]
        for i in range(l):
            for j in range(l):
                if i == j:
                    oddGraph[i][j] = self.inf
                else:
                    oddGraph[i][j] = self.graphWeight[oddList[i]][oddList[j]]
        edgeSet = self.cvx(oddGraph)
        pointSet = set()
        for edgeIndex in edgeSet:
            # if row_ind[i] > col_ind[i]:
            #     continue
            edge = (oddList[edgeIndex[0]], oddList[edgeIndex[1]])
            pointSet.add(edge[0])
            pointSet.add(edge[1])
            # print(edge)
            if edge in self.edges:
                self.edgeDouble.append(edge)
            else:
                r = edge[0]
                c = edge[1]
                temp = self.floydPath[r][c]
                self.edgeDouble.append((r, temp))
                while temp != c:
                    temp2 = temp
                    temp = self.floydPath[temp][c]
                    self.edgeDouble.append((temp2, temp))
        #print('edgeDouble size: ', len(edgeDouble))
        for edge in self.edgeDouble:
            if (edge[1], edge[0]) not in self.edgeDouble:
                print('非对称解')
                raise Exception('非对称解')
                break

    def dfs(self, graphCopy, start, path: list):
        for i in range(len(graphCopy[start])):
            if graphCopy[start][i] is not None:
                graphCopy[start][i] = None
                if (start, i) not in self.edgeDouble:
                    graphCopy[i][start] = None
                self.dfs(graphCopy, i, path)
                path.append((start, i))

    def getBestPath(self, seListBest):
        for se in seListBest:
            oddList = []
            pd = list(self.pointDegree)
            startPointIndex = se[0]
            endPointIndex = se[1]
            pd[startPointIndex] = 2
            pd[endPointIndex] = 2
            for i in range(len(pd)):
                if pd[i] % 2 == 1:
                    oddList.append(i)
            self.edgeDouble.clear()
            if len(oddList) != 0:
                self.addEdge(oddList)
            path = []
            graphCopy = self.optGraph(self.graph)
            # graphCopy[startPointIndex][endPointIndex] = e
            # graphCopy[endPointIndex][startPointIndex] = e
            self.dfs(graphCopy, startPointIndex, path)
            rightPath = True
            for i in range(len(path)-1):
                if path[i][0] != path[i+1][1]:
                    rightPath = False
                    break
            if rightPath:
                self.startPoint = self.pointList[se[0]]
                print('最佳起止点坐标：', self.pointList[se[0]], self.pointList[se[1]])
                print('需要加', len(self.edgeDouble)/2, '条边')
                path.reverse()
                return path

    def move(self, start, end, direction: bool):
        start = list(start)
        end = list(end)
        try:
            k = (end[1]-start[1])/(end[0]-start[0])
        except:
            k = 10000
        if direction:
            # 重边第一次出现，往一个方向移位
            d = self.delta/2
        else:
            # 重边第二次出现，往令一个方向移位
            d = -self.delta/2
        if abs(k) < 0.001:
            start[1] += d
            end[1] += d
        elif abs(k) > 1000:
            start[0] += d
            end[0] += d
        else:
            theta = math.atan(-1/k)
            start[0] += d*math.cos(theta)
            start[1] += d*math.sin(theta)
            end[0] += d*math.cos(theta)
            end[1] += d*math.sin(theta)
        return start, end

    def reshapeEntities(self, path: list) -> list:
        drawed = set()
        entityList = []
        """ for i in range(len(path)-1):
            edge1 = path[i]
            edge2 = path[i+1]
            if edge1[0] == edge2[1]:  # 立马返回
                start = self.pointList[edge1[0]]
                end = self.pointList[edge1[1]]
                d = self.delta2 if end[1] > start[1] or end[0] > start[0] else -self.delta2
                if end[0]-start[0] == 0:
                    x = end[0]
                    y = end[1] + d
                else:
                    k = (end[1]-start[1])/(end[0]-start[0])
                    angle = math.atan(k)
                    x = end[0] + d*math.cos(angle)
                    y = end[1] + d*math.sin(angle)
                p = (x, y, 0)
                # self.graph[edge1[0]][edge1[1]].end = p
                # self.graph[edge2[0]][edge2[1]].start = p
                self.pointList[edge1[1]] = p """
        i = 0
        # for edge in path:
        while i < len(path):
            edge = path[i]
            e = self.graph[edge[0]][edge[1]]
            if e.dxftype == 'LINE':
                start = list(self.pointList[edge[0]])
                end = list(self.pointList[edge[1]])
                """ if edge not in drawed:
                    drawed.add(edge)
                    drawed.add((edge[1], edge[0]))
                else:
                    # 重边移位
                    try:
                        k = (end[1]-start[1])/(end[0]-start[0])
                    except:
                        k = 10000
                    if abs(k) < 0.001:
                        start[1] += self.delta
                        end[1] += self.delta
                    elif abs(k) > 1000:
                        start[0] += self.delta
                        end[0] += self.delta
                    else:
                        theta = math.atan(-1/k)
                        start[0] += self.delta*math.cos(theta)
                        start[1] += self.delta*math.sin(theta)
                        end[0] += self.delta*math.cos(theta)
                        end[1] += self.delta*math.sin(theta) """
                if edge in self.edgeDouble:
                    # 属于重边
                    if i < len(path)-1 and edge[0] == path[i+1][1]:  # 立马返回
                        xStart = round(start[0], 3)
                        yStart = round(start[1], 3)
                        xEnd = round(end[0], 3)
                        yEnd = round(end[1], 3)
                        #d = self.delta2 if end[1] > start[1] or end[0] > start[0] else -self.delta2
                        d = self.delta2 if xEnd > xStart or (
                            xEnd == xStart and yEnd > yStart) else -self.delta2
                        if xEnd-xStart == 0:
                            x = end[0]
                            y = end[1] + d
                        elif yEnd-yStart == 0:
                            x = end[0]+d
                            y = end[1]
                        else:
                            k = (end[1]-start[1])/(end[0]-start[0])
                            angle = math.atan(k)
                            x = end[0] + d*math.cos(angle)
                            y = end[1] + d*math.sin(angle)
                        p = [x, y, 0]  # 延长后的终点
                        start1, end1 = self.move(start, p, True)
                        start2, end2 = self.move(p, start, False)
                        e1 = Line(start=start, end=start1)
                        e2 = Line(start=start1, end=end1)
                        e3 = Line(start=end1, end=start2)
                        e4 = Line(start=start2, end=end2)
                        e5 = Line(start=end2, end=start)
                        entityList.extend([e1, e2, e3, e4, e5])
                        i += 1
                    else:
                        # 不是立马返回的情况
                        if edge not in drawed:
                            # 重边第一次出现，往一个方向移位
                            start1, end1 = self.move(start, end, True)
                            drawed.add(edge)
                            drawed.add((edge[1], edge[0]))
                        else:
                            # 重边第二次出现，往令一个方向移位
                            start1, end1 = self.move(start, end, False)
                        e1 = Line(start=start, end=start1)
                        e2 = Line(start=start1, end=end1)
                        e3 = Line(start=end1, end=end)
                        entityList.extend([e1, e2, e3])
                else:
                    # 不是重边，直接加入
                    e = Line(start=start, end=end)
                    entityList.append(e)
            else:
                # 不是直线，直接加入
                entityList.append(e)
            i += 1
        return entityList

    def getELLIPSEPoint(self, e):
        xa = e.major_axis[0]
        ya = e.major_axis[1]
        xb = e.minor_axis[0]
        yb = e.minor_axis[1]
        sn = (xa*yb-ya*xb)
        c = (xa**2+ya**2)**0.5
        d = c*e.ratio
        t = math.atan(ya/xa)
        spx = round(e.start[0], 4)
        spy = round(e.start[1], 4)
        epx = round(e.end[0], 4)
        epy = round(e.end[1], 4)
        x0, y0 = round(e.center[0], 4), round(e.center[1], 4)
        eangle1 = sympy.Symbol('eangle1')
        re = sympy.solve((x0 + c*sympy.cos(t)*sympy.cos(eangle1) -
                         d*sympy.sin(t)*sympy.sin(eangle1))-spx, eangle1)
        # print(re)
        eangle1 = round(re[0], 4)
        if abs((y0 + c*sympy.sin(t)*sympy.cos(eangle1)+d*sympy.cos(t)*sympy.sin(eangle1))-spy) < 0.001:
            pass
        else:
            eangle1 = round(re[1], 4)
        sa = eangle1/math.pi*180
        eangleend = sympy.Symbol('eangleend')
        re = sympy.solve((y0 + c*sympy.sin(t)*sympy.cos(eangleend) +
                         d*sympy.cos(t)*sympy.sin(eangleend))-epy, eangleend)
        eangleend = round(re[0], 4)
        if abs((x0 + c*sympy.cos(t)*sympy.cos(eangleend)-d*sympy.sin(t)*sympy.sin(eangleend))-epx) < 0.001:
            pass
        else:
            eangleend = round(re[1], 4)
        ea = eangleend/math.pi*180
        a = 2  # 分割度数
        xList = []
        yList = []
        if ea < sa:
            ea = ea+360

        if sn < 0:
            sa += 360
            sa, ea = ea, sa

        for angle in range(int(sa), int(ea)+a, a):
            angle = math.pi*(angle/180)  # 角度转弧度
            c1 = e.center[0] + c * \
                math.cos(t)*math.cos(angle)-d*math.sin(t)*math.sin(angle)
            c2 = e.center[1] + c * \
                math.sin(t)*math.cos(angle)+d*math.cos(t)*math.sin(angle)
            xList.append(c1)
            yList.append(c2)
        return xList, yList

    def plot(self, entities, title):
        plt.rcParams['font.sans-serif'] = ['FangSong']
        plt.plot()
        plt.title(title)
        plt.xlim(0, self.xMax+10)
        plt.ylim(0, self.yMax+10)
        l = len(entities)
        for i in range(l):
            e = entities[i]
            if e.dxftype == 'LINE':
                x = [e.start[0], e.end[0]]
                y = [e.start[1], e.end[1]]
            elif e.dxftype == 'ARC':
                start_ang = math.pi * (e.start_angle/180)
                end_ang = math.pi*(e.end_angle/180)
                if start_ang > end_ang:
                    end_ang += 2*math.pi
                angle = np.linspace(start_ang, end_ang)
                x = e.center[0]+e.radius*np.cos(angle)
                y = e.center[1]+e.radius*np.sin(angle)
            elif e.dxftype == 'ELLIPSE':
                x, y = self.getELLIPSEPoint(e)
            plt.plot(x, y)
            plt.pause(0.1)
        plt.show()

    def plot2(self, entities, ax):
        ax.plot()
        # ax.set_title(title)
        ax.set_xlim(0, self.xMax+10)
        ax.set_ylim(0, self.yMax+10)
        l = len(entities)
        for i in range(l):
            e = entities[i]
            if e.dxftype == 'LINE':
                x = [e.start[0], e.end[0]]
                y = [e.start[1], e.end[1]]
            elif e.dxftype == 'ARC':
                start_ang = math.pi * (e.start_angle/180)
                end_ang = math.pi*(e.end_angle/180)
                if start_ang > end_ang:
                    end_ang += 2*math.pi
                angle = np.linspace(start_ang, end_ang)
                x = e.center[0]+e.radius*np.cos(angle)
                y = e.center[1]+e.radius*np.sin(angle)
            elif e.dxftype == 'ELLIPSE':
                x, y = self.getELLIPSEPoint(e)
            ax.plot(x, y)
            plt.pause(0.1)
        # plt.show()

    def readDxf(self, filename: str) -> list:
        dxf = ezdxf.readfile(filename)
        print('开始分析dxf文件...')
        size = len(dxf.entities)*2
        graphTemp = [[None]*size for i in range(size)]
        print(len(dxf.entities), '条边')
        entities = []
        for e in dxf.entities:
            # print(type(e))
            # print(e.dxftype, e.layer)
            if e.dxf.dxftype == 'LINE':
                e = Line(e)
            elif e.dxf.dxftype == 'ARC':
                e = Arc(e)
            elif e.dxf.dxftype == 'ELLIPSE':
                e = Ellipse(e)
            else:
                raise Exception('dxf中有不能识别的实体, 只能识别线段、圆弧、椭圆弧')
            self.xMax = max(self.xMax, e.start[0], e.end[0])
            self.yMax = max(self.yMax, e.start[1], e.end[2])
            i = self.getPointIndex(e.start)
            j = self.getPointIndex(e.end)
            graphTemp[i][j] = e
            graphTemp[j][i] = e
            self.edges.append((i, j))
            self.edges.append((j, i))
            entities.append(e)
            # print(i, j)
        plt.rcParams['font.sans-serif'] = ['FangSong']
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
        ax[0].set_title('原图')
        self.plot2(entities, ax[0])
        self.pointCount = len(self.pointList)
        print(self.pointCount, '个点')
        print('dxf文件分析完毕。')
        print('开始加边转为欧拉图...')
        self.graph = self.optGraph(graphTemp)
        self.graphWeight = self.weightGraph()
        self.floyd()
        print('寻找最优起点...')
        seListBest = self.getStartEndIndex()
        print('寻找最优路径...')
        path = self.getBestPath(seListBest)
        entities = self.reshapeEntities(path)
        print('寻找最优路径完毕。')
        ax[1].set_title('加边转为欧拉图')
        self.plot2(entities, ax[1])
        plt.show()
        return entities


if __name__ == '__main__':
    euler = Euler()
    euler.readDxf("005.dxf")
