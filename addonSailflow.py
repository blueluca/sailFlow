bl_info = {
    "name": "Create Sailprofile",
    "description": "Creates a profile for a sail",
    "author": "blueluca",
    "version": (0, 2, 4),
    "blender": (2, 7, 1),
    "api": 33411,  # Not certain on the API version
    "location": "View3D > Add > Mesh > Airfoil",
    "warning": "",
    "category": "Add Mesh"}

import imp
from math import *
import math
import sys
import time
import bpy
import re

from bpy.props import IntProperty, EnumProperty, BoolProperty, StringProperty
from bpy.types import Operator
from fpdf import FPDF
from mathutils import Vector, Euler, geometry
import mathutils
from mathutils.geometry import interpolate_bezier

# import pydevd


# PYDEV_SOURCE_DIR = 'C:/eclipse/plugins/org.python.pydev_3.4.1.201403181715/pysrc'

# if sys.path.count(PYDEV_SOURCE_DIR) < 1:
#   sys.path.append(PYDEV_SOURCE_DIR)


custom_profile = None

Vxs = []
FlatVxs = []
Flattened = []
VxFlat = []
flatObj = None
bpy.selection = []


def extractPerimeterEdges(m):
    el = []
    perifEdgeList = []
    for p in m.data.polygons:
        for e in p.edge_keys:
            if e[0] > e[1]:
                el.append((e[1], e[0]))
            else:
                el.append(e)
    while el:
        e = el.pop()
        if e in el:
            el.remove(e)
        else:
            perifEdgeList.append(e)
    return perifEdgeList


def profile(x, m, p):
    # =====================================================================================
    # The profile calculator
    #
    # inputs
    #   x  : point in the range 0 to 1 to calculate
    #   m : camber percentage
    #   p : camber position percentage
    # =====================================================================================
    if (x < 0):
        return 0
    if (x > 1):
        return 0
    if (x < p):
        y = (m / (p * p)) * (2 * p * x - x * x);
    else:
        y = (m / ((1 - p) * (1 - p))) * ((1 - 2 * p) + 2 * p * x - x * x);
    # debug(("Profile return Z="+str(y))
    return y


def curveProfile(x, vxs):
    for i in vxs:
        if i[0] > x:
            return i[1]
    return 0


class MyPoly:
    l1 = 0  # p1 - p3
    l2 = 0  # p1 - p2
    l3 = 0  # p2 - p3

    def clockWise(self, co1, co2, co3):
        A = (co2.x - co1.x) * (co2.y + co1.y)
        A = A + (co3.x - co2.x) * (co3.y + co2.y)
        A = A + (co1.x - co3.x) * (co1.y + co3.y)
        return A > 0

    def doubleCircle(self, d, R, r):
        x = (d ** 2 - r ** 2 + R ** 2) / (2 * d)
        y = 0.5 * (1 / d) * sqrt(4 * (d ** 2) * (R ** 2) - (d ** 2 - r ** 2 + R ** 2) ** 2)
        return x, y

    def __init__(self, poly, idx=-1):
        global Vxs

        a1 = poly.vertices[0]
        a2 = poly.vertices[1]
        a3 = poly.vertices[2]
        co1 = Vxs[a1].co
        co2 = Vxs[a2].co
        co3 = Vxs[a3].co

        #    print("Init ",a1,"-",co1,a2,"-",co2,a3,"-",co3)
        if self.clockWise(co1, co2, co3):
            self.p1 = a1
            self.p2 = a2
            self.p3 = a3
        elif self.clockWise(co1, co3, co2):
            self.p1 = a1
            self.p2 = a3
            self.p3 = a2
        elif self.clockWise(co2, co1, co3):
            self.p1 = a2
            self.p2 = a1
            self.p3 = a3
        elif self.clockWise(co2, co3, co1):
            self.p1 = a2
            self.p2 = a3
            self.p3 = a1
        elif self.clockWise(co3, co1, co2):
            self.p1 = a3
            self.p2 = a1
            self.p3 = a2
        else:
            self.p1 = a3
            self.p2 = a2
            self.p3 = a1

        self.l1 = abs((Vxs[self.p1].co - Vxs[self.p3].co).length)
        self.l2 = abs((Vxs[self.p1].co - Vxs[self.p2].co).length)
        self.l3 = abs((Vxs[self.p2].co - Vxs[self.p3].co).length)
        self.idx = idx
        self.vertices = [self.p1, self.p2, self.p3]

    #    print("Init end ",self)

    def __iter__(self):
        return self

    def __str__(self):
        return "MyPoly [" + str(self.idx) + "] " + str(self.p1) + "," + str(self.p2) + "," + str(
            self.p3)  # +" L:"+str(self.l1)+","+str(self.l2)+","+str(self.l3)

    def x(self, n):
        if n == 1:
            return Vxs[self.p1].co.x
        elif n == 2:
            return Vxs[self.p2].co.x
        elif n == 3:
            return Vxs[self.p3].co.x
        return 0

    def y(self, n):
        if n == 1:
            return Vxs[self.p1].co.y
        elif n == 2:
            return Vxs[self.p2].co.y
        elif n == 3:
            return Vxs[self.p3].co.y
        raise 0

    def z(self, n):
        if n == 1:
            return Vxs[self.p1].co.z
        elif n == 2:
            return Vxs[self.p2].co.z
        elif n == 3:
            return Vxs[self.p3].co.z
        raise 0

    def setx(self, n, val):
        if n == 1:
            Vxs[self.p1].co.x = val
        elif n == 2:
            Vxs[self.p2].co.x = val
        elif n == 3:
            Vxs[self.p3].co.x = val

    def sety(self, n, val):
        if n == 1:
            Vxs[self.p1].co.y = val
        elif n == 2:
            Vxs[self.p2].co.y = val
        elif n == 3:
            Vxs[self.p3].co.y = val

    def zeroZ(self):
        global VxFlat

        Vxs[self.p1].co.z = 0
        Vxs[self.p2].co.z = 0
        Vxs[self.p3].co.z = 0
        if self.p1 not in VxFlat:
            VxFlat.append(self.p1)
        if self.p2 not in VxFlat:
            VxFlat.append(self.p2)
        if self.p3 not in VxFlat:
            VxFlat.append(self.p3)

    def co(self, n):
        if n == 1:
            return Vxs[self.p1].co
        if n == 2:
            return Vxs[self.p2].co
        return Vxs[self.p3].co

    def flatAllVertices(self):
        global Vxs

        #   print("    Flat start",self)
        self.zeroZ()
        if self.x(1) < self.x(2):
            d = self.l2
            R = self.l1
            r = self.l3
            x, y = self.doubleCircle(d, R, r)
            self.setx(3, x)
            self.sety(3, -y)
            self.setx(1, 0)
            self.sety(1, 0)
            self.setx(2, self.l2)
            self.sety(2, 0)
        else:
            d = self.l2
            R = self.l3
            r = self.l1
            x, y = self.doubleCircle(d, R, r)
            self.setx(3, x)
            self.sety(3, y)
            self.setx(2, 0)
            self.sety(2, 0)
            self.setx(1, self.l2)
            self.sety(1, 0)

            #   print("    Flat end ["+str(self.idx)+ "] "+str(Vxs[self.p1].co)+str(Vxs[self.p2].co)+str(Vxs[self.p3].co))
        # debug("    L1: old="+str(self.l1)+" new=" +str((Vxs[self.p1].co - Vxs[self.p3].co).length))
        # debug("    L2: old="+str(self.l2)+" new="+str((Vxs[self.p1].co - Vxs[self.p2].co).length))
        # debug("    L3: old="+str(self.l3)+" new="+str((Vxs[self.p2].co - Vxs[self.p3].co).length))

        if not self.clockWise(self.co(1), self.co(2), self.co(3)):
            print("Error clockwise at:", self)

    def edgeLength(self, a, b):
        if a in [1, 2] and b in [1, 2]:
            return self.l2
        if a in [2, 3] and b in [2, 3]:
            return self.l3
        if a in [1, 3] and b in [1, 3]:
            return self.l1

    def NEWsetToFlat(self, a, b, flatted):
        xa = self.x(a)
        ya = self.y(a)
        xb = self.x(b)
        yb = self.y(b)
        print(self, end='')
        print("SetToFlat (a=%d,b=%d,f=%d) xa=%2.2f,ya=%2.2f,xb=%2.2f,yb=%2.2f" % (a, b, flatted, xa, ya, xb, yb))
        d = self.edgeLength(a, b)
        ra = self.edgeLength(a, flatted)
        rb = self.edgeLength(b, flatted)
        print("          d=%2.2f ra=%2.2f rb=%2.2f" % (d, ra, rb))
        k = 0.25 * math.sqrt((((ra + rb) * (ra + rb) - (d * d)) * (d * d - (ra - rb) * (ra - rb))))
        x1 = 0.5 * (xb + xa) + 0.5 * (xb - xa) * (ra * ra - rb * rb) / (d * d) + 2 * (yb - ya) * k / (d * d)
        x2 = 0.5 * (xb + xa) + 0.5 * (xb - xa) * (ra * ra - rb * rb) / (d * d) - 2 * (yb - ya) * k / (d * d)
        y1 = 0.5 * (yb + ya) + 0.5 * (yb - ya) * (ra * ra - rb * rb) / (d * d) - 2 * (xb - xa) * k / (d * d)
        y2 = 0.5 * (yb + ya) + 0.5 * (yb - ya) * (ra * ra - rb * rb) / (d * d) + 2 * (xb - xa) * k / (d * d)
        print("          x1=%2.2f y1=%2.2f x2=%2.2f y2=%2.2f" % (x1, y1, x2, y2))

        va = Vector([xa, ya])
        vb = Vector([xb, yb])
        vc1 = Vector([x1, y1])
        if self.clockWise(va, vb, vc1):
            self.setx(flatted, x1)
            self.sety(flatted, y1)
        else:
            self.setx(flatted, x2)
            self.sety(flatted, y2)

            # =======================================================================
            # print("  p1-p3 diff %2.2f perc %2.2f"%(self.l1- (Vxs[self.p1].co - Vxs[self.p3].co).length,
            #                                        (self.l1- (Vxs[self.p1].co - Vxs[self.p3].co).length)/self.l1),end='')
            # print("  p1-p2 diff %2.2f perc %2.2f"%(self.l2 - (Vxs[self.p1].co - Vxs[self.p2].co).length,
            #                                        (self.l2 - (Vxs[self.p1].co - Vxs[self.p2].co).length)/self.l2),end='')
            # print("  p2-p3 diff %2.2f perc %2.2f"%(self.l3 - (Vxs[self.p2].co - Vxs[self.p3].co).length,
            #                                        (self.l3 - (Vxs[self.p2].co - Vxs[self.p3].co).length)/self.l3))
            #
            # =======================================================================

    def setToFlat(self, master, slave, flatted):
        if self.x(master) > self.x(slave):
            x0 = self.x(slave)
            y0 = self.y(slave)
            x1 = self.x(master)
            y1 = self.y(master)
            d = self.edgeLength(master, slave)
            R = self.edgeLength(slave, flatted)
            r = self.edgeLength(master, flatted)
        else:
            x0 = self.x(master)
            y0 = self.y(master)
            x1 = self.x(slave)
            y1 = self.y(slave)
            d = self.edgeLength(master, slave)
            R = self.edgeLength(master, flatted)
            r = self.edgeLength(slave, flatted)

        x, y = self.doubleCircle(d, R, r)
        if (x1 - x0) < 0.000000001:
            k = 100000000
        else:
            k = (y1 - y0) / (x1 - x0)
        teta = atan(k)
        dx1 = x * cos(teta)
        dy1 = x * sin(teta)
        dx2 = y * cos(pi / 2 - abs(teta))
        dy2 = abs(y * sin(pi / 2 - abs(teta)))

        if self.x(master) > self.x(slave):
            if k > 0:
                self.setx(flatted, x0 + dx1 - dx2)
                self.sety(flatted, y0 + dy1 + dy2)
            else:
                self.setx(flatted, x0 + dx1 + dx2)
                self.sety(flatted, y0 + dy1 + dy2)
        else:
            if k > 0:
                self.setx(flatted, x0 + dx1 + dx2)
                self.sety(flatted, y0 + dy1 - dy2)
            else:
                self.setx(flatted, x0 + dx1 - dx2)
                self.sety(flatted, y0 + dy1 - dy2)

    def flatten(self):
        global Vxs
        global VxFlat

        if (self.p1 in VxFlat) and (self.p2 in VxFlat) and (self.p3 in VxFlat):
            #            print(self.idx," already flat ")
            return 0
        elif (((self.p1 not in VxFlat) and (self.p2 not in VxFlat)) or
                  ((self.p1 not in VxFlat) and (self.p3 not in VxFlat)) or
                  ((self.p2 not in VxFlat) and (self.p3 not in VxFlat))):
            self.flatAllVertices()
        else:
            if self.p3 not in VxFlat:
                self.setToFlat(1, 2, 3)
            elif self.p1 not in VxFlat:
                self.setToFlat(2, 3, 1)
            elif self.p2 not in VxFlat:
                self.setToFlat(3, 1, 2)

        self.zeroZ()

    def adjacentVx(self, n):
        if n == self.p1:
            return [[self.p2, self.l2], [self.p3, self.l1]]
        elif n == self.p2:
            return [[self.p1, self.l2], [self.p3, self.l3]]
        elif n == self.p3:
            return [[self.p1, self.l1], [self.p2, self.l3]]
        else:
            return None, None


class enerVertex:
    def __init__(self, i):
        self.idx = i
        self.adj = []
        self.energy = 0.0
        self.evx = 0

    #        print("Create enVX ",i)

    def addAdjacent(self, iandl):
        if not iandl in self.adj:
            self.adj.append(iandl)
        #            print("added adj ",iandl)

    def calcEnergy(self):
        global Vxs

        self.energy = 0
        #        print("..node ",self.idx)
        for v in self.adj:
            dl = (Vxs[v[0]].co - Vxs[self.idx].co).length - v[1]
            #            print("...with ",v[0]," dl=",dl," over l=",v[1])
            self.energy = self.energy + dl * dl / v[1]
        #            print("...calcEnergy [",self.idx,"-",v[0],"] dl=",dl," en=",self.energy)
        #        print("...calcEnery total ", self.idx,self.energy)
        return self.energy

    def isAdjacent(self, i):
        for a in self.adj:
            if a[0] == i:
                return True
        return False

    def calcDeltaEnergy(self, delta):
        global Vxs

        origEnergy = self.energy
        originalCo = Vxs[self.idx].co
        Vxs[self.idx].co.x = originalCo.x + delta
        self.overallpdx = self.calcEnergy()
        Vxs[self.idx].co.x = originalCo.x - delta
        self.overallmdx = self.calcEnergy()
        Vxs[self.idx].co.x = originalCo.x
        Vxs[self.idx].co.y = originalCo.y + delta
        self.overallpdy = self.calcEnergy()
        Vxs[self.idx].co.y = originalCo.y - delta
        self.overallmdy = self.calcEnergy()
        Vxs[self.idx].co = originalCo
        self.energy = origEnergy


class Flattener(bpy.types.Operator):
    bl_idname = "mesh.flattener"
    bl_label = "Flat surface"
    bl_description = "bla bla bla"

    def findAdjacentNotFlat(self, Polys, p):
        # debug("findAdj of :"+str(p),3)
        for testP in Polys:
            # debug("...checking "+str(testP),4)
            if (testP.p1 in p.vertices) and (testP.p2 in p.vertices):
                return testP
            elif (testP.p2 in p.vertices) and (testP.p3 in p.vertices):
                return testP
            elif (testP.p3 in p.vertices) and (testP.p1 in p.vertices):
                return testP
        return None

    def updateNodesMovToGain(self, e, delta):
        minEnergy = min(e.overallpdx, e.overallmdx, e.overallpdy, e.overallmdy)
        if e.energy > minEnergy:
            if e.overallpdx == minEnergy:
                self.nodesMovToGain[e.evx] = [delta, 0, e.energy - e.overallpdx]
            elif e.overallmdx == minEnergy:
                self.nodesMovToGain[e.evx] = [-delta, 0, e.energy - e.overallmdx]
            elif e.overallpdy == minEnergy:
                self.nodesMovToGain[e.evx] = [0, delta, e.energy - e.overallpdy]
            else:
                self.nodesMovToGain[e.evx] = [0, -delta, e.energy - e.overallmdy]
        else:
            self.nodesMovToGain[e.evx] = [0, 0, 0.0]

    def minimizeEnergy(self, F, maxDeformation, deltaDeformation):
        global Vxs

        vlist = []
        EVs = []
        self.nodesMovToGain = []

        print("=" * 80)
        print(" " * 40, "start", " " * 40)
        print("=" * 80)
        for p in F:
            vlist = vlist + p.vertices
        # remove duplicates
        vlist = list(set(vlist))
        # vlist is all the vertices of the flatten area, with no dup.
        # Build the database of vertices and their connected
        for v in vlist:
            # create a new energy vertix
            e = enerVertex(v)
            # get the triangles with the vertex vix in common
            C = []
            for p in F:
                if v in p.vertices:
                    C.append(p)
            # C now contains the list of polys with the vertix index vix
            for c in C:
                # get the other two vertices that belongs to the common poly
                l = c.adjacentVx(v)
                # if the first is not yet in the list added
                e.addAdjacent(l[0])
                # if the second is not yet in the list added
                e.addAdjacent(l[1])
            e.calcEnergy()
            e.evx = len(EVs)
            EVs.append(e)
        # end database build
        # print("end database build")

        MAXCOUNT = 5000
        delta = 0.1
        while delta >= 10 ** -deltaDeformation:
            print("---------------> NEW Loop minimize delta=", delta, "minimum decrease", maxDeformation)
            # Build the database
            for e in EVs:
                e.calcDeltaEnergy(delta)
            maxGain = maxDeformation
            maxIdx = None
            self.nodesMovToGain = [[0, 0, 0.0] for x in range(len(EVs))]
            # Search the node with the max decrease in energy
            # and populate the database of gains self.nodesMovToGain             
            for evIdx in range(len(EVs)):
                e = EVs[evIdx]
                if e.energy > maxDeformation:
                    self.updateNodesMovToGain(e, delta)
                    # if the gain is higher than the threshold or the best record the node

                    # if self.nodesMovToGain[evIdx][2]>0:
                    # print("Node ",e.idx,"energy=",e.energy," dec energy=",self.nodesMovToGain[evIdx][2])

                    if self.nodesMovToGain[evIdx][2] > maxGain:
                        # print("..new max decr energy for",e.idx,"dec energy=",self.nodesMovToGain[evIdx][2])
                        maxGain = self.nodesMovToGain[evIdx][2]
                        maxIdx = evIdx
                        maxE = e
            # for evIdx in range(len(EVs)):
            # - movToGain is now updated for the delta and for all nodes
            # all nodes with suitable decrease of energy are in the database
            # for the others the entry is 0

            # Go recursiverly in the database applying all changes for a max loops
            # or until no further improvement can be found (MAXCOUNT or no maxIdx)
            for count in range(MAXCOUNT):
                # if there one node that can minimize energy    
                if maxIdx:
                    # print("Selected node ",maxE.idx,"energy=",maxE.energy," dec energy=",self.nodesMovToGain[maxE.evx][2])
                    # apply the change in position
                    Vxs[maxE.idx].co.x += self.nodesMovToGain[maxIdx][0]
                    Vxs[maxE.idx].co.y += self.nodesMovToGain[maxIdx][1]
                    # update the energy content and the deltas for the winnig vertex
                    maxE.calcEnergy()
                    maxE.calcDeltaEnergy(delta)
                    # update the database of nodeMovToGain for the winning vertex       
                    self.updateNodesMovToGain(maxE, delta)
                    # print("..updating database for adjacent")
                    for ev in EVs:
                        if maxE.isAdjacent(ev.idx):
                            # print("... update",ev.idx)
                            ev.calcEnergy()
                            if ev.energy > maxDeformation:
                                ev.calcDeltaEnergy(delta)
                                self.updateNodesMovToGain(ev, delta)
                            else:
                                self.nodesMovToGain[ev.evx] = [0.0, 0.0, 0.0]

                    maxIdx = None
                    maxGain = maxDeformation
                    # print("scanning database again")
                    for evIdx in range(len(EVs)):
                        nmtg = self.nodesMovToGain[evIdx]
                        e = EVs[evIdx]
                        # if nmtg[2] > 0:
                        # print("..node ",EVs[evIdx].idx,nmtg)
                        if nmtg[2] > maxGain:
                            maxGain = nmtg[2]
                            maxIdx = evIdx
                            maxE = EVs[maxIdx]
                else:
                    # no more improvement exit the loop
                    break
                    # if maxIdx
            # for count in range(MAXCOUNT):

            # decrease the movement step and repeat   
            delta = delta / 10

        maxRE = 0
        idx = 0
        for e in EVs:
            if e.energy > maxRE:
                maxRE = e.energy
                idx = e.idx
        print("exit minimize energy, max residual energy node=", idx, " energy=", maxRE)
        return maxRE

    def makeFlattened(self, obj, energyMinimizer, maxDeformation, deltaDeformation, polSeed):
        global Vxs
        global Flattened
        global VxFlat

        Adjacents = []
        ToBeFlattened = []
        Vxs = []
        VxFlat = []
        me = obj.data
        Vxs = me.vertices[:]
        for p in me.polygons:
            if p.select:
                ToBeFlattened.append(MyPoly(p, idx=p.index))

        if polSeed != -1:
            for poly in ToBeFlattened:
                if poly.idx == polSeed:
                    Adjacents.append(poly)
                    ToBeFlattened.remove(poly)
                    break

        while ToBeFlattened or Adjacents:
            if Adjacents:
                poly = Adjacents.pop(0)
            else:
                poly = ToBeFlattened.pop()
            # Collect all the adjacent triangle
            # and put it in Active collection "A"
            while 1:
                adjPoly = self.findAdjacentNotFlat(ToBeFlattened, poly)
                if adjPoly:
                    Adjacents.append(adjPoly)
                    ToBeFlattened.remove(adjPoly)
                else:
                    break
            poly.flatten()
            Flattened.append(poly)

        if energyMinimizer:
            return (self.minimizeEnergy(Flattened, maxDeformation, deltaDeformation))
        else:
            return 0

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and context.active_object.type == 'MESH'

    def execute(self, context):
        global Vxs
        global Flattened
        global flatObj

        sce = bpy.context.scene
        # Record the selected faces
        bpy.context.active_object.update_from_editmode()
        # Check if all selected faces are triangles
        # before start
        obj = bpy.context.active_object
        polys = obj.data.polygons
        for p in polys:
            if p.select and len(p.vertices) > 3:
                self.report({'ERROR'}, 'Not all faces are triangles STOPPED')
                return {'FINISHED'}

        if Flattened:
            del Flattened[:]

        if sce.sailflow_model.useSeed == False:
            sce.sailflow_model.resEnergy =self.makeFlattened(bpy.context.active_object,
                                                             sce.sailflow_model.energyMinimizer,
                                                             sce.sailflow_model.maxDeformation,
                                                             sce.sailflow_model.deltaDeformation,
                                                             -1)
        else:
            sce.sailflow_model.resEnergy = self.makeFlattened(bpy.context.active_object,
                                                              sce.sailflow_model.energyMinimizer,
                                                              sce.sailflow_model.maxDeformation,
                                                              sce.sailflow_model.deltaDeformation,
                                                              sce.sailflow_model.polSeed)
        # Complet list of vertices used
        vIdxList = []
        faces = []
        vertices = []
        for p in Flattened:
            vIdxList.append(p.p1)
            vIdxList.append(p.p2)
            vIdxList.append(p.p3)
        vIdxList.sort()
        # The following remove duplicates
        vIdxList = list(set(vIdxList))
        for p in Flattened:
            p.p1 = vIdxList.index(p.p1)
            p.p2 = vIdxList.index(p.p2)
            p.p3 = vIdxList.index(p.p3)
            faces.append([p.p3, p.p2, p.p1])
        for vidx in vIdxList:
            vertices.append([Vxs[vidx].co.x, Vxs[vidx].co.y, Vxs[vidx].co.z])

        fmesh = bpy.data.meshes.new("panel")
        fmesh.from_pydata(vertices, [], faces)
        # Update the displayed mesh with the new data
        fmesh.update()
        fobj = bpy.data.objects.new("panel", fmesh)
        #        fobj.matrix_world = bpy.context.active_object.matrix_world
        fobj.scale = bpy.context.active_object.scale
        scene = bpy.context.scene
        scene.objects.link(fobj)
        fobj.location = Vector([0,0,0])
        flatObj = fobj
        return {'FINISHED'}


class VIEW3D_PT_airprofile_parameters(bpy.types.Panel):
    bl_label = "Profile Parameters"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Sailflow Design"

    def draw(self, context):
        sce = bpy.context.scene
        layout = self.layout
        col = layout.column(align=True)
        col.prop(sce.sailflow_model, "t")

        if sce.sailflow_model.t == 'NACA':
            col = layout.column(align=True)
            col.prop(sce.sailflow_model, "m")
            col.prop(sce.sailflow_model, "p")

        elif sce.sailflow_model.t == "CUSTOM":
            col.operator("mesh.load_library")
            col = layout.column(align=True)

        elif sce.sailflow_model.t == "CURVE":
            col.operator("mesh.curve_aquire")
            col = layout.column(align=True)

        elif sce.sailflow_model.t == "DAT":
            col.operator("mesh.curve_load")
            col.prop(sce.sailflow_model, "resolution")
            col = layout.column(align=True)
            col.operator("mesh.loft")
            col.prop(sce.sailflow_model, "steps")
            col.prop(sce.sailflow_model, "spans")
            col = layout.column(align=True)

        if sce.sailflow_model.t != "DAT":
            row = col.row(align=True)
            row.prop(sce.sailflow_model, "twist")
            row.prop(sce.sailflow_model, "tw")
            row = col.row(align=True)
            row.prop(sce.sailflow_model, "ellipDis")
            row = col.row(align=True)
            row.prop(sce.sailflow_model, "ellipAmount")
            row = col.row(align=True)
            row.prop(sce.sailflow_model, "ellipCenter")
            row = col.row(align=True)
            row.prop(sce.sailflow_model, "shrink")

        #
        #       col = layout.column(align=True)
        if sce.sailflow_model.t != "DAT":
            col.label(text="Operation:")
            col.operator("mesh.airprof")

        if sce.sailflow_model.t == 'NACA':
            box = layout.box()
            angleIn = atan(2 * sce.sailflow_model.m / sce.sailflow_model.p) * 180 / 3.14159
            angleOut = -atan(2 * sce.sailflow_model.m / 100 / (sce.sailflow_model.p / 100 - 1)) * 180 / 3.14159
            box.label(text="Angle in  " + str(round(angleIn, 2)))
            box.label(text="Angle out " + str(round(angleOut, 2)))


class VIEW3D_PT_flattener_parameters(bpy.types.Panel):
    bl_label = "Flattener Parameters"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Sailflow Design"
    bpy.types.Object.obj_property = bpy.props.FloatProperty(name="ObjectProperty")

    def areaSelectedFaces(self, obj):
        if obj == None:
            return 0
        if obj.type != 'MESH':
            return 0
        pols = obj.data.polygons
        polsSel = [p for p in pols if p.select]

        area = 0.0
        for p in polsSel:
            area += p.area
        return area

    def draw(self, context):
        sce = bpy.context.scene
        layout = self.layout
        col = layout.column(align=True)

        col.prop(sce.sailflow_model, "energyMinimizer")
        col.prop(sce.sailflow_model, "deltaDeformation")
        col.prop(sce.sailflow_model, "maxDeformation")
        col.prop(sce.sailflow_model, "useSeed")
        col.prop(sce.sailflow_model, "polSeed")
        col.prop(sce.sailflow_model, "resEnergy")

        col = layout.column(align=True)
        col.operator("mesh.flattener")
        box = layout.box()

        areaS = self.areaSelectedFaces(bpy.context.active_object)
        areaF = self.areaSelectedFaces(flatObj)

        box.label(text="Area panel on sail m^2 =" + str(round(areaS, 5)))
        box.label(text="Area panel flatten m^2 =" + str(round(areaF, 5)))
        box.label(text="Total area diff   cm^2 =" + str(round((areaF - areaS) * 10 ** 4, 1)))


class VIEW3D_PT_airprofile_print(bpy.types.Panel):
    bl_label = "Printout Generation"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Sailflow Design"

    def draw(self, context):
        sce = bpy.context.scene
        layout = self.layout
        col = layout.column(align=True)

        col.prop(sce.sailflow_model, "paperFormat")
        col.prop(sce.sailflow_model, "freeText")
        col.prop(sce.sailflow_model, "multiPages")
        col.prop(sce.sailflow_model, "margin")
        col.prop(sce.sailflow_model, "overlap")

        if sce.sailflow_model.paperFormat == 'Other':
            col.prop(sce.sailflow_model, "paperWidth")
            col.prop(sce.sailflow_model, "paperHeight")

        col = layout.column(align=True)
        col.operator("mesh.print_pdf")
        #
        # Print in ASCII
        #
        col = layout.column(align=True)
        col.prop(sce.sailflow_model, "asciiDx")
        col.operator("mesh.print_ascii")


class VIEW3D_PT_analyse(bpy.types.Panel):
    bl_label = "Analyse sail"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Sailflow Design"

    def draw(self, context):
        sce = bpy.context.scene
        layout = self.layout
        col = layout.column(align=True)
        col.operator("mesh.colorit")


class colorIt(bpy.types.Operator):
    bl_idname = "mesh.colorit"
    bl_label = "Color selected sail"
    bl_description = "Give color depending depth"

    def execute(self, context):

        if context.mode in 'MESH_EDIT':
            bpy.ops.object.mode_set(mode='OBJECT')
        obj = context.selected_objects[0]
        scn = context.scene
        maxz = max([v.co.z for v in obj.data.vertices])
        minz = min([v.co.z for v in obj.data.vertices])
        mesh = obj.data
        if mesh.vertex_colors:
            vcol_layer = mesh.vertex_colors.active
        else:
            vcol_layer = mesh.vertex_colors.new()
        colors = []
        for i in range(1, 102):
            colors.append((float(i) / 100.0, 1.0 - float(i) / 100.0, 1.0 - float(i) / 100.0))
        print(colors)

        for poly in mesh.polygons:
            for loop_index in poly.loop_indices:
                loop_vert_index = mesh.loops[loop_index].vertex_index
                entry = (mesh.vertices[loop_vert_index].co.z - minz) / (maxz - minz)
                vcol_layer.data[loop_index].color = colors[int(entry * 100.0)]
        return {'FINISHED'}


class DATLoad(bpy.types.Operator):
    bl_idname = "mesh.curve_load"
    bl_label = "Load airfoil data as mesh"
    bl_description = "create mesh from airfoils coordinates"

    filepath = bpy.props.StringProperty(subtype="FILE_PATH")

    def createMesh(self, objname, Vert, Edges=[], Faces=[]):
        """Helper Function to Create Meshes"""
        me = bpy.data.meshes.new(objname)
        ob = bpy.data.objects.new(objname, me)
        bpy.context.scene.objects.link(ob)

        me.from_pydata(Vert, Edges, Faces)
        me.update(calc_edges=True)

    def CubicInterpolate(self, y0, y1, y2, y3, mu):
        mu2 = mu * mu
        a0 = -0.5 * y0 + 1.5 * y1 - 1.5 * y2 + 0.5 * y3
        a1 = y0 - 2.5 * y1 + 2 * y2 - 0.5 * y3
        a2 = -0.5 * y0 + 0.5 * y2
        a3 = y1
        return (a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3)

    def upload(self, Foil, Resolution=250, interp_method="l"):
        FF = open(Foil, 'r')
        data = FF.readlines()

        # Create the regexp for finding the data
        r = re.compile('[.\S]*\.[0-9]*')
        # Get FoilName from the file
        FoilName = Foil.strip()
        # Create the coordintes from the regler exp
        FoilCoords = [r.findall(x) for x in data[1:]]
        # Convert the strings to Floats
        RawPoints = [(float(x[0]), float(x[1])) for x in FoilCoords if len(x) == 2]
        # Ensure the First point is not the Point Count that some DAT files include

        if RawPoints[0][0] > 1: RawPoints.remove(RawPoints[0])
        """Process to divide the foildata to upper and lower sections"""
        FoilGrad = [(RawPoints[i][0] - RawPoints[i + 1][0]) for i in range(len(RawPoints) - 1)]

        for i in range(len(FoilGrad) - 1):
            if FoilGrad[i] >= 0. >= FoilGrad[i + 1]:
                if FoilGrad[i + 1] <= 0. <= FoilGrad[i + 2]:
                    splitloc = i + 1
                else:
                    splitloc = i
                break
            elif FoilGrad[i] <= 0. <= FoilGrad[i + 1]:
                if FoilGrad[i + 1] >= 0. >= FoilGrad[i + 2]:
                    splitloc = i + 1
                else:
                    splitloc = i
                break
        # Split the airfoil along chord
        upper = RawPoints[:splitloc + 1]
        x = []
        y = []
        for u in upper:
            x.append(u[0])
            y.append(u[1])

        der = (y[1] - y[0]) / (x[1] - x[0])
        yp = y[0] - der * (x[1] - x[0])
        x.insert(0, x[0] - (x[1] - x[0]))
        y.insert(0, yp)
        yp = y[0] - der * (x[1] - x[0])
        x.insert(0, x[0] - (x[1] - x[0]))
        y.insert(0, yp)
        # end
        l = len(y) - 1
        der = (y[l] - y[l - 1]) / (x[l] - x[l - 1])
        yp = y[l] + der * (x[l] - x[l - 1])
        x.append(x[l] + (x[l] - x[l - 1]))
        y.append(yp)
        l = len(y) - 1
        yp = y[l] + der * (x[l] - x[l - 1])
        x.append(x[l] + (x[l] - x[l - 1]))
        y.append(yp)

        iy = []
        ix = []
        if Resolution < len(x):
            Resolution = len(x)
        numpoints = int(Resolution / (len(x)-1))
        for idx in range(1, len(y) - 2):
            for pp in range(0, numpoints):
                ppx = float(pp) / numpoints
                iy.append(self.CubicInterpolate(y[idx - 1], y[idx], y[idx + 1], y[idx + 2], ppx))
                ix.append(x[idx] + (x[idx + 1] - x[idx]) * ppx)
        # idx = len(y)-1
        # for pp in range(0, numpoints):
        #     ppx = float(pp) / numpoints
        #     iy.append(self.CubicInterpolate(y[idx - 3], y[idx-2], y[idx-1], y[idx], ppx))
        # ix.append(x[idx] - (x[idx-1] - x[idx]) * ppx)

        return (ix, iy, RawPoints)

    def execute(self, context):
        global custom_profile
        print("File Path: ", self.filepath)

        sce = bpy.context.scene
        InterpX, InterpY, RawPoints = self.upload(self.filepath, Resolution=sce.sailflow_model.resolution)

        verts = [(InterpX[idx], 0, InterpY[idx]) for idx in range(0, len(InterpX) - 1)]
        for i,v in enumerate(verts):
            print (i,v)
        # if the point 0 is not there we must create it.
        #if (min(InterpX) > 0.0) and (min(InterpY) > 0.0):
        #   verts.insert(0, (0, 0, 0))

        edges = []
        for i in range(1, len(verts)):
            edges.append([i - 1, i])
        print(edges)
        self.createMesh("airfDAT", Edges=edges, Vert=verts)

        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class LoftDAT(bpy.types.Operator):
    bl_idname = "mesh.loft"
    bl_label = "Loft curves or meshes"
    bl_description = "Loft between mesh or curves"

    @staticmethod
    def cubic(p, t):
        return p[0] * (1.0 - t) ** 3.0 + 3.0 * p[1] * t * (1.0 - t) ** 2.0 \
               + 3.0 * p[2] * (t ** 2.0) * (1.0 - t) + p[3] * t ** 3.0

    def getbezpoints(self, spl, mt, seg=0):
        points = spl.bezier_points
        p0 = mt * points[seg].co
        p1 = mt * points[seg].handle_right
        p2 = mt * points[seg + 1].handle_left
        p3 = mt * points[seg + 1].co
        return p0, p1, p2, p3

    def getnurbspoints(self, spl, mw):
        pts = []
        ws = []
        for p in spl.points:
            v = Vector(p.co[0:3]) * mw
            pts.append(v)
            ws.append(p.weight)
        return pts, ws

    def knots(self, n, order, type=0):  # 0 uniform 1 endpoints 2 bezier

        kv = []

        t = n + order
        if type == 0:
            for i in range(0, t):
                kv.append(1.0 * i)

        elif type == 1:
            k = 0.0
            for i in range(1, t + 1):
                kv.append(k)
                if i >= order and i <= n:
                    k += 1.0
        elif type == 2:
            if order == 4:
                k = 0.34
                for a in range(0, t):
                    if a >= order and a <= n: k += 0.5
                    kv.append(floor(k))
                    k += 1.0 / 3.0

            elif order == 3:
                k = 0.6
                for a in range(0, t):
                    if a >= order and a <= n: k += 0.5
                    kv.append(floor(k))

        ##normalize the knot vector
        for i in range(0, len(kv)):
            kv[i] = kv[i] / kv[-1]

        return kv

    def B(self, i, k, t, knots):
        ret = 0
        if k > 0:
            n1 = (t - knots[i]) * self.B(i, k - 1, t, knots)
            d1 = knots[i + k] - knots[i]
            n2 = (knots[i + k + 1] - t) * self.B(i + 1, k - 1, t, knots)
            d2 = knots[i + k + 1] - knots[i + 1]
            if d1 > 0.0001 or d1 < -0.0001:
                a = n1 / d1
            else:
                a = 0
            if d2 > 0.0001 or d2 < -0.0001:
                b = n2 / d2
            else:
                b = 0
            ret = a + b
            # print "B i = %d, k = %d, ret = %g, a = %g, b = %g\n"%(i,k,ret,a,b)
        else:
            if knots[i] <= t and t <= knots[i + 1]:
                ret = 1
            else:
                ret = 0
        return ret

    def C(self, t, order, points, weights, knots):
        # c = Point([0,0,0])
        c = Vector()
        rational = 0
        i = 0
        while i < len(points):
            b = self.B(i, order, t, knots)
            p = points[i] * (b * weights[i])
            c = c + p
            rational = rational + b * weights[i]
            i = i + 1

        return c * (1.0 / rational)

    # Return the coordinate of a point at t percentage
    # of the entire line
    def calct(self, obj, t):

        if obj.type == 'CURVE':
            spl = None
            mw = obj.matrix_world
            if obj.data.splines.active is None:
                if len(obj.data.splines) > 0:
                    spl = obj.data.splines[0]
            else:
                spl = obj.data.splines.active

            if spl is None:
                return False

            if spl.type == "BEZIER":
                points = spl.bezier_points
                nsegs = len(points) - 1

                d = 1.0 / nsegs
                seg = int(t / d)
                t1 = t / d - int(t / d)

                if t == 1:
                    seg -= 1
                    t1 = 1.0

                p = self.getbezpoints(spl, mw, seg)

                coord = self.cubic(p,t1)
                return coord

            elif spl.type == "NURBS":
                data = self.getnurbspoints(spl, mw)
                pts = data[0]
                ws = data[1]
                order = spl.order_u
                n = len(pts)
                ctype = spl.use_endpoint_u
                kv = self.knots(n, order, ctype)

                coord = self.C(t, order - 1, pts, ws, kv)

                return coord
        elif obj.type == 'MESH':
            if t == 1:
                t = 0.999
            mw = obj.matrix_world
            vNum = len(obj.data.vertices)
            vidx = int(vNum * t)
            t1 = vNum * t - int(vNum * t)
            print(
                "CALCT: The object", obj.name, " has num vertices=", vNum, ' index calculated is=', vidx, "with a t=",
                t,
                "and t1=", t1)
            if vidx == 0:
                coord = mw * (
                    (obj.data.vertices[vidx + 1].co - obj.data.vertices[vidx].co) * t1 + obj.data.vertices[vidx].co)
            else:
                coord = mw * (
                    (obj.data.vertices[vidx].co - obj.data.vertices[vidx - 1].co) * t1 + obj.data.vertices[vidx - 1].co)
            print("CALCT: returning the coordinate", coord)
            return coord
        else:
            assert false

    def intc(self, objs, i, t, tr, tipo=3, tension=0.0, bias=0.0):

        ncurves = len(objs)

        # if 2 curves go to linear interpolation regardless the one you choose
        if ncurves < 3:
            return self.intl(objs, i, t, tr)
        else:

            # calculates the points to be interpolated on each curve
            if i == 0:
                p0 = self.calct(objs[i], t)
                p1 = p0
                p2 = self.calct(objs[i + 1], t)
                p3 = self.calct(objs[i + 2], t)
            else:
                if ncurves - 2 == i:
                    p0 = self.calct(objs[i - 1], t)
                    p1 = self.calct(objs[i], t)
                    p2 = self.calct(objs[i + 1], t)
                    p3 = p2
                else:
                    p0 = self.calct(objs[i - 1], t)
                    p1 = self.calct(objs[i], t)
                    p2 = self.calct(objs[i + 1], t)
                    p3 = self.calct(objs[i + 2], t)

        # calculates the interpolation between those points
        # i used methods from this page: http://paulbourke.net/miscellaneous/interpolation/

        if tipo == 0:
            # linear
            return self.intl(objs, i, t, tr)
        elif tipo == 1:
            # natural cubic
            t2 = tr * tr
            a0 = p3 - p2 - p0 + p1
            a1 = p0 - p1 - a0
            a2 = p2 - p0
            a3 = p1
            return a0 * tr * t2 + a1 * t2 + a2 * tr + a3
        elif tipo == 2:
            # catmull it seems to be working. ill leave it for now.
            t2 = tr * tr
            a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3
            a1 = p0 - 2.5 * p1 + 2 * p2 - 0.5 * p3
            a2 = -0.5 * p0 + 0.5 * p2
            a3 = p1
            return a0 * tr * tr + a1 * t2 + a2 * tr + a3

        elif tipo == 3:
            # hermite
            tr2 = tr * tr
            tr3 = tr2 * tr
            m0 = (p1 - p0) * (1 + bias) * (1 - tension) / 2
            m0 += (p2 - p1) * (1 - bias) * (1 - tension) / 2
            m1 = (p2 - p1) * (1 + bias) * (1 - tension) / 2
            m1 += (p3 - p2) * (1 - bias) * (1 - tension) / 2
            a0 = 2 * tr3 - 3 * tr2 + 1
            a1 = tr3 - 2 * tr2 + tr
            a2 = tr3 - tr2
            a3 = -2 * tr3 + 3 * tr2

            return a0 * p1 + a1 * m0 + a2 * m1 + a3 * p2

    def intl(self, objs, i, t, tr):
        p1 = self.calct(objs[i], t)
        p2 = self.calct(objs[i + 1], t)

        r = p1 + (p2 - p1) * tr

        return r

    def loft(self, objs, steps, spans, interpolation=0, tension=0.0, bias=0.5):
        verts = []

        # for each object
        for i in range(0, len(objs)):
            # For each step
            for j in range(0, steps + 1):
                # t = percentage of steps according to j
                t = 1.0 * j / steps
                # verts filled in with the coordinate at the
                # point t of the curve
                verts.append(self.calct(objs[i], t))

            temp2 = []
            if i < len(objs) - 1:
                for l in range(1, spans):
                    tr = 1.0 * l / spans
                    for k in range(0, steps + 1):
                        t = 1.0 * k / steps
                        if interpolation:
                            pos = self.intc(objs, i, t, tr, tipo=interpolation, tension=tension, bias=bias)
                        else:
                            pos = self.intl(objs, i, t, tr)

                        temp2.append(pos)
                verts.extend(temp2)
        return verts

    def execute(self, context):
        print("Called Loft")
        sce = context.scene
        objs = bpy.selection
        print(objs)
        spans = sce.sailflow_model.spans
        steps = sce.sailflow_model.steps

        intype = 2  # no interpolation

        verts = self.loft(objs, steps, spans, intype)

        nfaces = steps * spans * (len(objs) - 1)
        faces = []
        for i in range(0, nfaces):
            d = int(i / steps)
            # f = [i + d, i + d + 1, i + d + steps + 2, i + d + steps + 1]
            # inverts normals
            f = [i + d, i + d + steps + 1, i + d + steps + 2, i + d + 1]
            faces.append(f)

        me = bpy.data.meshes.new("Loft")
        me.from_pydata(verts, [], faces)
        me.update()
        newobj = bpy.data.objects.new("Loft", me)
        # newobj.data = me
        scn = bpy.context.scene
        scn.objects.link(newobj)
        scn.objects.active = newobj
        newobj.select = True
        bpy.ops.object.shade_smooth()

        return {'FINISHED'}


class CurveAquire(bpy.types.Operator):
    bl_idname = "mesh.curve_aquire"
    bl_label = "acquire Bezier or Mesh Curve"
    bl_description = "Generate the sail custom profile from curve"

    def get_points(self, sp, clean=True):

        knots = sp.bezier_points
        if len(knots) < 2:
            return

        # verts per segment
        r = sp.resolution_u + 1

        # segments in spline
        segments = len(knots)

        if not sp.use_cyclic_u:
            segments -= 1

        master_point_list = []
        for i in range(segments):
            inext = (i + 1) % len(knots)

            knot1 = knots[i].co
            handle1 = knots[i].handle_right
            handle2 = knots[inext].handle_left
            knot2 = knots[inext].co

            bezier = knot1, handle1, handle2, knot2, r
            points = interpolate_bezier(*bezier)
            master_point_list.extend(points)

        # some clean up to remove consecutive doubles, this could be smarter...
        if clean:
            old = master_point_list
            good = [v for i, v in enumerate(old[:-1]) if not old[i] == old[i + 1]]
            good.append(old[-1])
            return good

        return master_point_list

    def execute(self, context):
        print("called acquire")
        obj = bpy.context.active_object
        sce = bpy.context.scene

        if obj.type == 'CURVE':
            points = self.get_points(bpy.context.active_object.data.splines[0])
            minx = min([p.x for p in points])
            maxx = max([p.x for p in points])
            miny = min([p.y for p in points])
        else:
            obj = bpy.context.active_object
            sce = bpy.context.scene

            points = [v.co for v in obj.data.vertices]
            minx = min([p.x for p in points])
            maxx = max([p.x for p in points])
            miny = min([p.z for p in points])

        l = maxx - minx

        del sce.sailflow_model.curvePoints[:]
        for v in sorted(points, key=lambda p: p.x):
            if obj.type == 'CURVE':
                sce.sailflow_model.curvePoints.append(((v.x - minx) / l, (v.y - miny) / l))
            else:
                sce.sailflow_model.curvePoints.append(((v.x - minx) / l, (v.z - miny) / l))

        return {'FINISHED'}


class LibraryLoader(bpy.types.Operator):
    bl_idname = "mesh.load_library"
    bl_label = "Load profiler"
    bl_description = "Generate the sail custom profile"

    filepath = bpy.props.StringProperty(subtype="FILE_PATH")

    def execute(self, context):
        global custom_profile
        print("File Path: ", self.filepath)

        custom_profile = imp.load_source('custom_profile', self.filepath)
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    @classmethod
    def poll(cls, context):
        return context.active_object and context.active_object.type == 'MESH'


class AirProfile(bpy.types.Operator):
    bl_idname = "mesh.airprof"
    bl_label = "Generate Profile"
    bl_description = "Generate the sail profile"

    def findWeight(self, v):

        vco = v
        vpco = []
        mw = bpy.data.objects["Profile"].matrix_world
        for vp in bpy.data.objects["Profile"].data.vertices:
            vpco.append(mw * vp.co)
        vpzl = []
        minim = 100
        for vp in vpco:
            if abs(vco.y - vp.y) < minim:
                minim = abs(vco.y - vp.y)
                vpMin = vp
            vpzl.append(vp.z)
        return (vpMin.z - min(vpzl)) / (max(vpzl) - min(vpzl))

    def getXinEdge(self, vl, e, yc):
        x1 = vl[e[0]].co.x
        x2 = vl[e[1]].co.x
        y1 = vl[e[0]].co.y
        y2 = vl[e[1]].co.y
        if (y1 == y2):
            # Horizontal line
            # print("gXE horizontal line")
            return (max(x1, x2))
        elif abs(x1 - x2) < 0.001:
            # Vertical line
            # print("gXE vert line yc=%2.2f e11[%d] el2[%d]" % (yc, e[0], e[1]))
            return x1
        else:
            # v = mathutils.geometry.intersect_line_line_2d(
            #         Vector((x1, y1, 0)),
            #         Vector((x2, y2, 0)),
            #         Vector((-10, yc, 0)),
            #         Vector((+10, yc, 0)))
            # print("gXE yc=%2.2f e11[%d] el2[%d]" % (yc, e[0], e[1]))
            a = (y1 - y2) / (x1 - x2)
            b = y2 - (a * x2)
            x = (yc - b) / a
            return x
            # if v:
            #     return v.x
            # else:
            #     return 0.0

    def getEdgesCrossing(self, vl, pe, y):
        ce = []
        c1 = (-1, 0)
        c2 = (-1, 0)

        for e in pe:
            #           print("--->Check y=%f e0y=%f e1y=%f"%(y,vl[e[0]].co.y,vl[e[1]].co.y))
            if (vl[e[0]].co.y >= y and vl[e[1]].co.y <= y) or (vl[e[0]].co.y <= y and vl[e[1]].co.y >= y):
                ce.append(e)
        # print("getEdgeCross: finished check found %d edges,scanning..."%(len(ce)))
        if len(ce) > 1:
            minx = 100.0
            maxx = -100.0
            for c in ce:
                if vl[c[0]].co.x < minx:
                    minx = vl[c[0]].co.x
                    c1 = c
                if vl[c[1]].co.x < minx:
                    minx = vl[c[1]].co.x
                    c1 = c
                if vl[c[0]].co.x > maxx:
                    maxx = vl[c[0]].co.x
                    c2 = c
                if vl[c[1]].co.x > maxx:
                    maxx = vl[c[1]].co.x
                    c2 = c
                    # print("            : %d:(%2.2f) %d:(%2.2f)"%(c[0],vl[c[0]].co.x,c[1],vl[c[1]].co.x))
                    # print("            : (%d,%d) minx=%2.4f maxx=%2.4f"%(c[0],c[1],minx,maxx))
        else:
            c1 = (-1, 0)
            c2 = (-1, 0)

        return (c1, c2)

    def extractPerimeterEdges(self, m):

        el = []
        perifEdgeList = []
        for p in m.data.polygons:
            for e in p.edge_keys:
                if e[0] > e[1]:
                    el.append((e[1], e[0]))
                else:
                    el.append(e)
        while el:
            e = el.pop()
            if e in el:
                el.remove(e)
            else:
                perifEdgeList.append(e)

        return perifEdgeList

    def makeTwist(self, a, maxTwist):
        pe = self.extractPerimeterEdges(a)
        miny = 100000
        maxy = -100000
        for v in a.data.vertices:
            if v.co.y < miny:
                miny = v.co.y
            if v.co.y > maxy:
                maxy = v.co.y
        for v in a.data.vertices:
            if v.co.y != 0 and v.co.x != 0:
                e1, e2 = self.getEdgesCrossing(a.data.vertices, pe, v.co.y)
                if e1[0] != -1 and e2[0] != -1:
                    x1 = self.getXinEdge(a.data.vertices, e1, v.co.y)
                    x2 = self.getXinEdge(a.data.vertices, e2, v.co.y)
                    leftx = min(x1, x2)
                    angle = ((v.co.y - miny) / (maxy - miny) * maxTwist) / 180 * 3.14159
                    a.data.vertices[v.index].co.z += sin(angle) * (v.co.x - leftx)
                    a.data.vertices[v.index].co.x = cos(angle) * (v.co.x - leftx) + leftx

    def camber(self, ctx=None, mp=0.1, pp=0.5, weightMode=False, profileMode=False, ctype="NACA", ellipDis=False,
               shrink=False):
        class localVertex:
            def __init__(self, x, y, z, index):
                self.x = x
                self.newx = x
                self.y = y
                self.z = z
                self.index = index

        sce = bpy.context.scene
        a = ctx.active_object

        if a is None:
            a = bpy.data.objects["Sail"]
        for v in a.data.vertices:
            v.co.z = 0
        pe = extractPerimeterEdges(a)
        vxsInPe = []
        for e in pe:
            vxsInPe.append(e[0])
            vxsInPe.append(e[1])
        miny = min([v.co.y for v in a.data.vertices])
        maxy = max([v.co.y for v in a.data.vertices])
        halfSpan = (maxy - miny) / 2
        halfSpanY = halfSpan + miny
        # Calculate the length of the profile by simple splitting
        # into line pieces and summing
        if shrink:
            steps = 1000
            dxf = 1.0 / float(steps)
            yp = 0  # lets assume all profile starts from 0 at x=0
            lineLength = 0.0
            for x in range(1, steps):
                xf = float(x) / float(steps)
                y = profile(xf, mp, pp)
                lineLength += sqrt((y - yp) * (y - yp) + dxf * dxf)
                yp = y
            XshrinkFactor = 1 / lineLength
        else:
            XshrinkFactor = 1.0
        verticesCopy = [localVertex(v.co.x, v.co.y, v.co.z, v.index) for v in a.data.vertices]
        for v in verticesCopy:
            e1, e2 = self.getEdgesCrossing(a.data.vertices, pe, v.y)
            if e1[0] != -1 and e2[0] != -1:
                x1 = self.getXinEdge(a.data.vertices, e1, v.y)
                x2 = self.getXinEdge(a.data.vertices, e2, v.y)
                leftx = min(x1, x2)
                rightx = max(x1, x2)
                x = (v.x - leftx) / (rightx - leftx + 0.0001)
                if v.index in vxsInPe:
                    y = 0.0
                elif ctype == 'NACA':
                    y = profile(x, mp, pp)
                elif ctype == 'CURVE':
                    y = curveProfile(x, sce.sailflow_model.curvePoints)
                elif ctype == 'CUSTOM':
                    heightPerc = (v.y - miny) / (maxy - miny)
                    y = custom_profile.profile(x, heightPerc, v.co.x, v.co.y, miny, maxy)
                if ellipDis:
                    y = y * sqrt(1 - ((v.y - halfSpanY) / halfSpan) ** 2)

                v.z = y * (rightx - leftx) * XshrinkFactor
                v.newx = leftx + (v.x - leftx) * XshrinkFactor
        vs = a.data.vertices
        for v in verticesCopy:
            vs[v.index].co.x = v.newx
            vs[v.index].co.z = v.z

    def execute(self, context):
        print("Called airprofile")
        sce = bpy.context.scene
        self.camber(ctx=context, mp=sce.sailflow_model.m / 100, pp=sce.sailflow_model.p / 100,
                    ctype=sce.sailflow_model.t,
                    ellipDis=sce.sailflow_model.ellipDis,
                    shrink=sce.sailflow_model.shrink)
        if sce.sailflow_model.twist:
            self.makeTwist(context.active_object, sce.sailflow_model.tw)
        return {'FINISHED'}

    @classmethod
    def poll(cls, context):
        return context.active_object and context.active_object.type == 'MESH' and context.mode == 'OBJECT'


class printPDF(bpy.types.Operator):
    bl_idname = "mesh.print_pdf"
    bl_label = "Print PDF"
    bl_description = "print PDF of selected panels"

    unitToMm = 1000

    filepath = bpy.props.StringProperty(subtype="FILE_PATH")

    def containedIn(self, v, w, h):
        return (v.x >= 0 and v.x <= w and v.y >= 0 and v.y <= h)

    def clipping(self, v1, v2, w, h):
        vout1 = None
        vout2 = None
        if self.containedIn(v1, w, h):
            # print(" clipping: V1 already in ")
            vout1 = v1
        if self.containedIn(v2, w, h):
            # print(" clipping: V2 already in ")
            vout1 = v2

        v = mathutils.geometry.intersect_line_line_2d(v1, v2, Vector((0, 0, 0)), Vector((w, 0, 0)))
        if v:
            # print(" clipping: h bassa")
            if vout1 == None:
                vout1 = v
            elif (v != vout1):
                vout2 = v
                return vout1, vout2
        v = mathutils.geometry.intersect_line_line_2d(v1, v2, Vector((w, 0, 0)), Vector((w, h, 0)))
        if v:
            # print(" clipping: v destra")
            if vout1 == None:
                vout1 = v
            elif (v != vout1):
                vout2 = v
                return vout1, vout2
        v = mathutils.geometry.intersect_line_line_2d(v1, v2, Vector((w, h, 0)), Vector((0, h, 0)))
        if v:
            # print(" clipping: h alta")
            if vout1 == None:
                vout1 = v
            elif (v != vout1):
                vout2 = v
                return vout1, vout2
        v = mathutils.geometry.intersect_line_line_2d(v1, v2, Vector((0, h, 0)), Vector((0, 0, 0)))
        if v:
            # print(" clipping: v sinistra")
            if vout1 == None:
                vout1 = v
            elif (v != vout1):
                vout2 = v
        return vout1, vout2

    def makePDF(self, a, pe, vxs, fmt, freeTxt='', mp=False, offSetSet=0, overlapSet=0, w=0, h=0):
        xs = []
        ys = []

        for e in pe:
            v0 = a.matrix_world * vxs[e[0]].co
            v1 = a.matrix_world * vxs[e[1]].co
            xs.append(v0.x)
            xs.append(v1.x)
            ys.append(v0.y)
            ys.append(v1.y)

        min_x = min(xs)
        min_y = min(ys)
        vxs_pix = []
        for v in vxs:
            vc = a.matrix_world * v.co
            vxs_pix.append(Vector([(vc.x - min_x) * self.unitToMm, (vc.y - min_y) * self.unitToMm]))

            # print("Vector")
            # for i in range(len(vxs_pix)):
            # print(i,vxs_pix[i])

        fname = self.filepath
        pdf = FPDF(format=fmt, pw=w, ph=h)
        offSet = pdf.dimension()[0] * offSetSet / 100
        # percentage of the offSet used for overlap between pages
        overlap = pdf.dimension()[0] * overlapSet / 100
        pdf.set_compression(False)
        if not mp:
            pdf.add_page()
            for e in pe:
                p1 = vxs_pix[e[0]]
                p2 = vxs_pix[e[1]]
                pdf.line(p1.x + offSet, (pdf.dimension()[1] - p1.y) - offSet, p2.x + offSet,
                         (pdf.dimension()[1] - p2.y) - offSet)

            # Make the header
            pdf.set_x(0)
            pdf.set_y(0)
            date = str(time.localtime().tm_mday) + "/" + str(time.localtime().tm_mon) + "/" + str(
                time.localtime().tm_year)
            pdf.set_font('Arial', 'B', 12)
            # pdf.cell(w=75, h=10, txt=a.name, border=1, ln=2, align='C')
            # pdf.set_font(family='Arial', size=12)
            pdf.cell(w=75, h=5, txt=date, border=1, ln=2, align='C')
            # pdf.cell(w=75, h=5, txt=freeTxt + " ", border=1, ln=0, align='C')
        else:
            clipSizeH = int(pdf.dimension()[0] - 2 * offSet)
            clipSizeV = int(pdf.dimension()[1] - 2 * offSet)
            minx = min([v.x for v in vxs_pix])
            maxx = max([v.x for v in vxs_pix])
            miny = min([v.y for v in vxs_pix])
            maxy = max([v.y for v in vxs_pix])
            numPagesH = int((maxx - minx) / (clipSizeH - overlap)) + 1
            numPagesV = int((maxy - miny) / (clipSizeV - overlap)) + 1

            for h in range(0, numPagesH):
                for v in range(0, numPagesV):
                    pdf.add_page()
                    line = 0
                    for e in pe:
                        # Translate the points
                        cross1 = Vector([vxs_pix[e[0]].x, vxs_pix[e[0]].y])
                        cross2 = Vector([vxs_pix[e[1]].x, vxs_pix[e[1]].y])
                        cross1.x = cross1.x - float(clipSizeH * h) + float(overlap * h)
                        cross2.x = cross2.x - float(clipSizeH * h) + float(overlap * h)
                        cross1.y = cross1.y - float(clipSizeV * v) + float(overlap * v)
                        cross2.y = cross2.y - float(clipSizeV * v) + float(overlap * v)

                        if not (self.containedIn(cross1, clipSizeH, clipSizeV) and
                                    self.containedIn(cross2, clipSizeH, clipSizeV)):
                            cross1, cross2 = self.clipping(cross1, cross2, clipSizeH, clipSizeV)

                        if cross1 != None and cross2 != None:
                            pdf.line(cross1.x + offSet,
                                     (pdf.dimension()[1] - cross1.y) - offSet,
                                     cross2.x + offSet,
                                     (pdf.dimension()[1] - cross2.y) - offSet)
                            line += 1

                    if (h + v != 0):
                        pdf.set_line_width(0.2)
                        pdf.line(0.0, pdf.dimension()[1] - offSet - overlap, pdf.dimension()[0],
                                 pdf.dimension()[1] - offSet - overlap)
                        pdf.line(offSet + overlap, 0.0, offSet + overlap, pdf.dimension()[1])

                    pdf.set_x(0)
                    pdf.set_y(0)
                    # date = str(time.localtime().tm_mday) + "/" + str(time.localtime().tm_mon) + "/" + str(time.localtime().tm_year)
                    # pdf.set_font('Arial', 'B', 24)
                    pdf.set_font(family='Arial', size=12)
                    pdf.cell(w=75, h=10, txt=a.name + "H:" + str(h) + " V:" + str(v), border=1, ln=2, align='C')
                    # pdf.cell(w=75, h=5, txt=date, border=1, ln=2, align='C')
                    # pdf.cell(w=75, h=5, txt=freeTxt + " ", border=1, ln=0, align='C')

        pdf.output(fname, 'F')

    def execute(self, context):
        sce = bpy.context.scene
        a = context.active_object
        el = []
        for p in a.data.polygons:
            for e in p.edge_keys:
                el.append((e[1], e[0]))
        vxs = a.data.vertices
        se = sce.sailflow_model
        self.makePDF(a, el, vxs, se.paperFormat, se.freeText, se.multiPages, se.margin, se.overlap, se.paperWidth,
                     se.paperHeight);
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    @classmethod
    def poll(cls, context):
        return context.active_object and context.active_object.type == 'MESH'


class outputAscii(bpy.types.Operator):
    bl_idname = "mesh.print_ascii"
    bl_label = "Print ASCII"
    bl_description = "print ASCII description of selected panel"
    filepath = bpy.props.StringProperty(subtype="FILE_PATH")

    def execute(self, context):
        sce = bpy.context.scene

        plotDX = sce.sailflow_model.asciiDx

        obj = context.active_object
        # finds the minimum and max X of the panel scanning all vertices
        vxs = [(v.co.x, v.co.y, idx) for idx, v in enumerate(obj.data.vertices)]
        xmin = min(vxs, key=lambda t: t[0])[0]
        xmax = max(vxs, key=lambda t: t[0])[0]
        ymin = min(vxs, key=lambda t: t[1])[1]
        ymax = max(vxs, key=lambda t: t[1])[1]

        # get all periferal edges
        pes = extractPerimeterEdges(obj)
        vxs = obj.data.vertices
        pnts = []
        # find corners
        corners = []
        for e in pes:
            for e2 in pes:
                if (e[0] == e2[0]) or (e[0] == e2[1]):
                    v1 = vxs[e[0]].co - vxs[e[1]].co
                    v2 = vxs[e2[0]].co - vxs[e2[1]].co
                    if radians(10) < v1.angle(v2) < radians(120):
                        if not vxs[e[0]].co in corners:
                            corners.append(vxs[e[0]].co)
                elif (e[1] == e2[0]) or (e[1] == e2[1]):
                    v1 = vxs[e[0]].co - vxs[e[1]].co
                    v2 = vxs[e2[0]].co - vxs[e2[1]].co
                    if radians(10) < v1.angle(v2) < radians(120):
                        if not vxs[e[1]].co in corners:
                            corners.append(vxs[e[1]].co)

        xscan = xmin
        while xscan <= xmax:
            for e in pes:
                iv = mathutils.geometry.intersect_line_line_2d(Vector((xscan, ymin - 1, 0)), \
                                                               Vector((xscan, ymax + 1, 0)), \
                                                               vxs[e[0]].co, vxs[e[1]].co)
                if iv:
                    pnts.append((xscan, iv[0], iv[1]))
            xscan = xscan + plotDX / 100

        f = open(self.filepath, 'w')
        f.write("---------------------------------\n")
        f.write(("Panel %s\n") % (obj.name))
        f.write("---------------------------------\n")
        f.write("Coordinates of the corners \n")
        for c in corners:
            f.write(("(%2.3f,%2.3f) ") % (c.x - xmin, c.y - ymin))

        f.write(("\nCoordinates %d cm at distance\n") % (plotDX))
        for p in pnts:
            s = ("[%2.3f] %2.3f %2.3f\n") % (p[0], p[1] - xmin, p[2] - ymin)
            f.write(s)
        f.close()
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    @classmethod
    def poll(cls, context):
        return context.active_object and context.active_object.type == 'MESH'


class AirFoilSettings(bpy.types.PropertyGroup):
    curveTypes = [
        ("NACA", "Naca 4", "", 1),
        ("CUSTOM", "Profile on routine", "", 2),
        ("DAT", "Load DAT and/or lofting", "", 3),
        ("CURVE", "Profile Bezier or Mesh", "", 4)
    ]
    m = IntProperty(
        name="% Max Camber",
        description="Maximum Camber (% of chord)",
        default=0,
        min=0,
        max=50)
    p = IntProperty(
        name="% Camber Pos",
        description="Position of Camber (% of chord)",
        default=0,
        min=20,
        max=70)

    t = EnumProperty(name="Curve", default="NACA", items=curveTypes)

    shrink = BoolProperty(name="Apply Shrink", default=False)
    curve = BoolProperty(name="Apply Curve", default=False)
    twist = BoolProperty(name="Apply Twist", default=False)
    ellipAmount = bpy.props.FloatProperty(name="Elliptical Amount", default=0.0, min=0.0, max = 1.0)
    ellipCenter = bpy.props.FloatProperty(name="Vertical Display", default=0.0, min=-1.0, max = 1.0)
    tw = IntProperty(name="Twist Angle", description="Value of angle", default=0, min=0, max=90)
    energyMinimizer = bpy.props.BoolProperty(name="Stress Relief", default=False)
    maxDeformation = bpy.props.FloatProperty(name="Min stress reduction", default=0.001, min=0.000001, max=0.1)
    deltaDeformation = bpy.props.IntProperty(name="accuracy", default=3, min=1, max=10)
    polSeed = bpy.props.IntProperty(name="Use start face", default=-1)
    useSeed = bpy.props.BoolProperty(name="Start from face", default=False)
    resEnergy = bpy.props.FloatProperty(name="Residual Strain", default=0)

    # paperFormat = bpy.props.StringProperty (name="Paper Size",description="others,4a0,2a0,a0,a1,a2,a3,a4",default='a0')
    freeText = bpy.props.StringProperty(name="Free Text", description="Anything appearing in the PDF",
                                        default='free text')
    margin = IntProperty(name="% of page margin",
                         default=10,
                         min=0,
                         max=25)
    overlap = IntProperty(name="% of overlap bw pages",
                          default=10,
                          min=0,
                          max=25)
    paperWidth = IntProperty(name="Width", max=1000)
    paperHeight = IntProperty(name="Height", max=3500)

    ellipDis = BoolProperty(name="Eliptic Distribution", default=False)

    paperSizes = [
        ("4A0", "4a0", "", 1),
        ("2A0", "2a0", "", 2),
        ("A0", "a0", "", 3),
        ("A1", "a1", "", 4),
        ("A2", "a2", "", 5),
        ("A3", "a3", "", 6),
        ("A4", "a4", "", 7),
        ("Other", "other", "", 8)]

    paperFormat = bpy.props.EnumProperty(name="Paper Size", default="A4", items=paperSizes)
    freeText = bpy.props.StringProperty(name="Free Text", description="Anything appearing in the PDF",
                                        default='free text')
    multiPages = bpy.props.BoolProperty(name="Multi pages", default=False)

    resolution = IntProperty(name="Resolution", min=10, max=200)
    asciiDx = IntProperty(name="X step cm", min=1, max=100)
    steps = IntProperty(name="Loft steps", min=5, max=100)
    spans = IntProperty(name="Loft spans", min=5, max=100)

    curvePoints = []
    curveName = bpy.props.StringProperty(name="object name")


def select():
    # print(bpy.context.mode)
    if bpy.context.mode == "OBJECT":
        obj = bpy.context.object
        sel = len(bpy.context.selected_objects)

        if sel == 0:
            bpy.selection = []
        else:
            if sel == 1:
                bpy.selection = []
                bpy.selection.append(obj)
            elif sel > len(bpy.selection):
                for sobj in bpy.context.selected_objects:
                    if (sobj in bpy.selection) == False:
                        bpy.selection.append(sobj)

            elif sel < len(bpy.selection):
                for it in bpy.selection:
                    if (it in bpy.context.selected_objects) == False:
                        bpy.selection.remove(it)

                        # on edit mode doesnt work well


# executes selection by order at 3d view
class Selection(bpy.types.Header):
    bl_label = "Selection"
    bl_space_type = "VIEW_3D"

    def __init__(self):
        # print("hey")
        select()

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Sel: " + str(len(bpy.selection)))


def register():
    bpy.utils.register_module(__name__)
    bpy.types.Scene.sailflow_model = bpy.props.PointerProperty(type=AirFoilSettings,
                                                               name="Airfoil Model",
                                                               description="Setting of the AirFoil")


#    bpy.utils.register_class(Selection)

def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.utils.unregister_class(Selection)


if __name__ == "__main__":
    register()
