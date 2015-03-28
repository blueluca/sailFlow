bl_info = {
    "name": "Create Sailprofile",
    "description": "Creates a profile for a sail",
    "author": "blueluca",
    "version": (0, 2, 3),
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

from bpy.props import IntProperty, EnumProperty, BoolProperty, StringProperty
from bpy.types import Operator
from fpdf import FPDF
from mathutils import Vector, Euler, geometry
import mathutils
#import pydevd


#PYDEV_SOURCE_DIR = 'C:/eclipse/plugins/org.python.pydev_3.4.1.201403181715/pysrc'
  
#if sys.path.count(PYDEV_SOURCE_DIR) < 1:
#   sys.path.append(PYDEV_SOURCE_DIR)
  

custom_profile = None   

Vxs = []
FlatVxs = []
F = []
VxFlat = []

def profile(x, m, p):
#=====================================================================================
# The profile calculator
#
# inputs 
#   x  : point in the range 0 to 1 to calculate
#   m : camber percentage
#   p : camber position percentage
#=====================================================================================
    if (x < 0):
        return 0
    if (x < p):
        y = (m / (p * p)) * (2 * p * x - x * x);
    else:
        y = (m / ((1 - p) * (1 - p))) * ((1 - 2 * p) + 2 * p * x - x * x);
    # debug(("Profile return Z="+str(y))
    return y

def curveProfile(x,vxs):

    for i in vxs:
        if i.x > x:
            found = True
            break       
    if found: 
        return i.y
    else:
        return 0
 
 
class MyPoly:
    l1 = 0  # p1 - p3
    l2 = 0  # p1 - p2
    l3 = 0  # p2 - p3

    def clockWise(self,co1, co2, co3):
        A = (co2.x - co1.x) * (co2.y + co1.y)
        A = A + (co3.x - co2.x) * (co3.y + co2.y)
        A = A + (co1.x - co3.x) * (co1.y + co3.y)
        return A > 0

    def doubleCircle(self,d, R, r):
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
 
    def __str__ (self):
        return "MyPoly [" + str(self.idx) + "] " + str(self.p1) + "," + str(self.p2) + "," + str(self.p3)  # +" L:"+str(self.l1)+","+str(self.l2)+","+str(self.l3)
 
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
 
    def NEWsetToFlat(self,a,b,flatted):
        xa = self.x(a)
        ya = self.y(a)
        xb = self.x(b)
        yb = self.y(b)
        print(self,end='')
        print("SetToFlat (a=%d,b=%d,f=%d) xa=%2.2f,ya=%2.2f,xb=%2.2f,yb=%2.2f"%(a,b,flatted,xa,ya,xb,yb))
        d = self.edgeLength(a, b)
        ra = self.edgeLength(a,flatted)
        rb = self.edgeLength(b,flatted)
        print("          d=%2.2f ra=%2.2f rb=%2.2f"%(d,ra,rb))
        k = 0.25*math.sqrt((((ra+rb)*(ra+rb)-(d*d))*(d*d-(ra-rb)*(ra-rb))))
        x1 = 0.5*(xb+xa)+0.5*(xb-xa)*(ra*ra - rb*rb)/(d*d) + 2*(yb-ya)*k/(d*d)
        x2 = 0.5*(xb+xa)+0.5*(xb-xa)*(ra*ra - rb*rb)/(d*d) - 2*(yb-ya)*k/(d*d)
        y1 = 0.5*(yb+ya)+0.5*(yb-ya)*(ra*ra - rb*rb)/(d*d) - 2*(xb-xa)*k/(d*d)
        y2 = 0.5*(yb+ya)+0.5*(yb-ya)*(ra*ra - rb*rb)/(d*d) + 2*(xb-xa)*k/(d*d)
        print("          x1=%2.2f y1=%2.2f x2=%2.2f y2=%2.2f"%(x1,y1,x2,y2))

        va = Vector([xa,ya])
        vb = Vector([xb,yb])
        vc1 = Vector([x1,y1])
        if self.clockWise(va,vb,vc1):
            self.setx(flatted,x1)
            self.sety(flatted,y1)
        else:
            self.setx(flatted,x2)
            self.sety(flatted,y2)  
        
        #=======================================================================
        # print("  p1-p3 diff %2.2f perc %2.2f"%(self.l1- (Vxs[self.p1].co - Vxs[self.p3].co).length,
        #                                        (self.l1- (Vxs[self.p1].co - Vxs[self.p3].co).length)/self.l1),end='')
        # print("  p1-p2 diff %2.2f perc %2.2f"%(self.l2 - (Vxs[self.p1].co - Vxs[self.p2].co).length,
        #                                        (self.l2 - (Vxs[self.p1].co - Vxs[self.p2].co).length)/self.l2),end='')
        # print("  p2-p3 diff %2.2f perc %2.2f"%(self.l3 - (Vxs[self.p2].co - Vxs[self.p3].co).length,
        #                                        (self.l3 - (Vxs[self.p2].co - Vxs[self.p3].co).length)/self.l3))
        #       
        #=======================================================================
        
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
 
    def flatten (self):
        global Vxs
        global VxFlat

        if (self.p1 in VxFlat) and (self.p2 in VxFlat) and (self.p3 in VxFlat) :
#            print(self.idx," already flat ") 
            return 0
        elif (((self.p1 not in VxFlat) and (self.p2 not in VxFlat)) or
             ((self.p1 not in VxFlat) and (self.p3 not in VxFlat)) or
             ((self.p2 not in VxFlat) and (self.p3 not in VxFlat))):
            self.flatAllVertices()
        else:
            if  self.p3 not in VxFlat:
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

    
    def __init__(self,i):
        self.idx = i
        self.adj = []
        self.energy = 0.0
        self.evx = 0
#        print("Create enVX ",i)
        
    def addAdjacent(self,iandl):
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
            self.energy = self.energy + dl*dl/ v[1]
#            print("...calcEnergy [",self.idx,"-",v[0],"] dl=",dl," en=",self.energy)
#        print("...calcEnery total ", self.idx,self.energy)
        return self.energy
    
    def isAdjacent(self,i):
        for a in self.adj:
            if a[0] == i:
                return True
        return False
    
    def calcDeltaEnergy(self,delta):
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
#        print("calcDeltaEnergy e,e+x,e-x,e+y,e-y ",self.idx,self.energy,self.overallpdx,self.overallmdx,self.overallpdy,self.overallmdy)
#        print("                e,+x,-x,+y,-y     ",self.idx,self.energy,self.energy-self.overallpdx,self.energy-self.overallmdx,self.energy-self.overallpdy,self.energy-self.overallmdy)

    
# =================================================================
#             
#             
#             FLATTENER
# 
#            
# ==================================================================            
    
class Flattener(bpy.types.Operator):
    bl_idname = "mesh.flattener"
    bl_label = "Flat surface"
    bl_description = "bla bla bla"
    
            
    def findAdjacentNonFlat(self, Polys, p):
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
             

    def updateNodesMovToGain(self,e,delta):
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
            self.nodesMovToGain[e.evx] = [0,0,0.0]        

    def minimizeEnergy(self, F, maxDeformation, deltaDeformation):
        global Vxs

        vlist = []
        EVs = []  
        self.nodesMovToGain = []
              
        print("="*80)
        print(" "*40,"start"," "*40)
        print("="*80)
        for p in F:
            vlist = vlist + p.vertices
        #remove duplicates 
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
        #print("end database build")
    
        MAXCOUNT = 5000
        delta = 0.1
        while delta >= 10**-deltaDeformation:
            print("---------------> NEW Loop minimize delta=",delta,"minimum decrease",maxDeformation)
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
                    self.updateNodesMovToGain(e,delta)
                    # if the gain is higher than the threshold or the best record the node
    
                    #if self.nodesMovToGain[evIdx][2]>0:
                        #print("Node ",e.idx,"energy=",e.energy," dec energy=",self.nodesMovToGain[evIdx][2])
                        
                    if self.nodesMovToGain[evIdx][2] > maxGain:
                        #print("..new max decr energy for",e.idx,"dec energy=",self.nodesMovToGain[evIdx][2])
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
                    #print("Selected node ",maxE.idx,"energy=",maxE.energy," dec energy=",self.nodesMovToGain[maxE.evx][2])
                    # apply the change in position
                    Vxs[maxE.idx].co.x += self.nodesMovToGain[maxIdx][0]
                    Vxs[maxE.idx].co.y += self.nodesMovToGain[maxIdx][1]
                    # update the energy content and the deltas for the winnig vertex
                    maxE.calcEnergy()
                    maxE.calcDeltaEnergy(delta)               
                    # update the database of nodeMovToGain for the winning vertex       
                    self.updateNodesMovToGain(maxE,delta)  
                    #print("..updating database for adjacent")                          
                    for ev in EVs:
                        if maxE.isAdjacent(ev.idx):
                            #print("... update",ev.idx)
                            ev.calcEnergy()
                            if ev.energy > maxDeformation:                        
                                ev.calcDeltaEnergy(delta)
                                self.updateNodesMovToGain(ev,delta)
                            else:
                                self.nodesMovToGain[ev.evx] = [0.0,0.0,0.0]
                                
                    maxIdx = None      
                    maxGain = maxDeformation
                    #print("scanning database again")              
                    for evIdx in range(len(EVs)):
                        nmtg = self.nodesMovToGain[evIdx]
                        e = EVs[evIdx] 
                        #if nmtg[2] > 0:  
                            #print("..node ",EVs[evIdx].idx,nmtg)
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

        maxRE=0
        idx=0
        for e in EVs:
            if e.energy > maxRE:
                maxRE = e.energy
                idx = e.idx
        print("exit minimize energy, max residual energy node=",idx," energy=",maxRE)
        return maxRE
                     
                
    def makeItFlat(self, obj, energyMinimizer, maxDeformation, deltaDeformation,polSeed):
        # Variables
        #
        # V available polygons
        # A active polygons
        # F flatten polygons
        # vt face of the available
        # ft is vt flatten version
        # at temporary active polygon
     
        global Vxs
        global F
        global VxFlat
     
        A = []
        V = []
        Vxs = []
        VxFlat = []
        me = obj.data
        Vxs = me.vertices[:]
        count = 0
#        pydevd.settrace(stdoutToServer=True, stderrToServer=True, suspend=True)
        s = None
        for p in me.polygons:
            if p.select:
                vt = MyPoly(p, idx=p.index)
                V.append(vt)
        
        if polSeed != -1:
            for vt in V:
                if vt.idx == polSeed:
                    A.append(vt)
                    V.remove(vt)
                    break;
         
        while V or A:
            if A:
                s = A.pop(0)
            else:
                s = V.pop()
#                print("==============================restart from ",s.idx)
            # Collect all the adjacent triangle
            # and put it in Active collection "A"
            found = True
            while(found):
                at = self.findAdjacentNonFlat(V, s)
                if at:
#                    print("..adj is",at.idx)
                    A.append(at)
                    V.remove(at) 
                    found = True
                else:
                    found = False
#            print(count," Flattening: " + str(s))
            count += 1
            s.flatten()
            F.append(s)
    
        if energyMinimizer:
            return(self.minimizeEnergy(F, maxDeformation, deltaDeformation))
        else:
            return 0
       
    @classmethod
    def poll(cls, context):
        return context.active_object and context.active_object.type == 'MESH'
    
    def execute(self, context): 
        global Vxs
        global F
        

        sce = bpy.context.scene  
        # Record the selected faces
        bpy.context.active_object.update_from_editmode()
            
        if F:
            del F[:]
        
        if sce.sailflow_model.useSeed == False:
            self.makeItFlat(bpy.context.active_object, sce.sailflow_model.energyMinimizer, sce.sailflow_model.maxDeformation, 
                            sce.sailflow_model.deltaDeformation,-1)
        else:        
            sce.sailflow_model.resEnergy = self.makeItFlat(bpy.context.active_object, sce.sailflow_model.energyMinimizer, sce.sailflow_model.maxDeformation, 
                                            sce.sailflow_model.deltaDeformation,sce.sailflow_model.polSeed)
        # Complet list of vertices used
        vIdxList = []
        faces = []
        vertices = []
        for p in F:
            vIdxList.append(p.p1)
            vIdxList.append(p.p2)
            vIdxList.append(p.p3)
        vIdxList.sort()
        # The following remove duplicates
        vIdxList = list(set(vIdxList))
        for p in F:
            p.p1 = vIdxList.index(p.p1)
            p.p2 = vIdxList.index(p.p2)
            p.p3 = vIdxList.index(p.p3)
            faces.append([p.p1, p.p2, p.p3])
        for vidx in vIdxList:
            vertices.append([Vxs[vidx].co.x, Vxs[vidx].co.y, Vxs[vidx].co.z])

#===============================================================================
# # ============================================================        
#         faces = []
#         vertices = []
#         count = 0
#         idx=0
#         while F:
#             for p in F:
#                 faces.append([p.p1, p.p2, p.p3])
#                 idx = idx + 1
#                 F.remove(p)
#         for v in Vxs:
#             vertices.append([v.co.x, v.co.y, v.co.z])
# 
# # =========================================================
#===============================================================================

        fmesh = bpy.data.meshes.new("panel")
        fmesh.from_pydata(vertices, [], faces)
        # Update the displayed mesh with the new data
        fmesh.update()
        fobj = bpy.data.objects.new("panel", fmesh)
#        fobj.matrix_world = bpy.context.active_object.matrix_world
        fobj.scale = bpy.context.active_object.scale
        scene = bpy.context.scene  
        scene.objects.link(fobj)  
        return {'FINISHED'} 

class PMemorize(bpy.types.Operator):
    bl_idname = "mesh.panelmem"
    bl_label = "Panel Memorize"
    bl_description = "bla bla bla"
    
    def execute(self,context):
        bpy.ops.object.mode_set(mode = 'OBJECT')
        bpy.ops.object.mode_set(mode = 'EDIT')
                
        sce = bpy.context.scene        
        faces = bpy.context.active_object.data.polygons
        for face in faces:
            if face.select:
                if sce.sailflow_model.panelNumber == 1:
                    sce.sailflow_model.panel1.append(face.index)
                elif sce.sailflow_model.panelNumber == 2:
                    sce.sailflow_model.panel2.append(face.index)
                elif sce.sailflow_model.panelNumber == 3:
                    sce.sailflow_model.panel3.append(face.index)
                elif sce.sailflow_model.panelNumber == 4:
                    sce.sailflow_model.panel4.append(face.index)
                elif sce.sailflow_model.panelNumber == 5:
                    sce.sailflow_model.panel5.append(face.index)
                elif sce.sailflow_model.panelNumber == 6:
                    sce.sailflow_model.panel6.append(face.index)
                elif sce.sailflow_model.panelNumber == 7:
                    sce.sailflow_model.panel7.append(face.index)
                elif sce.sailflow_model.panelNumber == 8:
                    sce.sailflow_model.panel8.append(face.index)
                elif sce.sailflow_model.panelNumber == 9:
                    sce.sailflow_model.panel9.append(face.index)
                elif sce.sailflow_model.panelNumber == 10:
                    sce.sailflow_model.panel10.append(face.index)
                    
        return {'FINISHED'} 

class PRecall(bpy.types.Operator):
    bl_idname = "mesh.panelrec"
    bl_label = "Panel Recall"
    bl_description = "bla bla bla"

    def execute(self,context):
        sce = bpy.context.scene        
        faces = bpy.context.active_object.data.polygons
        bpy.ops.mesh.select_all(action='DESELECT')        
        bpy.ops.object.mode_set(mode = 'OBJECT')
        if sce.sailflow_model.panelNumber == 1:
            p = sce.sailflow_model.panel1
        elif sce.sailflow_model.panelNumber == 2:
            p = sce.sailflow_model.panel2
        elif sce.sailflow_model.panelNumber == 3:
            p = sce.sailflow_model.panel3
        elif sce.sailflow_model.panelNumber == 4:
            p = sce.sailflow_model.panel4
        elif sce.sailflow_model.panelNumber == 5:
            p = sce.sailflow_model.panel5
        elif sce.sailflow_model.panelNumber == 6:
            p = sce.sailflow_model.panel6
        elif sce.sailflow_model.panelNumber == 7:
            p = sce.sailflow_model.panel7
        elif sce.sailflow_model.panelNumber == 8:
            p = sce.sailflow_model.panel8
        elif sce.sailflow_model.panelNumber == 9:
            p = sce.sailflow_model.panel9
        elif sce.sailflow_model.panelNumber == 10:
            p = sce.sailflow_model.panel10  
        
        for fidx in p:
            faces[fidx].select = True
            print("set face ",fidx)
        bpy.ops.object.mode_set(mode = 'EDIT')

        return {'FINISHED'} 
                    
    
class VIEW3D_PT_airprofile_print(bpy.types.Panel):
    bl_label = "PDF Generation"
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
            col.prop(sce.sailflow_model,"paperWidth")
            col.prop(sce.sailflow_model,"paperHeight")
 
        
        col = layout.column(align=True)
        col.operator("mesh.print_pdf")
        
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
             
        elif sce.sailflow_model.t == "THREE":
            col = layout.column(align=True)
            col.label(text="High Section:")
            col.prop(sce.sailflow_model, "sec3M")
            col.prop(sce.sailflow_model, "sec3P")
            
            col = layout.column(align=True)
            col.label(text="Middle Section:")
            col.prop(sce.sailflow_model, "sec2M")
            col.prop(sce.sailflow_model, "sec2P")
            col.prop(sce.sailflow_model, "sec2H")

            col = layout.column(align=True)            
            col.label(text="Low Section:")
            col.prop(sce.sailflow_model, "sec1M")
            col.prop(sce.sailflow_model, "sec1P")
            
        elif sce.sailflow_model.t == "CUSTOM":
            col.operator("mesh.load_library")
            col = layout.column(align=True)  
            
        elif sce.sailflow_model.t == "CURVE":
            col.operator("mesh.bezier_aquire")
            col.operator("mesh.curve_aquire")
            col = layout.column(align=True)  

        col = layout.column(align=True)
        col.prop(sce.sailflow_model, "weight")
        
        row = col.row(align=True)
        row.prop(sce.sailflow_model, "twist")
        row.prop(sce.sailflow_model, "tw")
        row = col.row(align=True)
        row.prop(sce.sailflow_model,"ellipDis")
# 
#       col = layout.column(align=True)    
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
    
    def draw(self, context):
        sce = bpy.context.scene
        layout = self.layout
        col = layout.column(align=True)    

        col.prop(sce.sailflow_model,"energyMinimizer")
        col.prop(sce.sailflow_model,"deltaDeformation")
        col.prop(sce.sailflow_model,"maxDeformation") 
        col.prop(sce.sailflow_model,"useSeed")
        col.prop(sce.sailflow_model,"polSeed")
        col.prop(sce.sailflow_model,"resEnergy")
        
        col = layout.column(align=True)    
        col.operator("mesh.flattener")
        
class VIEW3D_PT_setPanel(bpy.types.Panel):
    bl_label = "Panels"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Sailflow Design"
        
    def draw(self, context):
        sce = bpy.context.scene
        layout = self.layout

        col = layout.column(align=True)    
        col.prop(sce.sailflow_model,"panelNumber")
        
        col = layout.column(align=True)    
        col.operator("mesh.panelmem")

        col = layout.column(align=True)    
        col.operator("mesh.panelrec")

class CurveAquire(bpy.types.Operator):
    bl_idname = "mesh.curve_aquire"
    bl_label = "Aquire Mesh Curve"
    bl_description = "Generate the sail custom profile from curve"    

    def execute (self, context):
        print("called aquire")
        obj = bpy.context.active_object
        sce = bpy.context.scene

        points = []
        for v in obj.data.vertices:
            points.append(v.co)
            
        maxx = -1000
        minx = 1000
        miny = 1000
        
        for p in points:
            if p.x < minx:
                minx = p.x
            elif p.x > maxx:
                maxx = p.x
            if p.y < miny:
                miny = p.y

        l = maxx-minx
                         
        del sce.sailflow_model.curvePoints[:]
        p = list(points)
        points = p
        for v in points:
            v.x = (v.x - minx)/l
            v.y = (v.y - miny)/l
            sce.sailflow_model.curvePoints.append(v)
        return {'FINISHED'}

from mathutils.geometry import interpolate_bezier

class BezierAquire(bpy.types.Operator):
    bl_idname = "mesh.bezier_aquire"
    bl_label = "Aquire Bezier Curve"
    bl_description = "Generate the sail custom profile from curve"    

    def get_points(self,sp, clean=True):
        
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
            good = [v for i, v in enumerate(old[:-1]) if not old[i] == old[i+1]]
            good.append(old[-1])
            return good
                
        return master_point_list

    def execute (self, context):
        print("called aquire")
        obj = bpy.context.active_object
        sce = bpy.context.scene

        points = self.get_points(bpy.context.active_object.data.splines[0])

        maxx = -1000
        minx = 1000
        miny = 1000
        
        for p in points:
            if p.x < minx:
                minx = p.x
            elif p.x > maxx:
                maxx = p.x
            if p.y < miny:
                miny = p.y

        l = maxx-minx
                         
        del sce.sailflow_model.curvePoints[:]
        for v in points:
            v.x = (v.x - minx)/l
            v.y = (v.y - miny)/l
            sce.sailflow_model.curvePoints.append(v)
        return {'FINISHED'}
    

      
class LibraryLoader(bpy.types.Operator):
    bl_idname = "mesh.load_library"
    bl_label = "Load profiler"
    bl_description = "Generate the sail custom profile"

    
    filepath = bpy.props.StringProperty(subtype="FILE_PATH")
    
    def execute(self, context):
        global custom_profile
        print("File Path: ",self.filepath)

        custom_profile = imp.load_source('custom_profile',self.filepath)        
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

    def findWeight(self,v):
        
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
        return(vpMin.z - min(vpzl)) / (max(vpzl) - min(vpzl))         
    
    def getXinEdge(self,vl, e, yc):
        x1 = vl[e[0]].co.x
        x2 = vl[e[1]].co.x
        y1 = vl[e[0]].co.y
        y2 = vl[e[1]].co.y
        if (y1 == y2):
            return(max(x1,x2))
        else:
            v = mathutils.geometry.intersect_line_line_2d(Vector((x1, y1, 0)),
                                                      Vector((x2, y2, 0)),
                                                      Vector((-10, yc, 0)),
                                                      Vector((+10, yc, 0)))
            if v: 
                return v.x
            else:
                return 0.0    
        
    def getEdgesCrossing(self,vl, pe, y):
        ce = []
        c1 = (-1, 0)
        c2 = (-1, 0)
        
        for e in pe:
 #           print("--->Check y=%f e0y=%f e1y=%f"%(y,vl[e[0]].co.y,vl[e[1]].co.y))
            if (vl[e[0]].co.y >= y and vl[e[1]].co.y <= y) or (vl[e[0]].co.y <= y and vl[e[1]].co.y >= y):
                ce.append(e)
 #       print("----> finished check ce len=",len(ce))
        if len(ce) > 1:
                minx = 100.0
                maxx = -100.0
                for c in ce:
                    if vl[c[0]].co.x < minx:
                        minx = vl[c[0]].co.x
                        c1 = c
                    if vl[c[0]].co.x > maxx:
                        maxx = vl[c[0]].co.x
                        c2 = c
        else:
            c1 = (-1, 0)
            c2 = (-1, 0)
            
        return (c1, c2)
    
    def extractPerimeterEdges(self,m):     
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
    
    def makeTwist(self,a, maxTwist):
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
    
    def camber(self,ctx=None, mp=0.1, pp=0.5, weightMode=False, profileMode=False, ctype="NACA",ellipDis=False): 
    
        sce = bpy.context.scene    
        a = ctx.active_object
        
        if a is None:   
            a = bpy.data.objects["Sail"]
        for v in a.data.vertices:
            v.co.z = 0
        pe = self.extractPerimeterEdges(a)
        pv = []
        for e in pe:
            pv.append(e[0])
            pv.append(e[1])

        miny = 100000
        maxy = -100000
        for v in a.data.vertices:
            if v.co.y < miny:
                miny = v.co.y
            if v.co.y > maxy:
                maxy = v.co.y 
        halfSpan  = (maxy - miny)/2
        halfSpanY = halfSpan+miny
        print("Half span = ", halfSpan, "full=",(maxy - miny) )
        print("miny=",miny," maxy=",maxy)
        if ctype == 'NACA':            
            for v in a.data.vertices:
                if True:
                    e1, e2 = self.getEdgesCrossing(a.data.vertices, pe, v.co.y)
                    if e1[0] != -1 and e2[0] != -1:
                        x1 = self.getXinEdge(a.data.vertices, e1, v.co.y)
                        x2 = self.getXinEdge(a.data.vertices, e2, v.co.y)
                        leftx = min(x1, x2)
                        rightx = max(x1, x2)
                        x = (v.co.x - leftx) / (rightx - leftx+0.0001)
                        y = profile(x, mp, pp)
                        if weightMode:
                            if a.data.vertices[v.index].groups:
                                print("weight=",a.data.vertices[v.index].groups[0].weight)
                                a.data.vertices[v.index].co.z = a.data.vertices[v.index].groups[0].weight * y * (rightx - leftx)
                            else:    
                                a.data.vertices[v.index].co.z = y * (rightx - leftx)
                        elif profileMode:
                            weight = self.findWeight(a.matrix_world * v.co)
                            a.data.vertices[v.index].co.z = y * (rightx - leftx) * weight
                        elif ellipDis:
                            print("y=",(v.co.y-miny),"x=",x,"camber=",y)
                            print("factor=",((v.co.y-halfSpanY)/halfSpan)**2)
                            y = y*sqrt(1-((v.co.y-halfSpanY)/halfSpan)**2)
                            a.data.vertices[v.index].co.z = y * (rightx - leftx)                            
                        else:
                            a.data.vertices[v.index].co.z = y * (rightx - leftx)

        elif ctype == 'THREE':
            mp = [sce.sailflow_model.sec1M/100.0,sce.sailflow_model.sec2M/100.0,sce.sailflow_model.sec3M/100.0] 
            pp = [sce.sailflow_model.sec1P/100.0,sce.sailflow_model.sec2P/100.0,sce.sailflow_model.sec3P/100.0]
            heights = [0.0, sce.sailflow_model.sec2H/100.0, 1.0]  # Percentage, first and last 0 and 1
            for v in a.data.vertices:
    #            if not v.index in pv:
    #            if True:
                if v.co.y != 0 and v.co.x != 0:
                    e1, e2 = self.getEdgesCrossing(a.data.vertices, pe, v.co.y)
                    if e1[0] != -1 and e2[0] != -1:
                        x1 = self.getXinEdge(a.data.vertices, e1, v.co.y)
                        x2 = self.getXinEdge(a.data.vertices, e2, v.co.y)
                        leftx = min(x1, x2)
                        rightx = max(x1, x2)
                        x = (v.co.x - leftx) / (rightx - leftx+0.0001)
                        heightPerc = (v.co.y - miny) / (maxy - miny+0.0001)
                        print("HeightPerc=",heightPerc," height")
                        if heightPerc <= heights[1]:
                            rmp = mp[0]+((mp[1]-mp[0])/heights[1])*heightPerc
                            rpp = pp[0]+((pp[1]-pp[0])/heights[1])*heightPerc
                        else:
                            rmp = mp[1]+((mp[2]-mp[1])/heights[2])*(heightPerc-heights[1])    
                            rpp = pp[1]+((pp[2]-pp[1])/heights[2])*(heightPerc-heights[1])  
                        y = profile(x, rmp, rpp)
                        if ellipDis:
                            print("y=",(v.co.y-miny),"x=",x,"camber=",y)
                            print("factor=",1-((v.co.y-halfSpanY)/halfSpan)**2)
                            y = y*sqrt(1-((v.co.y-halfSpanY)/halfSpan)**2)
                            a.data.vertices[v.index].co.z = y * (rightx - leftx)                            
                        else:
                            a.data.vertices[v.index].co.z = y * (rightx - leftx)

        elif ctype == 'CUSTOM':
            print("Custom profile")
            for v in a.data.vertices:
                if not v.index in pv:
                    e1, e2 = self.getEdgesCrossing(a.data.vertices, pe, v.co.y)
                    if e1[0] != -1 and e2[0] != -1:
                        x1 = self.getXinEdge(a.data.vertices, e1, v.co.y)
                        x2 = self.getXinEdge(a.data.vertices, e2, v.co.y)
                        leftx = min(x1, x2)
                        rightx = max(x1, x2)
                        x = (v.co.x - leftx) / (rightx - leftx+0.0001)
                        heightPerc = (v.co.y - miny) / (maxy - miny)
                        y = custom_profile.profile(x, heightPerc,v.co.x,v.co.y,miny,maxy)
                        if ellipDis:
                            print("y=",(v.co.y-miny),"x=",x,"camber=",y)
                            print("factor=",1-((v.co.y-halfSpanY)/halfSpan)**2)
                            y = y*sqrt(1-((v.co.y-halfSpanY)/halfSpan)**2)
                            a.data.vertices[v.index].co.z = y * (rightx - leftx)                            
                        else:    
                            a.data.vertices[v.index].co.z = y * (rightx - leftx)
                            
        elif ctype == 'CURVE':
            print("Curve profile option")
            for v in a.data.vertices:
                if not v.index in pv:
                    e1, e2 = self.getEdgesCrossing(a.data.vertices, pe, v.co.y)
                    if e1[0] != -1 and e2[0] != -1:
                        x1 = self.getXinEdge(a.data.vertices, e1, v.co.y)
                        x2 = self.getXinEdge(a.data.vertices, e2, v.co.y)
                        leftx = min(x1, x2)
                        rightx = max(x1, x2)
                        x = (v.co.x - leftx) / (rightx - leftx+0.0001)
                        heightPerc = (v.co.y - miny) / (maxy - miny)
                        y = curveProfile(x, sce.sailflow_model.curvePoints)
                        if ellipDis:
                            print("y=",(v.co.y-miny),"x=",x,"camber=",y)
                            print("factor=",1-((v.co.y-halfSpanY)/halfSpan)**2)
                            y = y*sqrt(1-((v.co.y-halfSpanY)/halfSpan)**2)
                            a.data.vertices[v.index].co.z = y * (rightx - leftx)                            
                        else:    
                            a.data.vertices[v.index].co.z = y * (rightx - leftx)
              
            
    def execute(self, context):
        print("Called airprofile")    
        sce = bpy.context.scene
        self.camber(context, sce.sailflow_model.m / 100, sce.sailflow_model.p / 100, sce.sailflow_model.weight, 
                    sce.sailflow_model.curve, sce.sailflow_model.t,sce.sailflow_model.ellipDis)
        if sce.sailflow_model.twist:
            self.makeTwist(context.active_object, sce.sailflow_model.tw)
        return {'FINISHED'}
        
    @classmethod
    def poll(cls, context):
        return context.active_object and context.active_object.type == 'MESH'
     
class printPDF(bpy.types.Operator):
    bl_idname = "mesh.print_pdf"
    bl_label = "Print PDF"
    bl_description = "print PDF of selected panels"

    unitToMm = 1000

    
    filepath = bpy.props.StringProperty(subtype="FILE_PATH")
  
    def containedIn(self,v,w,h):
        return (v.x >= 0 and v.x <= w and v.y >= 0 and v.y <= h)
          
    def clipping(self,v1,v2,w,h):
        vout1 = None
        vout2 = None
        if self.containedIn(v1,w,h):
            #print(" clipping: V1 already in ")
            vout1 = v1
        if self.containedIn(v2,w,h):
            #print(" clipping: V2 already in ")
            vout1 = v2

        v = mathutils.geometry.intersect_line_line_2d(v1,v2,Vector((0,0,0)),Vector((w,0,0)))
        if v:
            #print(" clipping: h bassa")
            if vout1 == None:
                vout1 = v
            elif (v != vout1):
                vout2 = v
                return vout1,vout2
        v = mathutils.geometry.intersect_line_line_2d(v1,v2,Vector((w,0,0)),Vector((w,h,0)))
        if v:
            #print(" clipping: v destra")
            if vout1 == None:
                vout1 = v
            elif (v != vout1):
                vout2 = v
                return vout1,vout2
        v = mathutils.geometry.intersect_line_line_2d(v1,v2,Vector((w,h,0)),Vector((0,h,0)))
        if v:
            #print(" clipping: h alta")
            if vout1 == None:
                vout1 = v
            elif (v != vout1):
                vout2 = v
                return vout1,vout2
        v = mathutils.geometry.intersect_line_line_2d(v1,v2,Vector((0,h,0)),Vector((0,0,0)))
        if v:
            #print(" clipping: v sinistra")            
            if vout1 == None:
                vout1 = v
            elif (v != vout1):
                vout2 = v
        return vout1,vout2
        
    def makePDF(self, a, pe, vxs, fmt, freeTxt='', mp=False,offSetSet=0,overlapSet=0,w=0,h=0):
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

        #print("Vector")
        #for i in range(len(vxs_pix)):
            #print(i,vxs_pix[i])
            
        fname = self.filepath   
        pdf = FPDF(format=fmt,pw=w,ph=h) 
        offSet = pdf.dimension()[0]*offSetSet/100 
        #percentage of the offSet used for overlap between pages
        overlap = pdf.dimension()[0]*overlapSet/100
        pdf.set_compression(False)
        if not mp:
            pdf.add_page()
            for e in pe:
                p1 = vxs_pix[e[0]]
                p2 = vxs_pix[e[1]]
                pdf.line(p1.x + offSet, (pdf.dimension()[1]-p1.y)-offSet , p2.x + offSet, (pdf.dimension()[1]-p2.y)-offSet)
            
            # Make the header
            pdf.set_x(0)
            pdf.set_y(0)
            date = str(time.localtime().tm_mday) + "/" + str(time.localtime().tm_mon) + "/" + str(time.localtime().tm_year)
            pdf.set_font('Arial', 'B', 12)
            #pdf.cell(w=75, h=10, txt=a.name, border=1, ln=2, align='C')
            #pdf.set_font(family='Arial', size=12)    
            pdf.cell(w=75, h=5, txt=date, border=1, ln=2, align='C')
            #pdf.cell(w=75, h=5, txt=freeTxt + " ", border=1, ln=0, align='C')
        else:
            clipSizeH = int(pdf.dimension()[0]-2*offSet) 
            clipSizeV = int(pdf.dimension()[1]-2*offSet)  
            minx = min([v.x for v in vxs_pix])
            maxx = max([v.x for v in vxs_pix])
            miny = min([v.y for v in vxs_pix])
            maxy = max([v.y for v in vxs_pix])
            numPagesH = int((maxx - minx)/(clipSizeH-overlap))+1
            numPagesV = int((maxy - miny)/(clipSizeV-overlap))+1  
 
            for h in range(0,numPagesH):
                for v in range(0,numPagesV):
                    pdf.add_page()
                    line=0
                    for e in pe:
                        # Translate the points
                        cross1 = Vector([vxs_pix[e[0]].x,vxs_pix[e[0]].y])
                        cross2 = Vector([vxs_pix[e[1]].x,vxs_pix[e[1]].y])
                        cross1.x = cross1.x - float(clipSizeH*h)+float(overlap*h)
                        cross2.x = cross2.x - float(clipSizeH*h)+float(overlap*h)
                        cross1.y = cross1.y - float(clipSizeV*v)+float(overlap*v)
                        cross2.y = cross2.y - float(clipSizeV*v)+float(overlap*v)

                        if not (self.containedIn(cross1,clipSizeH, clipSizeV) and
                                self.containedIn(cross2,clipSizeH, clipSizeV)):
                            cross1, cross2 = self.clipping(cross1,cross2,clipSizeH, clipSizeV)
 
                        if cross1 != None and cross2 != None:
                            pdf.line(cross1.x+offSet,
                                     (pdf.dimension()[1]-cross1.y)-offSet,
                                     cross2.x+offSet,
                                     (pdf.dimension()[1]-cross2.y)-offSet)
                            line +=1
                    
                    if (h+v!=0):                 
                        pdf.set_line_width(0.2)        
                        pdf.line(0.0,pdf.dimension()[1]-offSet-overlap,pdf.dimension()[0],pdf.dimension()[1]-offSet-overlap)
                        pdf.line(offSet+overlap,0.0,offSet+overlap,pdf.dimension()[1] )
                    
                    pdf.set_x(0)
                    pdf.set_y(0)   
                    #date = str(time.localtime().tm_mday) + "/" + str(time.localtime().tm_mon) + "/" + str(time.localtime().tm_year)
                    #pdf.set_font('Arial', 'B', 24)
                    pdf.set_font(family='Arial', size=12)    
                    pdf.cell(w=75, h=10, txt=a.name+"H:"+str(h)+" V:"+str(v), border=1, ln=2, align='C')
                    #pdf.cell(w=75, h=5, txt=date, border=1, ln=2, align='C')
                    #pdf.cell(w=75, h=5, txt=freeTxt + " ", border=1, ln=0, align='C')
                     
                    
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
        self.makePDF(a, el, vxs, se.paperFormat, se.freeText,se.multiPages,se.margin,se.overlap,se.paperWidth,se.paperHeight);  
        return {'FINISHED'}
    
    def invoke(self, context,event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
    
    @classmethod
    def poll(cls, context):
        return context.active_object and context.active_object.type == 'MESH'                                          

class AirFoilSettings(bpy.types.PropertyGroup):

    curveTypes = [
      ("NACA", "Naca 4", "", 1),
      ("THREE", "Three Sections", "", 2),
      ("CUSTOM", "Custom Profile", "", 3),
      ("CURVE" , "Curve", "", 4),
      ]
    m = IntProperty (
        name="% Max Camber",
        description="Maximum Camber (% of chord)",
        default=0,
        min=0,
        max=50)
    p = IntProperty (
        name="% Camber Pos",
        description="Position of Camber (% of chord)",
        default=0,
        min=20,
        max=70)   
    
    t = EnumProperty(name="Curve",default="NACA",items=curveTypes)
    
    weight = BoolProperty(name="Apply Weight",default=False)
    curve = BoolProperty(name="Apply Curve",default=False)
    twist = BoolProperty(name="Apply Twist",default=False)
    tw = IntProperty (name="Twist Angle",description="Value of angle", default=0,min=0,max=90)
    energyMinimizer = bpy.props.BoolProperty(name="Stress Relief", default=False)
    maxDeformation = bpy.props.FloatProperty(name="Min stress reduction", default=0.001, min=0.000001, max=0.1)
    deltaDeformation = bpy.props.IntProperty(name="accuracy", default=3, min=1, max=10)
    polSeed = bpy.props.IntProperty(name="Use start face",default=-1)
    useSeed =  bpy.props.BoolProperty(name="Start from face", default=False)
    resEnergy = bpy.props.FloatProperty(name="Residual Strain",default=0)
    
    #paperFormat = bpy.props.StringProperty (name="Paper Size",description="others,4a0,2a0,a0,a1,a2,a3,a4",default='a0')
    freeText = bpy.props.StringProperty (name="Free Text", description="Anything appearing in the PDF",default='free text')
    margin = IntProperty (name="% of page margin",
                          default=10,
                          min=0,
                          max=25)
    overlap = IntProperty (name="% of overlap bw pages",
                          default=10,
                          min=0,
                          max=25)
    paperWidth = IntProperty(name="Width",max=1000)
    paperHeight = IntProperty(name="Height", max=3500)
    
    ellipDis = BoolProperty(name="Eliptic Distribution", default=False)
    
    sec1M = IntProperty(name="Section 1 % Camber percentage",min=0,max=70)
    sec1P = IntProperty(name="Section 1 % Camber position",min=0,max=70)
    sec2M = IntProperty(name="Section 2 % Camber percentage",min=0,max=70)
    sec2P = IntProperty(name="Section 2 % Camber position",min=0,max=70)
    sec2H = IntProperty(name="Section 2 % Position Height",min=0,max=70)
    sec3M = IntProperty(name="Section 3 % Camber percentage",min=0,max=70)
    sec3P = IntProperty(name="Section 3 % Camber position",min=0,max=70)
    
    paperSizes = [
      ("4A0","4a0","",1),
      ("2A0","2a0","",2),
      ("A0","a0","",3),
      ("A1","a1","",4),
      ("A2","a2","",5),
      ("A3","a3","",6),
      ("A4","a4","",7),
      ("Other","other","",8)]
    
    paperFormat = bpy.props.EnumProperty (name = "Paper Size", default="A4",items=paperSizes)   
    freeText = bpy.props.StringProperty (name = "Free Text", description = "Anything appearing in the PDF", default = 'free text')
    multiPages = bpy.props.BoolProperty(name="Multi pages",default=False)
    panel1 = []
    panel2 = []
    panel3 = []
    panel4 = []
    panel5 = []
    panel6 = []
    panel7 = []
    panel8 = []
    panel9 = []
    panel10 = []
    panelNumber = IntProperty(name="Panel Number",min=1,max=10)
    
    curvePoints = []
                                    
def register():
    bpy.utils.register_module(__name__)
    bpy.types.Scene.sailflow_model = bpy.props.PointerProperty(type=AirFoilSettings,
                                        name="Airfoil Model",
                                        description="Setting of the AirFoil")

def unregister():
    bpy.utils.unregister_module(__name__)
    
if __name__ == "__main__":
    register()
