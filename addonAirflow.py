bl_info = {
    "name": "Create Sailprofile",
    "description": "Creates a profile for a sail",
    "author": "Luca Stagnaro",
    "version": (0, 1, 0),
    "blender": (2, 6, 9),
    "api": 33411,  # Not certain on the API version
    "location": "View3D > Add > Mesh > Airfoil",
    "warning": "",
    "category": "Add Mesh"}

import bpy
import sys
import mathutils
import math

from math import *
from mathutils import Vector, Euler,geometry
  
from bpy.types import Operator
from bpy.props import IntProperty, EnumProperty, BoolProperty, StringProperty

from fpdf import FPDF
import time
import imp

PYDEV_SOURCE_DIR = 'C:/eclipse/plugins/org.python.pydev_3.4.1.201403181715/pysrc'
  
if sys.path.count(PYDEV_SOURCE_DIR) < 1:
   sys.path.append(PYDEV_SOURCE_DIR)
import pydevd
  

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
     
    def minimizeEnergy(self, F, maxDeformation, deltaDeformation):
        global Vxs
        
      #  print("Energy minimizer ", maxDeformation)
        MAXCOUNT = 5000
        delta = 0.1
        while delta > 10**-deltaDeformation:
            print("---------------> Loop minimize delta=",delta)
            for count in range(MAXCOUNT):
                nodesMovToGain = [0 for x in range(len(Vxs))]            
                for vix in range(len(Vxs)): 
                    # get the triangles with the vertex vix in common
                    C = []
                    energy = 0
                    for p in F:
                        if vix in p.vertices:
                            C.append(p)
                    # C now contains the list of polys with the vertix index vix
                    vandl = []
                    for c in C:
                        l = c.adjacentVx(vix)
                        if not l[0] in vandl:
                            vandl.append(l[0])
                        if not l[1] in vandl:
                            vandl.append(l[1])
                    #vandl contains the list of points and original length 
                    #of the verctices not vix
                    # energy accumulate the difference of energy
                    for v in vandl: 
                        dl = (Vxs[v[0]].co - Vxs[vix].co).length - v[1]
                        #dl = difference between original length and new length
                        energy = energy + dl*dl/ v[1]
                    originalCo = Vxs[vix].co
                    overallpdx = 0.0
                    overallmdx = 0.0
                    overallpdy = 0.0
                    overallmdy = 0.0
                    Vxs[vix].co.x = originalCo.x + delta
                    for v in vandl:
                        dl = (Vxs[v[0]].co - Vxs[vix].co).length - v[1]
                        overallpdx = overallpdx + dl* dl / v[1]
                        # move the vertix by a -dx and calculate the overallmdx
                    Vxs[vix].co.x = originalCo.x - delta
                    for v in vandl:
                        dl = (Vxs[v[0]].co - Vxs[vix].co).length - v[1]
                        overallmdx = overallmdx + dl* dl / v[1]
                        # move the vertix by a +dy and calculate the overallpdy
                    Vxs[vix].co.x = originalCo.x
                    Vxs[vix].co.y = originalCo.y + delta
                    for v in vandl:
                        dl = (Vxs[v[0]].co - Vxs[vix].co).length - v[1]
                        overallpdy = overallpdy + dl*dl / v[1]
                        # move the vertix by a -dy and calculate the overallpdx
                    Vxs[vix].co.y = originalCo.y - delta
                    for v in vandl:
                        dl = (Vxs[v[0]].co - Vxs[vix].co).length - v[1]
                        overallmdy = overallmdy + dl*dl / v[1]
                    Vxs[vix].co = originalCo
                    # the overallxxx are the energy with a small change in length
                    if energy > min(overallpdx, overallmdx, overallpdy, overallmdy):
                        if overallpdx == min(overallpdx, overallmdx, overallpdy, overallmdy):
                            nodesMovToGain[vix] = [delta,0,energy - overallpdx]
                        elif overallmdx == min(overallpdx, overallmdx, overallpdy, overallmdy):
                            nodesMovToGain[vix] = [-delta,0,energy - overallmdx]
                        elif overallpdy == min(overallpdx, overallmdx, overallpdy, overallmdy):
                            nodesMovToGain[vix] = [0,delta,energy - overallpdy]
                        else:
                            nodesMovToGain[vix] = [0,-delta,energy - overallmdy]
                    else:
                            nodesMovToGain[vix] = [0,0,0.0]
                
                maxGain = maxDeformation
                maxIdx = None
                #search of the node with the max gain in term of energy
                for n in range(len(nodesMovToGain)):
                    if nodesMovToGain[n][2] > maxGain:
                        maxGain = nodesMovToGain[n][2]
                        maxIdx = n
                if maxIdx:
#                    print("%d) Max reduction for node %d, Denergy %.6f coord %2.4f %2.4f"%(count,maxIdx,nodesMovToGain[maxIdx][2],Vxs[maxIdx].co.x,Vxs[maxIdx].co.y))
                    Vxs[maxIdx].co.x += nodesMovToGain[maxIdx][0]
                    Vxs[maxIdx].co.y += nodesMovToGain[maxIdx][1]
                else:
#                    print("No more gain, residual energy:", energy)
                    break
            delta = delta / 10
          
    def makeItFlat(self, obj, energyMinimizer, maxDeformation, deltaDeformation):
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
        
#        pydevd.settrace(stdoutToServer=True, stderrToServer=True, suspend=True)
        s = None
        for p in me.polygons:
            if p.select:
                vt = MyPoly(p, idx=p.index)
                V.append(vt)
         
        while V or A:
            if A:
                s = A.pop()
            else:
                s = V.pop()
            # Collect all the adjacent triangle
            # and put it in Active collection "A"
            found = True
            while(found):
                at = self.findAdjacentNonFlat(V, s)
                if at:
                    # debug("..adj is"+str(at))
                    A.append(at)
                    V.remove(at) 
                    found = True
                else:
                    found = False
            # debug("Flattening: " + str(s))
            s.flatten()
            F.append(s)
    
        if energyMinimizer:
            self.minimizeEnergy(F, maxDeformation, deltaDeformation)
       
    @classmethod
    def poll(cls, context):
        return context.active_object and context.active_object.type == 'MESH'
    
    def execute(self, context): 
        global Vxs
        global F
        

        sce = bpy.context.scene      
        if F:
            del F[:]
            
        self.makeItFlat(bpy.context.active_object, sce.airflow_model.energyMinimizer, sce.airflow_model.maxDeformation, sce.airflow_model.deltaDeformation)
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

class VIEW3D_PT_airprofile_print(bpy.types.Panel):
    bl_label = "PDF Generation"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Airflow Design"
    
    def draw(self, context):
        sce = bpy.context.scene
        layout = self.layout
        col = layout.column(align=True)

        col.prop(sce.airflow_model, "paperFormat")
        col.prop(sce.airflow_model, "freeText")
        col.prop(sce.airflow_model, "multiPages")
        
        col = layout.column(align=True)
        col.operator("mesh.print_pdf")
        
class VIEW3D_PT_airprofile_parameters(bpy.types.Panel):
    bl_label = "Parameters"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Airflow Design"
    
    def draw(self, context):
        sce = bpy.context.scene
        layout = self.layout
        col = layout.column(align=True)

        if sce.airflow_model.t == 'NACA':
            col.prop(sce.airflow_model, "m")
            col.prop(sce.airflow_model, "p")
             
        elif sce.airflow_model.t == "THREE":
            col = layout.column(align=True)
            col.label(text="Section 1:")
            col.prop(sce.airflow_model, "sec1M")
            col.prop(sce.airflow_model, "sec1P")

            col = layout.column(align=True)
            col.label(text="Section 2:")
            col.prop(sce.airflow_model, "sec2M")
            col.prop(sce.airflow_model, "sec2P")
            col.prop(sce.airflow_model, "sec2H")

            col = layout.column(align=True)            
            col.label(text="Section 3:")
            col.prop(sce.airflow_model, "sec3M")
            col.prop(sce.airflow_model, "sec3P")
            
        elif sce.airflow_model.t == "CUSTOM":
            col.operator("mesh.load_library")
            
        col = layout.column(align=True)
        col.prop(sce.airflow_model, "t")
        col.prop(sce.airflow_model, "weight")
        
        col = layout.column(align=True)
        col.label(text="Twist:")
        row = col.row(align=True)
        row.prop(sce.airflow_model, "twist")
        row.prop(sce.airflow_model, "tw")

        col = layout.column(align=True)    
        col.label(text="Operators:")
        col.operator("mesh.airprof")
        col.prop(sce.airflow_model,"energyMinimizer")
        col.prop(sce.airflow_model,"deltaDeformation")
        col.prop(sce.airflow_model,"maxDeformation") 
        col = layout.column(align=True)    
        col.operator("mesh.flattener")
            
        if sce.airflow_model.t == 'NACA':
            box = layout.box()
            angleIn = atan(2 * sce.airflow_model.m / sce.airflow_model.p) * 180 / 3.14159
            angleOut = -atan(2 * sce.airflow_model.m / 100 / (sce.airflow_model.p / 100 - 1)) * 180 / 3.14159
            box.label(text="Angle in  " + str(round(angleIn, 2)))
            box.label(text="Angle out " + str(round(angleOut, 2)))

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
    #    x,y = lineIntersection2D(x1,y1,x2,y2,-10,yc,+10,yc) 
        v = mathutils.geometry.intersect_line_line_2d(Vector((x1, y1, 0)),
                                                      Vector((x2, y2, 0)),
                                                      Vector((-10, yc, 0)),
                                                      Vector((+10, yc, 0)))
        if v: 
            return v.x
        else:
            return 0.0    
        
    def getEdgesCrossing(self,vl, pe, y):
        c1 = (-1, 0)
        c2 = (0, 0)
        for e in pe:
            if vl[e[0]].co.y >= y and vl[e[1]].co.y <= y:
                if c1[0] == -1:
                    c1 = e
                else:
                    if not ((c1[0] in e) or (c1[1] in e)):
                        c2 = e
            elif vl[e[0]].co.y <= y and vl[e[1]].co.y >= y:
                if c1[0] == -1:
                    c1 = e
                else:
                    if not ((c1[0] in e) or (c1[1] in e)):
                        c2 = e
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
                       # print("X1 %2.2f  X2 %2.2f\n"%(x1,x2))
                        leftx = min(x1, x2)
                        angle = ((v.co.y - miny) / (maxy - miny) * maxTwist) / 180 * 3.14159
                        a.data.vertices[v.index].co.z += sin(angle) * (v.co.x - leftx)
                        a.data.vertices[v.index].co.x = cos(angle) * (v.co.x - leftx) + leftx   

    def camber(self,ctx=None, mp=0.1, pp=0.5, weightMode=False, profileMode=False, ctype="NACA"): 
    
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

        if ctype == 'NACA':            
            for v in a.data.vertices:
#                if not (v.index in pv) or True:
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
                                a.data.vertices[v.index].co.z = a.data.vertices[v.index].groups[0].weight * y * (rightx - leftx)
                            else:    
                                a.data.vertices[v.index].co.z = y * (rightx - leftx)
                        elif profileMode:
                            weight = self.findWeight(a.matrix_world * v.co)
                            a.data.vertices[v.index].co.z = y * (rightx - leftx) * weight
                        else:
                            a.data.vertices[v.index].co.z = y * (rightx - leftx)
        elif ctype == 'THREE':
            miny = 100000
            maxy = -100000
            for v in a.data.vertices:
                if v.co.y < miny:
                    miny = v.co.y
                if v.co.y > maxy:
                    maxy = v.co.y 
            mp = [sce.airflow_model.sec1M/100.0,sce.airflow_model.sec2M/100.0,sce.airflow_model.sec3M/100.0] 
            pp = [sce.airflow_model.sec1P/100.0,sce.airflow_model.sec2P/100.0,sce.airflow_model.sec3P/100.0]
            heights = [0.0, sce.airflow_model.sec2H/100.0, 1.0]  # Percentage, first and last 0 and 1
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
                        x = (v.co.x - leftx) / (rightx - leftx)
                        heightPerc = (v.co.y - miny) / (maxy - miny)
                        print("HeightPerc=",heightPerc," height")
                        if heightPerc <= heights[1]:
                            rmp = mp[0]+((mp[1]-mp[0])/heights[1])*heightPerc
                            rpp = pp[0]+((pp[1]-pp[0])/heights[1])*heightPerc
                        else:
                            rmp = mp[1]+((mp[2]-mp[1])/heights[2])*(heightPerc-heights[1])    
                            rpp = pp[1]+((pp[2]-pp[1])/heights[2])*(heightPerc-heights[1])  
                        y = profile(x, rmp, rpp)
                        a.data.vertices[v.index].co.z = y * (rightx - leftx) 
        elif ctype == 'CUSTOM':
            print("Custom profile")
            miny = 100000
            maxy = -100000
            for v in a.data.vertices:
                if v.co.y < miny:
                    miny = v.co.y
                if v.co.y > maxy:
                    maxy = v.co.y 
            for v in a.data.vertices:
                if not v.index in pv:
                    e1, e2 = self.getEdgesCrossing(a.data.vertices, pe, v.co.y)
                    if e1[0] != -1 and e2[0] != -1:
                        x1 = self.getXinEdge(a.data.vertices, e1, v.co.y)
                        x2 = self.getXinEdge(a.data.vertices, e2, v.co.y)
                        leftx = min(x1, x2)
                        rightx = max(x1, x2)
                        x = (v.co.x - leftx) / (rightx - leftx)
                        heightPerc = (v.co.y - miny) / (maxy - miny)
                        y = custom_profile.profile(x, heightPerc)
                        a.data.vertices[v.index].co.z = y * (rightx - leftx)
            
            
            
    def execute(self, context):
        print("Called airprofile")    
        sce = bpy.context.scene    
        self.camber(context, sce.airflow_model.m / 100, sce.airflow_model.p / 100, sce.airflow_model.weight, sce.airflow_model.curve, sce.airflow_model.t)
        if sce.airflow_model.twist:
            self.makeTwist(context.active_object, sce.airflow_model.tw)
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
        
    def makePDF(self, a, pe, vxs, fmt, freeTxt='', mp=False):
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

        print("Vector")
        for i in range(len(vxs_pix)):
            print(i,vxs_pix[i])
            
        fname = self.filepath   
        pdf = FPDF(format=fmt)
        self.offSet = pdf.dimension()[0]*0.01 
        pdf.set_compression(False)
        if not mp:
            pdf.add_page()
            for e in pe:
                p1 = vxs_pix[e[0]]
                p2 = vxs_pix[e[1]]
                pdf.line(p1.x + self.offSet, (pdf.dimension()[1]-p1.y)-self.offSet , p2.x + self.offSet, (pdf.dimension()[1]-p2.y)-self.offSet)
            
            # Make the header
            pdf.set_x(0)
            pdf.set_y(0)
            date = str(time.localtime().tm_mday) + "/" + str(time.localtime().tm_mon) + "/" + str(time.localtime().tm_year)
            pdf.set_font('Arial', 'B', 24)
            pdf.cell(w=75, h=10, txt=a.name, border=1, ln=2, align='C')
            pdf.set_font(family='Arial', size=16)    
            pdf.cell(w=75, h=5, txt=date, border=1, ln=2, align='C')
            pdf.cell(w=75, h=5, txt=freeTxt + " ", border=1, ln=0, align='C')
        else:
            maxx = -100000
            maxy = -100000
            miny = 100000
            minx = 100000
            clipSizeH = int(pdf.dimension()[0]-2*self.offSet) 
            clipSizeV = int(pdf.dimension()[1]-2*self.offSet)  
            for v in vxs_pix:
                if minx > v.x:
                    minx = v.x
                if maxx < v.x:
                    maxx = v.x
                if maxy < v.y:
                    maxy = v.y
                if miny > v.y:
                    miny = v.y
            numPagesH = int((maxx - minx)/clipSizeH)+1
            numPagesV = int((maxy - miny)/clipSizeV)+1  
            print("-"*80)
            print("Num pages H & V,",clipSizeH,clipSizeV," pages -->",numPagesH,numPagesV)      
            print("-"*80)
 
            for h in range(0,numPagesH):
                for v in range(0,numPagesV):
                    pdf.add_page()
                    print("="*40)
                    print("New page h=%d v=%d"%(h,v))
                    print("="*40)
                    line=0
                    for e in pe:
                        # Translate the points
                        cross1 = Vector([vxs_pix[e[0]].x,vxs_pix[e[0]].y])
                        cross2 = Vector([vxs_pix[e[1]].x,vxs_pix[e[1]].y])
                        cross1.x = cross1.x - float(clipSizeH*h)
                        cross2.x = cross2.x - float(clipSizeH*h)
                        cross1.y = cross1.y - float(clipSizeV*v)
                        cross2.y = cross2.y - float(clipSizeV*v)
                                                
                        if (self.containedIn(cross1,clipSizeH, clipSizeV) and
                            self.containedIn(cross2,clipSizeH, clipSizeV)):
                            #print("  All contained",e[0],e[1],cross1,cross2)
                            pdf.line(cross1.x+self.offSet,(pdf.dimension()[1]-cross1.y)-self.offSet,
                                     cross2.x+self.offSet,(pdf.dimension()[1]-cross2.y)-self.offSet)
                            line +=1
                            
                        else:
                            #print("  clipping",e[0],e[1],cross1,cross2)
                            cross1, cross2 = self.clipping(cross1,cross2,clipSizeH, clipSizeV)
                            #print("  return clipping",e[0],e[1],cross1,cross2)
                            if cross1 != None and cross2 != None:
                                pdf.line(cross1.x+self.offSet,(pdf.dimension()[1]-cross1.y)-self.offSet,
                                         cross2.x+self.offSet,(pdf.dimension()[1]-cross2.y)-self.offSet)
                                line +=1
                                
                    #pdf.set_line_width(0.2)        
                    #pdf.rect(self.offSet/2,self.offSet/2,pdf.dimension()[0]-2*self.offSet/2,pdf.dimension()[1]-2*self.offSet/2)

                    pdf.set_x(0)
                    pdf.set_y(0)   
                    date = str(time.localtime().tm_mday) + "/" + str(time.localtime().tm_mon) + "/" + str(time.localtime().tm_year)
                    pdf.set_font('Arial', 'B', 24)
                    pdf.cell(w=75, h=10, txt=a.name+"V:"+str(h)+" H:"+str(v), border=1, ln=2, align='C')
                    pdf.set_font(family='Arial', size=16)    
                    pdf.cell(w=75, h=5, txt=date, border=1, ln=2, align='C')
                    pdf.cell(w=75, h=5, txt=freeTxt + " ", border=1, ln=0, align='C')
                     
                    
        pdf.output(fname, 'F')    
     
    def execute(self, context):
        sce = bpy.context.scene    
        a = context.active_object
        el = []
        for p in a.data.polygons:
            for e in p.edge_keys:
                el.append((e[1], e[0]))
        vxs = a.data.vertices
        self.makePDF(a, el, vxs, sce.airflow_model.paperFormat, sce.airflow_model.freeText,sce.airflow_model.multiPages);  
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
      ("CUSTOM", "Custom Profile", "", 3)
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
    paperFormat = bpy.props.StringProperty (name="Paper Size",description="4a0,2a0,a0,a1,a2,a3,a4",default='a0')
    freeText = bpy.props.StringProperty (name="Free Text", description="Anything appearing in the PDF",default='free text')
    
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
      ("A4","a4","",7)]
    paperFormat = bpy.props.EnumProperty (name = "Paper Size", default="A4",items=paperSizes)   
    freeText = bpy.props.StringProperty (name = "Free Text", description = "Anything appearing in the PDF", default = 'free text')
    multiPages = bpy.props.BoolProperty(name="Multi pages",default=False)
                                    
def register():
    bpy.utils.register_module(__name__)
    bpy.types.Scene.airflow_model = bpy.props.PointerProperty(type=AirFoilSettings,
                                        name="Airfoil Model",
                                        description="Setting of the AirFoil",
                                        options={'SKIP_SAVE'})

def unregister():
    bpy.utils.unregister_module(__name__)
    
if __name__ == "__main__":
    register()
