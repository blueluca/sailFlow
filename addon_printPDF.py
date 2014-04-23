bl_info = {
    "name": "Print PDF",
    "description": "Export PDF",
    "author": "Luca Stagnaro",
    "version": (1, 0, 0),
    "blender": (2, 6, 9),
    "api": 33411, #Not certain on the API version
    "location": "View3D > Object > Print PDF",
    "warning": "",
    "category": "Import-Export"}


import bpy
import sys
import mathutils
import math
from fpdf import FPDF
import time

from mathutils import Vector, Euler    
from bpy.types import Operator
from bpy.props import IntProperty,EnumProperty,BoolProperty,StringProperty

unitToMm = 1000
offSet = 20  

def mainPart(a,pe,vxs,fmt,freeTxt=''):
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
        vc = a.matrix_world*v.co        
        vxs_pix.append(Vector([(vc.x-min_x)*unitToMm,(vc.y-min_y)*unitToMm]))
    
    fname = '/Users/lstagnar/Desktop/'+a.name+'.pdf'    
    pdf=FPDF(format=fmt)
    pdf.set_compression(False)
    pdf.add_page()
#    pdf.set_font('Arial','B',16)
#    pdf.cell(40,10,'Hello World!')
    #pdf.cell(40,10,a.name)
    for e in pe:
        p1 = vxs_pix[e[0]]
        p2 = vxs_pix[e[1]]
        print(fmt,p1,p2)
        pdf.line(p1.x+offSet,p1.y+offSet,p2.x+offSet,p2.y+offSet)
    
    # Make the header
    pdf.set_x(0)
    pdf.set_y(0)
    date = str(time.localtime().tm_mday)+"/"+str(time.localtime().tm_mon)+"/"+str(time.localtime().tm_year)
    pdf.set_font('Arial','B',24)
    pdf.cell(w=75,h=10,txt=a.name,border=1,ln=2,align='C')
    pdf.set_font(family='Arial',size=16)    
    pdf.cell(w=75,h=5,txt=date,border=1,ln=2,align='C')
    pdf.cell(w=75,h=5,txt=freeTxt+" ",border=1,ln=0,align='C')

    pdf.output(fname,'F')
            
class printPDF(bpy.types.Operator):
    """ PrintSVG Operator"""
    bl_idname = "object.print_pdf"
    bl_label = "Print PDF"
    f = bpy.props.StringProperty (
        name = "Paper Size",
        description = "4a0,2a0,a0,a1,a2,a3,a4",
        default = 'a0')
    t = bpy.props.StringProperty (
        name = "Free Text",
        description = "Anything appearing in the PDF",
        default = 'free text')
     
    def execute(self, context):  
#        bpy.ops.object.dialog_svg('INVOKE_DEFAULT')
        a = context.active_object
        el = []
        for p in a.data.polygons:
            for e in p.edge_keys:
                el.append((e[1],e[0]))
        vxs = a.data.vertices
        mainPart(a,el,vxs,self.f,self.t)
        return {'FINISHED'}
    
    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)        

def add_object_button(self,context):
    self.layout.operator(printPDF.bl_idname,text="Print in PDF")

def register():
    bpy.utils.register_class(printPDF)
    bpy.types.VIEW3D_MT_object.append(add_object_button)

def unregister():
    bpy.utils.unregister_class(printPDF)
    bpy.types.VIEW3D_MT_object.remove(add_object_button)
    
if __name__ == "__main__":
    a = bpy.context.active_object
    el = []
    for p in a.data.polygons:
        for e in p.edge_keys:
            el.append((e[1],e[0]))
    vxs = a.data.vertices
    mainPart(a,el,vxs,'a4',freeTxt='Ciccio Ciccio')
    #register()
    