import bpy
import bmesh
from math import radians,degrees
from mathutils import *
import mathutils.geometry as mg

def flatFace(verts):
    #print("    Flatface: original vertices %s %s %s"%(verts[0],verts[1],verts[2]))
    quatr = mg.normal(verts).rotation_difference(Vector((0, 0, 1)))
    eul = quatr.to_euler()
    eul.z = 0.0
    #print("    Flatface called Euler is x=%2.2g deg y=%2.2g deg"%(degrees(eul.x),degrees(eul.y)))
    for v in verts:
        v.rotate(eul)
        v.z = 0
    #print("    Flatface: rotated  vertices %s %s %s"%(verts[0],verts[1],verts[2]))

def faceTranslate(dv,verts):
    for i,v in enumerate(verts):
        verts[i] = v + dv

def addFace(v,f,verts,faceIndex,bm_verts):
    global facesDictionary
    global vertsDictionary

    if bm_verts[0] in vertsDictionary:
        v1 = vertsDictionary[bm_verts[0]]
    else:
        v.append(verts[0])
        v1 = len(v)-1
        vertsDictionary[bm_verts[0]] = v1
    if bm_verts[1] in vertsDictionary:
        v2 = vertsDictionary[bm_verts[1]]
    else:
        v.append(verts[1])
        v2 = len(v)-1
        vertsDictionary[bm_verts[1]] = v2
    if bm_verts[2] in vertsDictionary:
        v3 = vertsDictionary[bm_verts[2]]
    else:
        v.append(verts[2])
        v3 = len(v)-1
        vertsDictionary[bm_verts[2]] = v3
    f.append((v1, v2, v3))
    facesDictionary[faceIndex] = len(f)-1

def addDEBUGFace(v,f,verts,faceIndex,bm_verts):
    global facesDictionary
    global vertsDictionary

#    print("    adding one face for face %d indexs %d %d %d "%(faceIndex,bm_verts[0],bm_verts[1],bm_verts[2]))
    v.append(verts[0])
    v1 = len(v)-1
    vertsDictionary[bm_verts[0]] = v1

    v.append(verts[1])
    v2 = len(v)-1
    vertsDictionary[bm_verts[1]] = v2

    v.append(verts[2])
    v3 = len(v)-1
    vertsDictionary[bm_verts[2]] = v3

    f.append((v1, v2, v3))
#    print("    added ",f[len(f)-1])
    facesDictionary[faceIndex] = len(f)-1
#    print("    Area new  face %g" % (mg.area_tri(verts[0], verts[1], verts[2])))

bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
vertices = []
faces = []
flatted = []
facesDictionary = {}
vertsDictionary = {}

print("-"*60)
print("Faces")

bfaces = [f for f in bm.faces if f.select]

# lets flat the first face
f = bfaces[0]
fverts = [Vector(v.co) for v in f.verts]
flatFace(fverts)
addFace(vertices, faces, fverts, f.index, [f.verts[0].index,f.verts[1].index,f.verts[2].index])
flatted.append(f.index)
refFace = []
lfverts = [Vector((0,0,0)),Vector((0,0,0)),Vector((0,0,0))]


while(bfaces):
    f = bfaces.pop(0)
    refFace.append(f)
    while (refFace):
        f = refFace.pop()
        #print("Handling face:", f.index)
        # go through the edges
        managedEdge = False
        for e in [ e for e in f.edges if len(e.link_faces)>1]:
            #print(" +-edge %d with edges %d %d"%(e.index,e.verts[0].index,e.verts[1].index))
            v1 = e.verts[0]
            v2 = e.verts[1]

            # get the adjacent face
            if e.link_faces[0] == f:
                lf = e.link_faces[1]
            else:
                lf  = e.link_faces[0]
            if lf.index not in flatted and lf.select:
                #print("  +-adj face %d vertices %d %d %d"%(lf.index,lf.verts[0].index,lf.verts[1].index,lf.verts[2].index))
                v = [v.index for v in lf.verts]
                #print("    list of index %s"%(v))
                i1 = v.index(v1.index)
                i2 = v.index(v2.index)
                #print("    i1=%d i2=%d"%(i1,i2))
                if (i1==0 and i2==1) or (i1==1 and i2==2) or (i1==2 and i2==0):
                    v1 = e.verts[0]
                    v2 = e.verts[1]
                else:
                    #print("    reverting verti  ces")
                    v1 = e.verts[1]
                    v2 = e.verts[0]
                lfverts[0]=Vector((v1.co))
                lfverts[1]=Vector((v2.co))
                nonCommonVIndex = list(set([v.index for v in lf.verts]) ^ set([v1.index,v2.index]))[0]
                lfverts[2]=Vector((bm.verts[nonCommonVIndex].co))
                areaOld = lf.calc_area()
                flatFace(lfverts)
                vectorAlF = vertices[vertsDictionary[v1.index]]-vertices[vertsDictionary[v2.index]]
                vectorTF = lfverts[0]-lfverts[1]
                q = vectorTF.rotation_difference(vectorAlF)
                for v in lfverts:
                    v.rotate(q)
                faceTranslate(vertices[vertsDictionary[v1.index]]-lfverts[0],lfverts)
                areaNew = mg.area_tri(lfverts[0], lfverts[1], lfverts[2])
                addFace(vertices, faces, lfverts, lf.index, [v1.index,v2.index,nonCommonVIndex])
                flatted.append(lf.index)
                #print("   flatted %d, old area=%2.2g new area=%2.2g"%(lf.index,areaOld,areaNew))
                refFace.append(lf)
                managedEdge = True

        #print("Finished managing edges")
        if not managedEdge:
            #print("+-- no adjacent to flat, remove from list bfaces")
            if f in bfaces:
                bfaces.remove(f)
        #print("List faces still to manage:%d"%(len(refFace)))
        #for i,face in enumerate(refFace):
            #print("       --+[%d] %d"%(i,face.index))

#print("final verts dict", vertsDictionary)
# print(" faces")
# for f in faces:
#     print(f)
# print(" stored %d vertices "%(len(vertices)))
# for i,v in enumerate(vertices):
#     print("[%d]: x=%2.2g y=%2.2g"%(i,v.x,v.y))

fmesh = bpy.data.meshes.new("Test mesh")
fmesh.from_pydata(vertices, [], faces)
fmesh.update()

bmn = bmesh.new()  # create an empty BMesh
bmn.from_mesh(fmesh)  # fill it in from a Mesh
#bmesh.ops.remove_doubles(bmn, verts=bmn.verts, dist=0.1)
bmn.to_mesh(fmesh)
bmn.free()  # free and prevent further access

fmesh.validate()
fmesh.update()

obj = bpy.data.objects.new("Test", fmesh)
#obj.scale = bpy.context.active_object.scale
bpy.context.scene.objects.link(obj)
