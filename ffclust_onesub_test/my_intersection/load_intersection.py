from BundleTools import bundleTools as BT
from BundleTools import bundleTools3 as BT3
import numpy as np
import os
import sys
import math
import statistics
import matplotlib.pyplot as plt
import visualizationTools as vt
import visual_tools as vs
from dipy.segment.metric import EuclideanMetric
from dipy.segment.metric import mdf
from dipy.segment.metric import dist
from subprocess import call
from os import path
import pandas as pd
from collections import defaultdict
import vtk
from dipy.segment.metric import CosineMetric


def read_intersection( infile ):

    f = open(infile, 'rb');

    total_triangles = np.frombuffer( f.read( 4 ), np.uint32 )[ 0 ];

    InTri = np.frombuffer( f.read( 4 * total_triangles), np.uint32 );
    FnTri = np.frombuffer( f.read( 4 * total_triangles), np.uint32 );

    InPoints = np.frombuffer( f.read( 4 * 3 * total_triangles), np.float32 ).reshape(-1,3);
    FnPoints = np.frombuffer( f.read( 4 * 3 * total_triangles), np.float32 ).reshape(-1,3);

    fib_index = np.frombuffer( f.read( 4 * total_triangles), np.uint32 );

    f.close();
    return InTri, FnTri, InPoints, FnPoints, fib_index;

def read_intersection2( infile ):

    f = open(infile, 'rb');

    total_triangles = np.frombuffer( f.read( 4 ), np.uint32 )[ 0 ];

    InTri = np.frombuffer( f.read( 4 * total_triangles), np.uint32 );
    FnTri = np.frombuffer( f.read( 4 * total_triangles), np.uint32 );

    InPoints = np.frombuffer( f.read( 4 * 3 * total_triangles), np.float32 ).reshape(-1,3);
    FnPoints = np.frombuffer( f.read( 4 * 3 * total_triangles), np.float32 ).reshape(-1,3);

    fib_index = np.frombuffer( f.read( 4 * total_triangles), np.uint32 );
    fiber_class = np.frombuffer( f.read( 4 * total_triangles), np.uint32);

    f.close();
    return InTri, FnTri, InPoints, FnPoints, fib_index, fiber_class;

def return_fiberclass(fib_class):
    fibs_lh = list(np.where(np.array(fib_class) == 1)[0])
    fibs_rh = list(np.where(np.array(fib_class) == 2)[0])
    fibs_inter_l_r = list(np.where(np.array(fib_class) == 3)[0])
    fibs_inter_r_l = list(np.where(np.array(fib_class) == 4)[0]) 
    
    return fibs_lh, fibs_rh, fibs_inter_l_r, fibs_inter_r_l

archivo_intersection = 'results/tractography-streamline-regularized-deterministic_001_21p.intersectiondata'


InTri, FnTri, InPoints, FnPoints, fib_index, fiber_class = read_intersection2(archivo_intersection)
fibs_lh, fibs_rh, fibs_inter_l_r, fibs_inter_r_l = return_fiberclass(fiber_class)

fibers = BT.read_bundle('set_tractografia/tractography-streamline-regularized-deterministic_001_21p.bundles')

#LEER MESHES.
#Lectura de datos y definición de objetos para VTK
render = vt.Render();
Rhemi_path = "meshes/rh.obj"
Rvertex, Rpolygons = BT.read_mesh_obj(Rhemi_path)
Rhemi = vt.Polygon(Rvertex, Rpolygons);
Rhemi.setOpacity(0.20)
render.AddActor(Rhemi);

#Lectura de datos y definición de objetos para VTK
Lhemi_path = "meshes/lh.obj"
Lvertex, Lpolygons = BT.read_mesh_obj(Lhemi_path)
Lhemi = vt.Polygon(Lvertex, Lpolygons);
Lhemi.setOpacity(0.20)
render.AddActor(Lhemi);


intris_visualize_lh = []
intris_visualize_rh = []
fntris_visualize_lh = []
fntris_visualize_rh = []

n = 122

for i in fibs_lh:#[n:n+10]:
    if False:
        print("fibs_lh: " + str(fibs_lh.index(i)) + '/' + str(len(fibs_lh)))
        fib = vt.Line(fibers[fib_index[i]])
        fib.setColor((1,0,0))
        #fib.setJetColormap()
        render.AddActor(fib)

        #append tris
        intris_visualize_lh.append(InTri[i])
        fntris_visualize_lh.append(FnTri[i])

for i in fibs_rh:#[n:n+10]:
    if False:
        print("fibs_rh: " + str(fibs_rh.index(i)) + '/' + str(len(fibs_rh)))
        print(fib_index[i])
        #print(fibers[fib_index[i]])
        fib = vt.Line(fibers[fib_index[i]])
        #fib = vt.Line(fibers[i])
        #fib.setJetColormap()
        fib.setColor((0,0,1))
        render.AddActor(fib)
        
        #render_points(fibers[fib_index[i]], i, InPoints, FnPoints, render)
        
        #append tris
        intris_visualize_rh.append(InTri[i])
        fntris_visualize_rh.append(FnTri[i])  

for i in fibs_inter_l_r:#[n:n+10]:
    if False:
        print("fibs_inter L to R: " + str(fibs_inter_l_r.index(i)) + '/' + str(len(fibs_inter_l_r)))
        fib = vt.Line(fibers[fib_index[i]])
        fib.setColor((1,0,1))
        #fib.setJetColormap()
        render.AddActor(fib)
       
    #append tris
    intris_visualize_lh.append(InTri[i])
    fntris_visualize_rh.append(FnTri[i])    
                
for i in fibs_inter_r_l:#[n:n+10]:
    if False:
        print("fibs_inter R to L: " + str(fibs_inter_r_l.index(i)) + '/' + str(len(fibs_inter_r_l)))
        fib = vt.Line(fibers[fib_index[i]])
        fib.setColor((1,0,1))
        #fib.setJetColormap()
        render.AddActor(fib)
   
    #append tris
    fntris_visualize_lh.append(FnTri[i])
    intris_visualize_rh.append(InTri[i])  

triangles_rh_dict = {'in': intris_visualize_rh, 'fn': fntris_visualize_rh}
triangles_lh_dict = {'in': intris_visualize_lh, 'fn': fntris_visualize_lh}

#Triangles plot

for k,v in triangles_lh_dict.items():
    print("tris lh")
    if k == 'in':
        tri = vt.Polygon(Lvertex, Lpolygons[v[0:len(v)]])
        tri.setColor((0,0,1))
        render.AddActor(tri)
        
    elif k == 'fn':
        tri = vt.Polygon(Lvertex, Lpolygons[v[0:len(v)]])
        tri.setColor((1,0,0))
        render.AddActor(tri)

for k,v in triangles_rh_dict.items():
    print("tris rh")
    if k == 'in':
        tri = vt.Polygon(Rvertex, Rpolygons[v[0:len(v)]])
        tri.setColor((0,0,1))
        render.AddActor(tri)
        
    elif k == 'fn':
        tri = vt.Polygon(Rvertex, Rpolygons[v[0:len(v)]])
        tri.setColor((1,0,0))
        render.AddActor(tri)
    
render.Start();
del render
