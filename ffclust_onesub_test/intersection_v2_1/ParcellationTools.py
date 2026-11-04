#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Agosto 10 20:09:15 2020

@author: fondecyt-1190701
"""
import numpy as np
import os
import operator
import subprocess as sproc
import networkx as nx
from collections import defaultdict
from time import time
from itertools import combinations
from random import random
from copy import deepcopy
from statistics import mode

from BundleTools import bundleTools as bt
from BundleTools import bundleTools3 as bt3

from Parcellation import visualizationTools as vt

import pickle

def load_restricted_triangles():
    #Existen ciertos triángulos que, por definición, no pueden representar parcelas. Estos corresponden a regiones posterior-inferior de la línea media...
    #...donde ambos hemisferios se unen, ya que en estricto rigor no corresponden a corteza, si no que a materia blanca/gris interna.
    
    #La parcelación obtenida por Lefranc indica cuales triángulos están "restringidos", los cuales están guardados en el archivo .pkl:
    #Lrestricted.pkl para el hemisferio izquierdo.
    #Rrestricted.pkl para el hemisferio derecho.
    
    with open('Parcellation/Lnope.pkl', 'rb') as f:
        Lrestricted = set(pickle.load(f))
    
    with open('Parcellation/Rnope.pkl', 'rb') as f:
        Rrestricted = set(pickle.load(f))
    
    return Lrestricted, Rrestricted

#Obtiene los triángulos por cada subparcela
def triangles_per_subparcel(hemi_dict, sp_dict):
    #Para cada triángulo...
    for tri in hemi_dict.keys():
        for k, v in hemi_dict[tri].items():
            
            #Si es el contador general del triángulo, no hacer nada.
            if k == 'n':
                continue
            
            #Si es el contador de una subparcela, agregarlo al diccionario de subparcelas
            else:
                sp_dict[k][tri] = v

#Elimina subparcelas pequeñas (tamaño medido en cantidad de triángulos)
def filter_small_parcels(hemi_dict, sp_dict, thr):
    #Para cada subparcela...
    for sp in list(sp_dict.keys()):
        #Si tiene menos triángulos que el umbral, eliminar la subparcela del diccionario de subparcelas...
        if len(sp_dict[sp]) < thr:
            small = sp_dict.pop(sp)
            
            #...y también del diccionario de triángulos. Además, reducir el contador "n".
            for tri in list(small.keys()):
                hemi_dict[tri]['n'] -= small[tri]
                del hemi_dict[tri][sp]

#Cálculo del centro de densidad de cada parcela
def calculate_density_center(hemi_dict, dc_dict, thr):
    #Para cada triángulo...
    for tri in hemi_dict.keys():
        for k, v in hemi_dict[tri].items():
            
            #Si es el contador general del triángulo, no hacer nada.
            if k == 'n':
                continue
            
            #De lo contrario, calcular la probabilidad de la subparcela como "contador_subparcela/contador_general".
            else:
                p = v/hemi_dict[tri]['n']
                 
                #Si la probabilidad es mayor que el umbral, guardar en diccionario de centros de densidad.
                if p >= thr:
                    dc_dict[k][tri] = p

#Traslape entre sub-parcelas
def subparcel_overlap(dc, dc_dict, thr):
    #Se crea un grafo no dirigido
    dc_graph = nx.Graph()
    
    #Y una matriz cuadrada que almacena el traslape (idc) para cada par de centros de densidad
    n_dc = len(dc)
    overlaps = np.zeros((n_dc,n_dc))
    
    #Todos contra todos...
    for sp1 in dc_dict.keys():
        for sp2 in dc_dict.keys():
            if sp1 != sp2:
                sp1_tri = set(dc_dict[sp1].keys())
                sp2_tri = set(dc_dict[sp2].keys())
                
                #Calcular la intersección entre ambos centros de densidad
                inter = sp1_tri.intersection(sp2_tri)
                sp_min = np.min((len(sp1_tri),len(sp2_tri)))
                
                #Calcular el idc = interseccion(dc)/min(dc)
                idc = len(inter)/sp_min
                
                #Si es mayor al umbral definido, agregar al grafo. Además, agregar su idc a la matriz.
                if idc >= idc_thr:
                    dc_graph.add_edge(sp1,sp2)
                    overlaps[dc.index(sp1),dc.index(sp2)] = idc
                    overlaps[dc.index(sp2),dc.index(sp1)] = idc
        
    return dc_graph, overlaps

#Fusiona subparcelas con traslape significativo
def fuse_subparcels(dc, dc_graph, overlaps, hemi_dict, sp_dict):
    #Obtiene una lista con los cliques maximales del grafo
    cliques = sorted(nx.find_cliques(dc_graph), key = len, reverse = True)
    
    #Ordena los cliques maximales... primero por largo, y luego por porcentaje de intersección promedio
    cliques_avg_idc = []
    
    #Para cada clique...
    for clique in cliques:
        total_idc = []

        #Se obtienen todas las combinaciones posibles del clique, y para cada par...
        for pair in combinations(clique, 2):
            #...se obtienen el traslape
            total_idc.append(overlaps[dc.index(pair[0]),dc.index(pair[1])])
        
        #Y luego se obtiene el promedio de los traslapes para ese clique.
        cliques_avg_idc.append(np.mean(total_idc))
    
    #Se ordenan los cliques maximales según el traslape promedio
    cliques_temp = sorted(zip(cliques, cliques_avg_idc), key = lambda x: (len(x[0]), x[1]), reverse = True)
    cliques_sorted, avg_idc_sorted = map(list, zip(*cliques_temp))

    #Se fusionan los cliques en orden descendente, descartando aquellos que ya fueron utilizados
    used = set()
    
    for clique in cliques_sorted:
        for sp in clique:
            if sp in used:
                continue
            else:
                ref = sp
                used.add(sp)
                break
            
        for sp in clique:
            if sp == ref or sp in used:
                continue
            
            else:
                for tri in list(sp_dict[sp].keys()):
                    hemi_dict[tri][ref] += sp_dict[sp][tri]
                    del hemi_dict[tri][sp]            
                used.add(sp)
                
#Obtiene parcelación dura (determinista, asignando a cada triángulo la subparcela con mayor probabilidad)
def hard_parcellation(hemi_dict, sp_hard):
    #Para cada triángulo...
    for tri in list(hemi_dict.keys()):
        #Eliminar el contador general. Además, si el triángulo queda sin parcelas, eliminar su entrada.
        del hemi_dict[tri]['n']
        if len(hemi_dict[tri]) == 0:
            del hemi_dict[tri]
    
    #Para cada triángulo, obtener la subparcela con mayor probabilidad,
    for tri in list(hemi_dict.keys()):        
        m = max(hemi_dict[tri].items(), key=operator.itemgetter(1))[0]    
        hemi_dict[tri] = m
        
    #Agregar la parcela más probable al diccionario de parcelación dura, junto con sus triángulos.
    for tri, sp in hemi_dict.items():
        sp_hard[sp].append(tri)
        
#Componentes Conectadas Principales
def PCC(Tri, polygons): #Autor: Felipe Silva
    if len(Tri) == 0: return Tri;
    edgesPolys = [];

    for ind in Tri:
        edgesPolys.append([polygons[ind,0], polygons[ind,1]]);
        edgesPolys.append([polygons[ind,0], polygons[ind,2]]);
        edgesPolys.append([polygons[ind,1], polygons[ind,2]]);
    
    edgesPolys = np.unique(edgesPolys, axis=0);
    
    G = nx.Graph();
    G.add_edges_from(edgesPolys);
    
    cc = list(nx.connected_components(G));
    len_cc = [len(comp) for comp in cc];
    regions_cc = list(cc[np.argmax(len_cc)]);
    
    ind = [np.where(polygons[Tri]==region_cc)[0] for region_cc in regions_cc];
    ind = np.unique(np.concatenate(ind));
    return Tri[ind]

#Proyecta la parcelación obtenida sobre un mallado cortical ARCHI
def visualize_parcellation(meshes_path, L_sp, R_sp, sub, seed = False):
    #Cargar trianguos restringidos
    Lrestricted, Rrestricted= load_restricted_triangles()
    
    #Semilla para colores aleatorios
    if seed != False:
        np.random.seed(seed)
    
    #Parcelas finales a graficar
    final_parcels = set()

    for k in L_sp.keys():
        final_parcels.add(k)
        
    for k in R_sp.keys():
        final_parcels.add(k)

    fp = list(final_parcels)    

    #Paleta de colores según cantidad de parcelas
    paleta = [(np.random.random(), np.random.random(), np.random.random()) for i in range(len(fp))]

    #Directorios de los mallados corticales
    Lhemi_path = meshes_path + sub + '/lh.obj'; # left hemisphere path
    Rhemi_path = meshes_path + sub + '/rh.obj'; # right hemisphere path

    #Lectura de mallados
    Lvertex, Lpolygons = bt.read_mesh_obj(Lhemi_path)
    Rvertex, Rpolygons = bt.read_mesh_obj(Rhemi_path)
    
    Lhemi = vt.Polygon(Lvertex, Lpolygons);
    Rhemi = vt.Polygon(Rvertex, Rpolygons);
    
    Lhemi.setOpacity(1);
    Rhemi.setOpacity(1);
    
    #Creación del render a visualizar
    render = vt.Render();
    
    #Se renderizan los mallados
    render.AddActor(Lhemi);
    render.AddActor(Rhemi);
    
    #Para cada parcela del hemisferio izquierdo...
    for k, v in L_sp.items():
        #Se selecciona un color de la paleta
        color = paleta[fp.index(k)]
        
        if len(v) == 0:
            continue
      
        #De todos los triángulos con parcelas, se eliminan aquellos que pertenecen al conjunto de triángulos restringidos.
        v_restricted = list(set(v).difference(Lrestricted))

        #Se renderizan los triángulos y polígonos.
        sp_tri = vt.Polygon(Lvertex, Lpolygons[v_restricted]);
        sp_tri.setColor((color[0], color[1], color[2]));
        render.AddActor(sp_tri);

    #Ídem para el derecho
    for k, v in R_sp.items():
        color = paleta[fp.index(k)]
    
        if len(v) == 0:
            continue
    
        v_restricted = list(set(v).difference(Rrestricted))
    
        sp_tri = vt.Polygon(Rvertex, Rpolygons[v_restricted]);
        sp_tri.setColor((color[0], color[1], color[2]));
        render.AddActor(sp_tri);
    
    render.Start();
    del render

# Determina las subparcelas preliminares por cada triángulo
def preliminar_subparcels (start, end, sub_selected, intersection_path, meshes_path):
    print('Parcelas preliminares')
    #Lectura de polígonos del mallado izquierdo y derecho.
    Lhemi_path = meshes_path + sub_selected + '/lh.obj'; 
    Rhemi_path = meshes_path + sub_selected + '/rh.obj'; 
    
    Lvertex, Lpolygons = bt.read_mesh_obj(Lhemi_path)
    Rvertex, Rpolygons = bt.read_mesh_obj(Rhemi_path)
    
    #Cálculo de vecindad para cada triángulo del mallado (triángulos que comparten una arista o un vértice).
    Lneighbors = bt.mesh_neighbors(Lpolygons)
    Rneighbors = bt.mesh_neighbors(Rpolygons)
    
    # Creación de diccionarios
    
    # Key: Value => Subparcela: [triángulos en los que aparece la subparcela]
    L_sp_dict = defaultdict(lambda: defaultdict(float))
    R_sp_dict = defaultdict(lambda: defaultdict(float))
    
    # Key: Value => Triángulo: [subparcelas que llegan a ese triángulo, y cuántas veces lo hace cada una]
    Lhemi_dict = defaultdict(lambda: defaultdict(int))
    Rhemi_dict = defaultdict(lambda: defaultdict(int))
    
    # (i.e. cuál subparcela aparece en cada triángulo, y cuántas veces)
    for i in range(start, end+1):
        
        t0 = time()
        
        #Nombre del sujeto
        if i < 10:
            sub = '00' + str(i)
            
        else:
            sub = '0' + str(i)
        
        print(str(i)+'/'+str(end))
        
        #Tipo de intersección
        int_type = ['L-L','R-R','Inter']
        
        #Se itera sobre todas las intersecciones
        for tipo in int_type:
            int_path = intersection_path + sub + '/Final/' + tipo + '/'
            int_list = os.listdir(int_path)
            
            #Conteo de subparcelas preliminares por triángulo.
                #Por definición, cada clúster define dos subparcelas:
                    #En su extremo inicial, la subparcela clúster_A
                    #En su extremo final, la subparcela clúster_B
                    
            for file in int_list:
                cluster_name = file.split('.')[0].split('_')[1]
                
                #Dependiendo del tipo de intersección (intra o inter), se define el nombre la subparcela:
                    #Para Intra-Left: Tanto el extremo A como el extremo B están en el diccionario de subparcelas izquierdas.
                    #Para Intra-Right: Tanto el extremo A como el extremo B están en el diccionario de subparcelas derechas.
                    #Para Ìnter: El extremo A està en el lado izquierdo, y el B en el derecho.
                    
                if tipo == 'L-L':
                    if cluster_name + 'A' not in L_sp_dict:
                        L_sp_dict[cluster_name + 'A'] = {}
                
                    if cluster_name + 'B' not in L_sp_dict:
                        L_sp_dict[cluster_name + 'B'] = {}
                        
                elif tipo == 'Inter':
                    if cluster_name + 'A' not in L_sp_dict:
                        L_sp_dict[cluster_name + 'A'] = {}
                
                    if cluster_name + 'B' not in R_sp_dict:
                        R_sp_dict[cluster_name + 'B'] = {}
                        
                elif tipo == 'R-R':
                    if cluster_name + 'A' not in R_sp_dict:
                        R_sp_dict[cluster_name + 'A'] = {}
                
                    if cluster_name + 'B' not in R_sp_dict:
                        R_sp_dict[cluster_name + 'B'] = {}
                
                #Lectura de intersección y obtención de triángulos iniciales y finales
                ix = bt.read_intersection(int_path + file)
                in_tri = ix[0]
                fn_tri = ix[1]
                
                #Por cada triángulo, aumentar contador para la subparcela correspondiente.
                #Además, aumentar el contador general "n" del triángulo.
                for tri in in_tri:            
                    if tipo == 'L-L' or tipo == 'Inter':
                        Lhemi_dict[tri]['n'] += 1 
                        Lhemi_dict[tri][cluster_name + 'A'] += 1
                        
                        nbors = Lneighbors[tri]
                        
                        for nbor in nbors:
                            if nbor in in_tri:
                                Lhemi_dict[tri]['n'] += 1 
                                Lhemi_dict[tri][cluster_name + 'A'] += 1
                        
                    elif tipo == 'R-R':
                        Rhemi_dict[tri]['n'] += 1 
                        Rhemi_dict[tri][cluster_name + 'A'] += 1
                        
                        nbors = Rneighbors[tri]
                        
                        for nbor in nbors:
                            if nbor in in_tri:
                                Rhemi_dict[tri]['n'] += 1
                                Rhemi_dict[tri][cluster_name + 'A'] += 1
                    
                for tri in fn_tri:
                    if tipo == 'L-L':
                        Lhemi_dict[tri]['n'] += 1 
                        Lhemi_dict[tri][cluster_name + 'B'] += 1
                        
                        nbors = Lneighbors[tri]
                        
                        for nbor in nbors:
                            if nbor in in_tri:
                                Lhemi_dict[tri]['n'] += 1 
                                Lhemi_dict[tri][cluster_name + 'B'] += 1
                        
                    elif tipo == 'Inter' or tipo == 'R-R':
                        Rhemi_dict[tri]['n'] += 1 
                        Rhemi_dict[tri][cluster_name + 'B'] += 1
                        
                        nbors = Rneighbors[tri]
                        
                        for nbor in nbors:
                            if nbor in in_tri:
                                Rhemi_dict[tri]['n'] += 1
                                Rhemi_dict[tri][cluster_name + 'B'] += 1
                                
    
        print(str(time() - t0) + '[s]')    
                
    #Conteo de triángulos para cada subparcela (i.e. en cuáles triángulos aparece X subparcela, y cuántas veces llega ahí)
            
    triangles_per_subparcel(Lhemi_dict, L_sp_dict)            
    triangles_per_subparcel(Rhemi_dict, R_sp_dict)
   
    #Filtrado de subparcelas pequeñas (por cantidad de triángulos insuficientes). Incluye actualización de conteo en primer diccionario.

    size_thr = len(Lpolygons)/1000.0
                    
    filter_small_parcels(Lhemi_dict, L_sp_dict, size_thr)
    filter_small_parcels(Rhemi_dict, R_sp_dict, size_thr)
    
    return Lhemi_dict, Rhemi_dict, L_sp_dict, R_sp_dict

#Carga parcelas ('hard','cc','final')    
def load_parcels (parcels, parcels_path): #hard cc final
    Lparcels_path =  parcels_path + '/'+ parcels +'_parcels/left/'
    Rparcels_path = parcels_path + '/'+ parcels +'_parcels/right/'

    Lparcels = dict()
    Rparcels = dict()

    for Dir in os.listdir(Lparcels_path):
        Lparcels[Dir.split('.')[0]] = bt.read_parcels(Lparcels_path + Dir)
        
    for Dir in os.listdir(Rparcels_path):
        Rparcels[Dir.split('.')[0]] = bt.read_parcels(Rparcels_path + Dir)
            
    return Lparcels, Rparcels

#Obtencion de parcelas duras
def hard_parcels (start, end, sub_selected, intersection_path, meshes_path, dc_thr, idc_thr, output):
    #Parcelas preliminares
    Lhemi_dict, Rhemi_dict, L_sp_dict, R_sp_dict = preliminar_subparcels (start, end, sub_selected, intersection_path, meshes_path)
    
    print('Parcelas duras')
    # Key: Value => Subparcela: [triángulos correspondientes al centro de densidad de la subparcela]
    L_dc_dict = defaultdict(lambda: defaultdict(float))
    R_dc_dict = defaultdict(lambda: defaultdict(float))
    
    print('Calculo de centro de densidad')
    #Cálculo de probabilidades y centros de densidad de cada subparcela.
    calculate_density_center(Lhemi_dict, L_dc_dict, dc_thr)
    calculate_density_center(Rhemi_dict, R_dc_dict, dc_thr)
    print('OK')
    
    print('Modelacion de traslape')
    # Modelación del traslape entre subparcelas.
    L_dc = list(L_dc_dict.keys())
    R_dc = list(R_dc_dict.keys())
    
    L_dc_graph, L_overlaps = subparcel_overlap(L_dc, L_dc_dict, idc_thr)
    R_dc_graph, R_overlaps = subparcel_overlap(R_dc, R_dc_dict, idc_thr)
    print('OK')
    
    print('Fusion de parcelas preliminares')
    # Fusión de subparcelas cuyo centro de densidad tenga intersección significativa.
    fuse_subparcels(L_dc, L_dc_graph, L_overlaps, Lhemi_dict, L_sp_dict)
    fuse_subparcels(R_dc, R_dc_graph, R_overlaps, Rhemi_dict, R_sp_dict)
    print('OK')
    
    # Obtención de parcelación dura (según probabilidades).
    #Se utilizan copias para mantener las probabilidades originales.
    Lhemi_hard = deepcopy(Lhemi_dict)
    Rhemi_hard = deepcopy(Rhemi_dict)
    
    Lparcels_hard = defaultdict(list)
    Rparcels_hard = defaultdict(list)
    print('Obtencion de parcelas duras')
    hard_parcellation(Lhemi_hard, Lparcels_hard)
    hard_parcellation(Rhemi_hard, Rparcels_hard)
    print('OK')
    print('Guardando Parcelas duras')
    # Creación de directorios para almacenar hard parcels
    Lhard_path =  output + '/hard_parcels/left/'
    Rhard_path = output + '/hard_parcels/right/'
    
    if not os.path.exists(Lhard_path):
        os.makedirs(Lhard_path)
    
    if not os.path.exists(Rhard_path):
        os.makedirs(Rhard_path)
    
    Lparcels_hard_c = deepcopy(Lparcels_hard)
    Rparcels_hard_c=deepcopy(Rparcels_hard)
    
    # Escritura de hard parcels
    for k, v in Lparcels_hard_c.items(): #revisar
        if 1000000 in v:
            del Lparcels_hard[k]
            print('Parcel '+k+ ' Eliminada') 
            
    for k, v in Rparcels_hard_c.items():
        if 1000000 in v:
            del Rparcels_hard[k]
            print('Parcel '+k+ ' Eliminada')
            
    for sp, tris in Lparcels_hard.items():
        bt.write_parcels(Lhard_path + sp + '.parcelsdata', tris)
    
    for sp, tris in Rparcels_hard.items():
        bt.write_parcels(Rhard_path + sp + '.parcelsdata', tris)
    print('Listo')
    return Lparcels_hard, Rparcels_hard

#Obtencioon de parcelas finales
def final_parcels (output, meshes_path, sub_selected, ero, dil):
    print('Parcelas finales')
    #Creacion de directorios par almacenar parcelas finales
    #Parcelación CC
    Lparcelcc_path = output + '/cc_parcels/left/'
    Rparcelcc_path = output + '/cc_parcels/right/'
    
    if not os.path.exists(Lparcelcc_path):
        os.makedirs(Lparcelcc_path)
    
    if not os.path.exists(Rparcelcc_path):
        os.makedirs(Rparcelcc_path)    
    
    #Parcelación final (post-procesada)
    Lfinal_path = output + '/final_parcels/left/'
    Rfinal_path = output + '/final_parcels/right/'
    
    if not os.path.exists(Lfinal_path):
        os.makedirs(Lfinal_path)
    
    if not os.path.exists(Rfinal_path):
        os.makedirs(Rfinal_path)
        
    #Lectura de hard parcels
    Lparcels_hard, Rparcels_hard = load_parcels('hard', output)
    
    print('Calculo de componente conexa mas grande')
    #Cálculo de componente conexa más grande de cada parcela
    Lparcel_cc = {}
    Rparcel_cc = {}
    Lhemi_path = meshes_path + sub_selected + '/lh.obj'; 
    Rhemi_path = meshes_path + sub_selected + '/rh.obj'; 
    
    Lvertex, Lpolygons = bt.read_mesh_obj(Lhemi_path)
    Rvertex, Rpolygons = bt.read_mesh_obj(Rhemi_path)
    
    for k,v in Lparcels_hard.items():
        Lparcel_cc[k] = PCC(np.array(list(v)), Lpolygons)
        #PCC(np.array(v), Lpolygons)
        
    for k,v in Rparcels_hard.items():
        Rparcel_cc[k] = PCC(np.array(list(v)), Rpolygons)
        #PCC(np.array(v), Rpolygons)
    print('OK')
    
    #Escritura de parcelas CC
    for sp, tris in Lparcel_cc.items():
        bt.write_parcels(Lparcelcc_path + sp + '.parcelsdata', tris)
    
    for sp, tris in Rparcel_cc.items():
        bt.write_parcels(Rparcelcc_path + sp + '.parcelsdata', tris)
    
    #Lectura de parcelas CC
    Lparcel_cc, Rparcel_cc = load_parcels('cc', output)
    print('Opening morfologico y guardado de parcelas finales')
    #Opening morfológico (erosión + dilatación)    
    sproc.call(['make']);
    sproc.call(['Parcellation/./felipe', Lhemi_path, Rhemi_path, str(ero), str(dil),  Lparcelcc_path, Rparcelcc_path, Lfinal_path, Rfinal_path]);

    #Lectura de parcelas finales
    Lparcels_final, Rparcels_final = load_parcels('final', output)    
    print('Listo')
    return Lparcel_cc, Rparcel_cc, Rparcels_final, Lparcels_final

#Formula del coeficiente DICE
def DSC(A,B):
    set_A = set(A)
    set_B = set(B)
    
    interx = set_A.intersection(set_B)
    
    return (2*len(interx))/((len(set_A)+len(set_B)))

#Calculo de coeficiente DICE para comparar parcelas con atlas
def dise_comparision (atlas_comparision, parcels_path, dice_thr):
    # Lectura de parcelas propias.    
    lh_mine_dict, rh_mine_dict = load_parcels('final', parcels_path)
        
    # Cargar parcelación del atlas seleccionado.
    with open(atlas_comparision + '_Lparcels.pkl', 'rb') as f:   
        lh_atlas_dict = pickle.load(f)
    
    with open(atlas_comparision + '_Rparcels.pkl', 'rb') as f:   
        rh_atlas_dict = pickle.load(f)    
        
    # Cálculo de DSC entre parcelas.
    # Diccionario parcela_propia:parcela_atlas, para aquellos casos que superan el umbral de DSC.
    lh_dice_dict = defaultdict(int)
    rh_dice_dict = defaultdict(int)

    # Lista que almacena las parcelas ya utilizadas.
    used = []

    # Iteración sobre todas las parcelas propias obtenidas.
    for parcel, tris in lh_mine_dict.items():
        #Se registra el mayor DSC obtenido entre todas las comparaciones, y el nombre de la parcela asociada.
        dice_max = 0 
        winner = ''
      
        #Se itera sobre las parcelas del atlas.
        for ref, rtris in lh_atlas_dict.items():
            
            #Si la parcela ya fue utilizada, ignorar.
            if ref in used:
                continue
            
            #Cálculo de DSC entre parcela propia y parcela del atlas.
            dice = DSC(tris,rtris)
            
            #Si es mayor al máximo obtenido, actualizar.
            if dice > dice_max:
                dice_max = dice
                winner = ref
        
        #Si el mayor valor obtenido supera el umbral de similitud, se agrega al diccionario de parcelas similares, y a la lista de parcelas ya utilizadas.
        if dice_max >= dice_thr:
            lh_dice_dict[parcel] = winner
            used.append(winner)
            
        # Si no supera el umbral, continuar con la siguiente parcela.
        else:
            continue
    
    #Se repite lo mismo para el hemisferio derecho.
    used = []
    
    for parcel, tris in rh_mine_dict.items():
        dice_max = 0
        winner = ''
      
        for ref, rtris in rh_atlas_dict.items():
            if ref in used:
                continue
            dice = DSC(tris,rtris)
            if dice > dice_max:
                dice_max = dice
                winner = ref
            
        if dice_max >= dice_thr:
            rh_dice_dict[parcel] = winner
            used.append(winner)
            
        else:
            continue  
    
    #Diccionario de parcelas similares, considerando los triángulos asociados en el caso de la parcelación propia.
    lh_mine_common = {k:lh_mine_dict[k] for k in list(lh_dice_dict.keys()) if k in lh_mine_dict}
    rh_mine_common = {k:rh_mine_dict[k] for k in list(rh_dice_dict.keys()) if k in rh_mine_dict}
    
    #Diccionario de parcelas similares, considerando los triángulos asociados en el caso de la parcelación del atlas.
    lh_atlas_common = {k:lh_atlas_dict[k] for k in list(lh_dice_dict.values()) if k in lh_atlas_dict}
    rh_atlas_common = {k:rh_atlas_dict[k] for k in list(rh_dice_dict.values()) if k in rh_atlas_dict}
    
    if not os.path.exists(parcels_path+'Dise/'):
        os.makedirs(parcels_path+'Dise/')
        
    save_txt(parcels_path+'Dise/lh_list_'+atlas_comparision.replace("/",""),lh_mine_common, dice_thr)
    save_txt(parcels_path+'Dise/rh_list_'+atlas_comparision.replace("/",""), rh_mine_common, dice_thr)
    
    return lh_mine_common, rh_mine_common, lh_atlas_common, rh_atlas_common

def save_txt (name, diccionary, dice_thr):
    count=0
    with open(name+'_'+str(int(dice_thr*100))+'.txt', 'w') as f:
        for line in diccionary.keys():
            f.write(line)
            f.write('\n')
            count+=1
        f.write('\nTotal=%d' % count)
   #%%     
t0 = time()
#Sujeto base. Puede ser cualquiera, ya que todos los mallados tienen triángulos correspondientes.
sub = '001'
theta_QB=20
intersection_path= 'Archi_y_codigos/ARCHI/subs/filtered_centroids6/intersection/intersection_QB_'+str(theta_QB)+'/'
meshes_path= 'Archi_y_codigos/ARCHI/subs/meshes_obj/' 
#0.1idc - 0,2dc
idc_thr = 0.3
dc_thr = 0.4
dice_thr=0.5
ero = 1
dil = 6
output_parcellation= 'Parcellation/QB_'+str(theta_QB)+'/IDC_DC_'+str(idc_thr*100)+'%_'+str(dc_thr*100)+'%/' #revisar
# Selección de atlas a comparar.
atlases = ['atlas/Lefranc','atlas/Brainnetome','atlas/Narciso']
semilla_visualizacion = 48

#Parcelas preliminares
#preliminar_subparcels (1, 79, sub, intersection_path, meshes_path)

#Parcelas duras
#Lparcels_hard, Rparcels_hard= hard_parcels (1, 79, sub, intersection_path, meshes_path, dc_thr, idc_thr, output_parcellation)

#Parcelas finales
#Lparcel_cc, Rparcel_cc, Rparcels_final, Lparcels_final=  final_parcels (output_parcellation, meshes_path, sub, ero, dil)

#Comparacion Dice
#lh_mine_common, rh_mine_common, lh_atlas_common, rh_atlas_common= dise_comparision(atlases[0], output_parcellation, dice_thr)

#Visualizacion de parcelas
#visualize_parcellation(meshes_path, Lparcels_hard, Rparcels_hard, '001', seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcel_cc, Rparcel_cc, '001', seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_final, Rparcels_final, '001', seed = semilla_visualizacion)

#Visualizacion para DISE
#
#visualize_parcellation(meshes_path, lh_atlas_common, rh_atlas_common, '001', seed = 47)
#visualize_parcellation(meshes_path, lh_mine_common, rh_mine_common, '001', seed = 47)

#Lparcels_final, Rparcels_final= load_parcels('final', output_parcellation)
#Lparcels_hard, Rparcels_hard= load_parcels('hard', output_parcellation)
# Lparcels_cc, Rparcels_cc= load_parcels('cc', output_parcellation)
#visualize_parcellation(meshes_path, Lparcels_final, Rparcels_final, sub, seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_hard, Rparcels_hard, '001', seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_cc, Rparcels_cc, '001', seed = semilla_visualizacion)
