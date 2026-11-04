#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 20:51:47 2021

@author: fondecyt-1190701
"""

import os
import numpy as np
import nibabel as nib
import random
from dipy.io.streamline import load_tractogram
from dipy.tracking.streamline import Streamlines
from dipy.segment.clustering import QuickBundles
from dipy.segment.metric import IdentityFeature
from dipy.segment.metric import AveragePointwiseEuclideanMetric
from dipy.segment.metric import mdf
from dipy.segment.metric import dist
from dipy.segment.metric import EuclideanMetric
from dipy.io.pickles import save_pickle
from dipy.io.pickles import load_pickle
from dipy.data import get_fnames
from dipy.viz import window, actor
import time
from BundleTools import bundleTools as bt
from BundleTools import bundleTools3 as bt3
import shutil 
from numpy.linalg import inv
import matplotlib.pyplot as plt

byte_order = "DCBA"

#Formula para alineamiento de clusters
def align_cluster(cluster, ref):
    
    aligned_cluster = []
    
    for fiber in cluster:
        d = dist(EuclideanMetric(), fiber, ref)/21
        f = mdf(fiber, ref)
        
        if f < d:
            aligned_cluster.append(fiber[::-1])
            
        else:
            aligned_cluster.append(fiber)
            
    return aligned_cluster

#Calculo de centroide del cluster ínter-sujeto a nivel intra-sujeto 
def calculate_centroid(cluster):
    return np.mean(cluster, axis = 0)

#Filtrado de centroides get_centroids(data, sub, 15, 1.5, 100)
def get_centroids(brain, sub, remove_little_ones=False, size=False, nearest=False):
    # remove_little_ones: Cantidad mínima de fibras por cluster (para eliminar clusters poco densos). 
    # size: Distancia mínima entre los dos primeros puntos del centroide (para eliminar fibras muy cortas).
    # nearest: Sólo se selecciona un nearest% (porcentaje) aleatorio de las fibras del cluster (para disminuir la carga computacional).
    
    centroids = []
    original_data = []
    
    og_centroids = np.array(bt3.read_bundle_severalbundles(path + sub + '/merged/ffclust4/centroids.bundles')[0])
    
    if(remove_little_ones):
        removed1 = 0
        removed2 = 0
        for i,cluster in enumerate(brain):
            sz = len(cluster)
            if(sz >= remove_little_ones):
                centroid = og_centroids[i][0]
                if np.linalg.norm(centroid[0]-centroid[1]) > size:
                    centroids.append(centroid)
                    selected = random.sample(range(sz),int(sz*(nearest/100)))
                    original_data.append(np.array(cluster)[selected])
                else:
                    removed2 +=1
            else:
                removed1+=1
        print('Centroide '+sub+' filtrado')
        #print(str(removed1)+' and '+str(removed2)+' removed from '+str(len(brain))+'. '+str(len(centroids))+' remaining centroids.'+str(len(np.concatenate(brain)))+' remaining fibers')
    else:
        for cluster in brain:
            centroids.append(get_centroid(cluster))
            original_data.append(np.array(cluster))
    return centroids, np.array(original_data)

#Creacion de paleta de colores
def random_palette(n):
    return [(random.randrange(0, 255), random.randrange(0, 255), random.randrange(0, 255)) for i in range(n)]

#Escritura de archivo jerárquico para visualización por clusters en Anatomist
def write_hie(path,names,colors):
    f = open(path,"w")
    f.write("# tree 1.0\n\n*BEGIN TREE hierarchy\ngraph_syntax RoiArg\n\n*BEGIN TREE fold_name\nname ALL\n\n")
    for i,n in enumerate(names):
        f.write("*BEGIN TREE fold_name\n")
        f.write("name "+n+"\n")
        f.write("color "+str(colors[i][0])+" "+str(colors[i][1])+" "+str(colors[i][2])+"\n\n")
        f.write("*END\n\n")
    f.write("*END\n\n*END\n\n*END\n\n")
    f.close()

#Guardado de centroides y fibras filtradas
def save_filtered_centroids (start, end, path):
    all_centroids = []    #Arreglo que contiene todos los centroides
    all_data = []         #Arreglo que contiene todas las fibras
    ids = []              #Arreglo que contiene todos los identificadores (a qué sujeto pertenece cada cluster)
    t0 = time.time()
    
    if not os.path.exists(path+ '/filtered_centroids6/'):
                os.makedirs(path+'/filtered_centroids6/')
                
    for x in range(start, end+1): #Iteración sobre todos los sujetos
        
        #Se asume que los sujetos están en carpetas desde "001" hasta "079"
        if x < 10:
            sub = '00' + str(x)
        else:
            sub = '0' + str(x)
        
        #Lectura de clusters
        data = np.array(bt3.read_bundle_severalbundles(path + sub + '/merged/ffclust4/finalClusters.bundles')[0])
        
        #Obtención de fibras y centroides, y guardado en arreglos
        centroids, data = get_centroids(data, sub, 15, 1.5, 100)   #En definición de get_centroids se explican los parámetors
                        #Cantidad min de fibras por cluster, Largo centroide, Carga Computacional
        all_centroids.append(centroids)
        all_data.append(data)
        
        #Cálculo de identificadores
        ids += [x]*len(data) # se guarda con un identificador (de 1 a 79) tantas veces como largo tenga el cluster filtrado 
    print('Tiempo Total filtrado de centroides: ' + str(time.time() - t0) + ' [s]')    
    print('Saving all_data...')
    
    #Guardado de fibras, centroides, e ids en archivos .npy (optimizado para arreglos numpy)
    np.save(path+'filtered_centroids6/all_centroids', np.concatenate(all_centroids))
    np.save(path+'filtered_centroids6/ids', np.array(ids))
    np.save(path+'filtered_centroids6/all_data', np.concatenate(all_data))
    print('Tiempo Total guardado de centroides: ' + str(time.time() - t0) + ' [s]')
    
    os._exit(00) #Restart Kernel for freeing memory
    
    print('Saved...')
   

#Carga de centroides filtrados
def load_filtered_centroids (centroids_path):
    #Carga de identificadores => elemento i indica a qué sujeto pertenece el centroide i.
    ids = np.load(centroids_path + 'ids.npy')
    
    #Carga de centroides intra-sujeto y conversión a ArraySequence (AS, formato utilizado por DIPY/QuickBundles)
    intra_centroids = np.load(centroids_path + 'all_centroids.npy')
    intra_centroids_AS = nib.streamlines.ArraySequence(intra_centroids)
    
    return ids, intra_centroids_AS

#Carga de clusters filtrados
def load_filtered_clusters (clusters_path):
    #Arreglo de fibras...
    fibers = np.empty((0))
    fibers= np.load(clusters_path + 'all_data.npy', allow_pickle = True)
    
    return fibers

#Ejecucion de algoritmo de interclustering QuickBundles
def interclustering (theta_threshold, centroids_data, output_data):
    #Como las fibras ya están muestreadas con 21 puntos, se evita el resampleo interno de QB...
    feature = IdentityFeature()
    metric = AveragePointwiseEuclideanMetric(feature=feature)
    
    #Umbral de distancia de QB en mm (notar que debe ser un float)
    theta = float(theta_threshold)

    #Creación de método QB
    qb = QuickBundles(threshold = theta, metric = metric)
    
    #Clustering
    t0 = time.time()
    inter_clusters = qb.cluster(centroids_data)
    
    print('Tiempo total interclustering: ' + str(time.time() - t0) + ' [s]')
    
    #Guardado de clustering en pkl
    save_pickle(output_data + 'centroids_QB_' + str(theta_threshold) + 'mm.pkl', inter_clusters)
    
    return inter_clusters

#Carga de los centroides intersujetos
def load_intercentroids (intercentroids_path, theta):
    inter_clusters = load_pickle(intercentroids_path + 'centroids_QB_' + str(theta) + 'mm.pkl')
    
    return inter_clusters

#Calculo de representatividad de sujetos en los clusters intersujetos
def representativity (centroid_path, intercentroids_path, theta):
    #Arreglo para almacenar la cantidad de sujetos que componen cada cluster.
    represent = []
    ids =load_filtered_centroids(centroid_path)[0]
    inter_clusters=load_intercentroids(intercentroids_path, theta)
    
    for cluster in inter_clusters: 
        
        subs = len(set(ids[cluster.indices]))     #Cantidad de ids únicos que aparecen en el cluster ínter-sujeto.
        
        represent.append(subs) #Agrego cantidad total de sujetos que están en el cluster de QB al arreglo represent
    
    represent = np.array(represent)
    
    #Creación de histograma
    hist = [0]*79
    
    for i in range(0,79):
        hist[i] = np.sum(represent == i+1) #sumo la cantidad de veces que hay "i" sujetos en todos los cluster
    
    #Gráfico del histograma    
    fig = plt.figure(figsize = (8,4), dpi = 70)
    ax = fig.add_axes([0,0,1,1])
    ax.bar(range(1,80), hist[:])
    
    plt.xlabel('Número de sujetos', fontsize = 24)
    plt.xticks([1,10,20,30,40,50,60,70,79], fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.ylabel('Número de clústers', fontsize = 24)
    plt.title('Representatividad de clustering ínter-sujeto para '+str(theta)+'[mm]', fontsize = 24)
    
    print('numero de cluster con 70 o más sujetos =', sum((hist[69:79]))) 

#Filtrado de clusters segun umbral de representatividad
def clusters_filtering(subs_quantity,centroid_path, intercentroids_path, theta): 
    keys = range(1,80)  #De 1 hasta 79 (sujetos).
    cont=0
    inter_clusters=load_intercentroids(intercentroids_path, theta)
    ids, intra_centroids_AS =load_filtered_centroids(centroid_path)
    fibers= load_filtered_clusters(centroid_path)
    
    #Diccionarios que almacenan los clusters, nombre de clusters, y centroides para cada sujeto.
    subs_clusters = {key: [] for key in keys}
    subs_centroids = {key: [] for key in keys}
    subs_clusters_names = {key: [] for key in keys}
    
    #Umbral mínimo de representatividad (cantidad mínima de sujetos distintos que componen el cluster).
    subs_thr = subs_quantity
    
    #Índices de clusters ínter-sujeto que superan el umbral de representatividad.
    qb_valid_idx = []
    
    #Para cada cluster ínter-sujeto...
    for j, cluster in enumerate(inter_clusters): #para todos los clusters QB de centroides intrasujetos (cluster intersujetos)
        
        #Sujetos únicos que lo componen
        subs = set(ids[cluster.indices])
        
        #Si la cantidad de sujetos no supera el umbral, descartar el cluster. De lo contrario seguir.
        if len(subs) < subs_thr:
            #print(j)
            cont+=1
            continue
        
        #Se guarda el índice del cluster válido.
        qb_valid_idx.append(j)
        
        #Para cada uno de los sujetos que componen el cluster...
        for sub in subs:
            #Se guarda el nombre del cluster, basado en el índice relativo a los clusters obtenidos con QB.
            subs_clusters_names[sub].append('cluster_' + str(j))
        
        #Se crean diccionarios temporales para cada sujeto que compone el cluster, una para las fibras y otra para los centroides.
        subs_temp = {key: [] for key in subs}
        subs_centroid_temp = {key: [] for key in subs}
        
        #Para cada índice (centroide intra-sujeto) que contiene el cluster...
        for idx in cluster.indices:
            #Se guarda el centroide en la lista de centroides, en el elemento correspondiente al sujeto que originó el centroide.
            subs_centroid_temp[ids[idx]].append(intra_centroids_AS[idx]) # idx tomara los valores de los centroides que componen el cluster j
                                                                        # luego estos valores se refenciarán en los ID registrados.
            #Se hace lo mismo para las fibras que componen el cluster intra-sujeto asociadas al centroide intra-sujeto.
            for fiber in fibers[idx]:
                subs_temp[ids[idx]].append(fiber)
    
        #Finalmente, se guardan los diccionarios temporales en diccionarios definitivos.
        for sub in subs_temp.keys():
            if subs_centroid_temp[sub] != []:
                subs_centroids[sub].append(subs_centroid_temp[sub])
            if subs_temp[sub] != []:
                subs_clusters[sub].append(subs_temp[sub])
            
 
    print('Total cluster QB con menos de '+str(subs_thr)+' sujetos = '+str(cont))
    #print(subs_centroids)
    return subs_clusters_names, subs_centroids, subs_clusters, qb_valid_idx       

#Conteo de centroides restantes    
def count_remaining_centroids (sub_quatity, theta, centroid_path, fibers ):
    inter_clusters=load_intercentroids(centroid_path, theta)
    #fibers= load_filtered_clusters(centroid_path)
    ids = load_filtered_centroids(centroid_path)[0]

    cuales = []
    n_intra_centroids = 0
    n_intra_fibers = 0
    
    for idx, cluster in enumerate(inter_clusters):
        n_subs = len(set(ids[cluster.indices]))
        
        if n_subs >= sub_quatity:
            cuales.append(idx)
            n_intra_centroids += len(cluster.indices)
            for j in range(len(fibers[cluster.indices])):
                n_intra_fibers += len(fibers[cluster.indices][j])

    print('Quedan ' + str(n_intra_centroids) + ' centroides intra-sujeto,')
    print('Asociados a ' + str(n_intra_fibers) + ' fibras intra-sujeto.')

#Guardado de clusters y centroides intersujetos filtrado por cantidad de sujetos que lo componen
def intersubject_clusters_saving(subs_thr,centroid_path, theta):
    #Lista con nombres de clusters inter-sujeto (cluster_ + índice del cluster, p.ej. cluster_15).
    t0 = time.time()
    inter_clusters=load_intercentroids(centroid_path, theta)
    inter_clusters_names = []
    subs_clusters_names, subs_centroids, subs_clusters, qb_valid_idx = clusters_filtering(subs_thr,centroid_path,centroid_path,theta)
    print('filtrado mayor a '+str(subs_thr)+' listo')
    print('Tiempo: ' + str(time.time() - t0) + ' [s]')
    
    for j, cluster in enumerate(inter_clusters):
        if len(cluster) == 1: #
            continue
        inter_clusters_names.append('cluster_' + str(j))
   
    paleta = random_palette(len(inter_clusters_names))
    inter_path='inter-clustered_QB_' + str(theta)+'/'
    
    #Creación y guardado del archivo jerárquico para visualización por clusters en Anatomist.
    if not os.path.exists(centroid_path + inter_path):
            os.makedirs(centroid_path + inter_path)
    write_hie(centroid_path + inter_path + 'QB_' + str(theta) + '.hie', inter_clusters_names, paleta)
        
    for sub in subs_clusters.keys():

        if sub < 10:
            name = '00' + str(sub)
        else:
            name = '0' + str(sub)
        
        if not os.path.exists(centroid_path + inter_path + name):
            os.makedirs(centroid_path + inter_path + name)

        bt3.write_bundle_severalbundles(centroid_path + inter_path + name + '/clusters.bundles', subs_clusters[sub], subs_clusters_names[sub])
        bt3.write_bundle_severalbundles(centroid_path + inter_path + name + '/centroids.bundles', subs_centroids[sub], subs_clusters_names[sub])
        print('Datos intersujetos para sujeto '+name+' guardados')
        
    qb_valid_clusters = []
    qb_valid_centroids = []
    qb_valid_names = []
    
    for idx in qb_valid_idx:
        
        temp = []
        for centroid in inter_clusters[idx]:
            temp.append(centroid)
        qb_valid_clusters.append(temp)
        
        qb_valid_centroids.append([inter_clusters[idx].centroid])
        qb_valid_names.append('cluster_' + str(idx))
        
    bt3.write_bundle_severalbundles(centroid_path + inter_path + 'QB_' + str(theta) + '_centroids.bundles', qb_valid_centroids, qb_valid_names)
    bt3.write_bundle_severalbundles(centroid_path + inter_path + 'QB_' + str(theta) + '_clusters.bundles', qb_valid_clusters, qb_valid_names)
    print('Fibras y centroides de clusters intersujetos guardadas')
    print('Tiempo Total guardado de clusters y centroides intersujetos: ' + str(time.time() - t0) + ' [s]')

#Alineamiento de centroides y guardado de clusters individuales
def cluster_alignment (start, end, centroid_path, theta):
    
    inter_path='inter-clustered_QB_' + str(theta)+'/'
    ref_centroids = bt3.read_bundle_severalbundles(centroid_path + inter_path + 'QB_' + str(theta) + '_centroids.bundles')
    t0 = time.time()
    
    #Para cada sujeto...
    for i in range(start,end+1):
        
        #Listas con nuevos clusters y centroides
        sub_clusters = []
        sub_centroids = []
        
        if i < 10:
            sub = '00' + str(i)
        else:
            sub = '0' + str(i)
        
        #Lectura de clusters de cada sujeto
        clusters_obj = bt3.read_bundle_severalbundles(centroid_path + inter_path + sub + '/clusters.bundles')
        clusters = clusters_obj[0]
        
        #Nombres de los clusters (con respecto a los ínter-sujeto)
        names = clusters_obj[1]
        
        #Lectura de centroides
        centroids = bt3.read_bundle_severalbundles(centroid_path + inter_path + sub + '/centroids.bundles')[0]
    
        #Para cada cluster...
        for idx, cluster in enumerate(clusters):
            
            #Se obtiene el nombre del cluster en cuestión
            cluster_name = names[idx]
            
            #Obtenemos el índice del cluster en la lista de nombres de centroides de referencia...
            ref_idx = ref_centroids[1].index(cluster_name)
            
            #...y con el índice obtenemos el centroide correspondiente (a nivel ínter-sujeto)
            ref = ref_centroids[0][ref_idx][0]
            
            #Tomamos todas las fibras del cluster y las alineamos con respecto al centroide ínter-sujeto (utilizando la distancia MDF)
            cluster_a = align_cluster(cluster, ref)
            
            #Teniendo todas las fibras alineadas, calculamos el nuevo centroide del cluster ínter-sujeto a nivel intra-sujeto (promedio geométrico del aporte del sujeto al cluster ínter-sujeto)
            centroid = calculate_centroid(cluster_a)
            
            sub_clusters.append(cluster_a)
            sub_centroids.append([centroid])
        
        inter_path='inter-clustered_QB_' + str(theta)+'/'
        individual_path= centroid_path + inter_path + sub + '/individual_clusters_aligned/'
        
        if not os.path.exists(individual_path):
            os.makedirs(individual_path)

        #Se guardan los clusters y centroides finales. Ahora Final_Centroids sí contiene sólo un centroide por cluster.
        for j, name in enumerate(names):
            bt3.write_bundle(individual_path + str(name) + '.bundles', sub_clusters[j],names[j])
        
        print(str(len(names))+' Clusters individuales alineados y guardados para sujeto ' + sub)

        
    print('Tiempo total alineamiento de centroides y guardado individual de clusters: ' + str(time.time()- t0) + ' [s]')
  
#Path

path='Archi_y_codigos/ARCHI/subs/'
centroid_path=path+'filtered_centroids6/'
#Parámetros
theta_QB= 21
subs_th= 70 #al menos el 89% de sujetos totales
start= 1
end= 79
#Funciones
#save_filtered_centroids(start, end, path)
# intra_centroids_AS = load_filtered_centroids(centroid_path)[1]
# interclustering(theta_QB, intra_centroids_AS, centroid_path)
#intersubject_clusters_saving(subs_th, centroid_path, theta_QB)
#cluster_alignment(start, end, centroid_path, theta_QB)
