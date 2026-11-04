# -*- coding: utf-8 -*-
"""
Created on Wed May 26 18:15:03 2021

@author: fondecyt-1190701
"""

from BundleTools import bundleTools as bt
from BundleTools import bundleTools3 as bt3
import os
import shutil 
import numpy as np
from numpy.linalg import inv
from subprocess import  call
import time

#Leer archivos de matrices .trm
def readtrm(tr):
	ar=open(tr,'r')
	t=[i.split(' ') for i in ar.read()[0:-1].split('\n')]

	t=np.array(t).astype('float32')
	return t

#Transforma las matrices trm a matriz afin en coordenadas homogeneas
def trm_to_std(trm):
    T0 = trm[0,0]
    T1 = trm[0,1]
    T2 = trm[0,2]
    A0 = trm[1,0]
    A1 = trm[1,1]
    A2 = trm[1,2]
    B0 = trm[2,0]
    B1 = trm[2,1]
    B2 = trm[2,2]
    C0 = trm[3,0]
    C1 = trm[3,1]
    C2 = trm[3,2]
    
    std = np.array(([A0,A1,A2,T0],
                    [B0,B1,B2,T1],
                    [C0,C1,C2,T2],
                    [0,0,0,1]))
    
    return std

#Transforma las matrices afines con coordenadas homogeneas a trm de acuerdo a las convenciones de Anatomist
def std_to_trm(std):
    A0 = std[0,0];
    A1 = std[0,1];
    A2 = std[0,2];
    B0 = std[1,0];
    B1 = std[1,1];
    B2 = std[1,2];
    C0 = std[2,0];
    C1 = std[2,1];
    C2 = std[2,2];
    T0 = std[0,3];
    T1 = std[1,3];
    T2 = std[2,3];
    
    trm = np.array(([T0,T1,T2],
                    [A0,A1,A2],
                    [B0,B1,B2],
                    [C0,C1,C2]))    
    return trm

#Se guarda matriz trm
def save_trm(filename, trm):
    with open(filename,'w') as f:
        f.write("%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f" % (trm[0,0], trm[0,1], trm[0,2],
                                                            trm[1,0], trm[1,1], trm[1,2],
                                                            trm[2,0], trm[2,1], trm[2,2],
                                                            trm[3,0], trm[3,1], trm[3,2]))  

#Borrar archivos y carpetas        
def delete_files_folders (start, end, main_path, target_path):
    for x in range(start,end+1):
     	if x<10:
             name='00'+ str(x)
     	else:
    		 name='0'+ str(x)
     	#shutil.rmtree(path+name+ target_path)
     	os.remove(main_path + name + target_path)

#Unir archivos Bundlesdata
def files_union (start, end, path):
    #path example: path='Archi_y_codigos/resampling/filtering/transform/ffclust/inter-subject_clustering/ARCHI/subs/'
    for x in range(start,end+1):
        if x < 10:
            sub = '00' + str(x)
        else:
            sub = '0' + str(x)
                    
        fibers1=bt.read_bundle(path+sub+'/backup/tractography-streamline-regularized-deterministic_1_3.bundles')
        fibers2=bt.read_bundle(path+sub+'/backup/tractography-streamline-regularized-deterministic_2_3.bundles')
        fibers3=bt.read_bundle(path+sub+'/backup/tractography-streamline-regularized-deterministic_3_3.bundles')
        fibers=fibers1+fibers2+fibers3
        
        if not os.path.exists(path+sub+'/merged/'):
            os.makedirs(path + sub+'/merged/')
        bt3.write_bundle(path+sub+'/merged/tractography_merged.bundles', fibers, "")
        print('Partes de tractografia ' +sub+ ' unida') 

#Crear directorios para preprocesamiento
def preprocessing_directories(start, end, path, folder_name):
    #path example: path='Archi_y_codigos/resampling/filtering/transform/ffclust/inter-subject_clustering/ARCHI/subs/'
    #folder name example: folder_name='resamplig'
    for x in range(start,end+1):
        if x < 10:
            sub = '00' + str(x)
        else:
            sub = '0' + str(x)    
        
        if not os.path.exists(path+sub+'/merged/'+folder_name+'/'):
            os.makedirs(path+sub+'/merged/'+folder_name+'/')

#Aplicar algoritmo de remuestreo de fibras
def resampling (start, end, path, points):
    #path example: path='Archi_y_codigos/resampling/filtering/transform/ffclust/inter-subject_clustering/ARCHI/subs/'
    #example points: points='21'
    t0 = time.time()
    #Si fuese necesario, compilar antes:
    #call(['gcc','resampling/BundleTools_sp.c','resampling/resampling.c','-o','resampling/resampling','-lm'])
    for x in range(start,end+1):
        if (x<10):
            sub='00'+str(x)
        else:
            sub='0'+str(x)
            
        if not os.path.exists(path+sub+'/merged/resampled3/'):
            os.makedirs(path+sub+'/merged/resampled3/')
        name= 'tractography_merged.bundles'
        arg= ['resampling/./resampling', path+sub+'/merged/'+name, path+sub+'/merged/resampled3/tractography_resampled.bundles', points]
        call(arg)
        print('Tractografia '+sub+' remuestreada en: '+points+' puntos')
    print('Tiempo Total ' + str(time.time() - t0) + ' [s]')

#Aplicar algoritmo de filtrado de fibras
def filtering (start, end, path, min_size):
    #path example: path='Archi_y_codigos/resampling/filtering/transform/ffclust/inter-subject_clustering/ARCHI/subs/'
    #min_size example: min_size='40'
    t0 = time.time()
    for x in range(start,end+1):
        if (x<10):
            sub='00'+str(x)
        else:
            sub='0'+str(x)
    
        if not os.path.exists(path+sub+'/merged/filtered3/'):
            os.makedirs(path+sub+'/merged/filtered3/')
        
        name= 'tractography_resampled.bundles'
        arg= ['python','filtering/filtering.py',path+sub+'/merged/resampled3/'+name, path + sub +
              '/merged/filtered3/tractography_filtered.bundles', min_size]
        call(arg)     
        print('Tractografia filtrada: '+sub)
    print('Tiempo Total ' + str(time.time() - t0) + ' [s]')

#Aplicar algoritmo de transformacion de tractografías
def transform (start, end, path, path_matrix):
    #path example: path='Archi_y_codigos/resampling/filtering/transform/ffclust/inter-subject_clustering/ARCHI/subs/'
    #path_matrix example: path_matrix='/transformmatrices/trmT1ToTal.trm'
    t0 = time.time()
    #Si fuese necesario, compilar antes:
    #call(['g++','-std=c++14','-O3','transform/transform.cpp','-o','transform/transform','-fopenmp','-ffast-math'])
    for x in range(start,end+1):
        if (x<10):
            sub='00'+str(x)
        else:
            sub='0'+str(x)
    
        if not os.path.exists(path+sub+'/merged/transformed3/'):
            os.makedirs(path+sub+'/merged/transformed3/')
            
        arg= ['transform/./transform',path + sub + '/merged/filtered3/', path + sub+'/merged/transformed3/',
              path + sub + path_matrix]
        call(arg)   
        print('Tractografia transformada: '+sub)
    print('Tiempo Total ' + str(time.time() - t0) + ' [s]')

def transform2 (start, end, path, path_matrix):
    #path example: path='Archi_y_codigos/resampling/filtering/transform/ffclust/inter-subject_clustering/ARCHI/subs/'
    #path_matrix example: path_matrix='/transformmatrices/trmT1ToTal.trm'
    t0 = time.time()
    #Si fuese necesario, compilar antes:
    #call(['g++','-std=c++14','-O3','transform/transform.cpp','-o','transform/transform','-fopenmp','-ffast-math'])
    for x in range(start,end+1):
        if (x<10):
            sub='00'+str(x)
        else:
            sub='0'+str(x)

        if not os.path.exists(path+sub+'/merged/transformedT1/'):
            os.makedirs(path+sub+'/merged/transformedT1/') 
        arg= ['transform/./transform',path + sub + '/merged/resampled3/', path + sub+'/merged/transformedT1/',
              path + sub + path_matrix]
        call(arg)   
        print('Tractografia transformada: '+sub)
    print('Tiempo Total ' + str(time.time() - t0) + ' [s]')


#Aplicar algoritmo de transformación para cada cluster
def transform_individual_clusters (start, end, path, path_matrix, theta): #!!!!!Modificado por Joaquín

    #path example: path='Archi_y_codigos/resampling/filtering/transform/ffclust/inter-subject_clustering/ARCHI/subs/'
    #path_matrix example: path_matrix='/transformmatrices/trmT1ToTal.trm'
    t0 = time.time()
    #Si fuese necesario, compilar antes:
    #call(['g++','-std=c++14','-O3','transform/transform.cpp','-o','transform/transform','-fopenmp','-ffast-math'])
    for x in range(start,end+1):
        if (x<10):
            sub='00'+str(x)
        else:
            sub='0'+str(x)
        
        path_individual= path +'/filtered_centroids6/inter-clustered_QB_'+str(theta)+'/'+ sub +'/individual_clusters_aligned/'
        if not os.path.exists(path_individual+'T1/'):
            os.makedirs(path_individual+'T1/')
        
        print(path_individual + "T2/")
        print(path_individual+'T1/')
        print(path + "/" + sub + path_matrix)
        print("\n")
        arg= ['transform/./transform', path_individual + "T2/", path_individual+'T1/', path + "/" +sub + path_matrix]
        call(arg)   
        print('Cluster individuales transformados para: '+sub)
    print('Tiempo Total ' + str(time.time() - t0) + ' [s]')

#Aplicar algoritmo de FFClust
def ffclust(start, end, path):
    #path example: path='Archi_y_codigos/resampling/filtering/transform/ffclust/inter-subject_clustering/ARCHI/subs/'
    #Si fuese necesario, compilar antes:
    #call(['gcc','-fPIC','-shared','-O3','-o','ffclust/segmentation_clust_v1.2/segmentation.so','ffclust/segmentation_clust_v1.2/segmentation.c','-fopenmp','-ffast-math'])   
    t0 = time.time()      
    for x in range(start,end+1):
        if (x<10):
            sub='00'+str(x)
        else:
            sub='0'+str(x)
        
        if not os.path.exists(path+sub+'/merged/ffclust3/'):
            os.makedirs(path+sub+'/merged/ffclust3/')
    
        arg= ['python3','ffclust/main.py','--infile', path + sub + '/merged/transformed3/tractography_filtered_trm.bundles', 
               '--outdir', path + sub + '/merged/ffclust4']
        #python3 main.py --infile ARCHI/subs/020/OverSampledFibers/transformed/tractography-streamline-regularized-deterministic_filtered_trm.bundles --outdir ARCHI/subs/020/OverSampledFibers/ffclust 
        # # #call(['gcc','BundleTools_sp.c','resampling.c','-o','resampling','-lm'])
        call(arg)
        print('FFclust aplicado a: '+sub)
    print('Tiempo Total ' + str(time.time() - t0) + ' [s]') 
        
     
#Path    
# path='Archi_y_codigos/ARCHI/subs/'
# path_matrix='/transformmatrices/T2_to_Tal_tr_tmp.trm'
# #Parámetros
# start=1
# end=79
# points='21'
# min_size= '40'
#Funciones
#files_union(start, end, path)
#resampling(start, end, path, points)
#filtering(start, end, path, min_size)
#transform(start, end, path, path_matrix)
#ffclust(start, end, path)
