#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 20:51:47 2021

@author: fondecyt-1190701
"""

import os
import subprocess as sp
from time import time
import numpy as np
from shutil import copyfile, move
import glob
from BundleTools import bundleTools as bt
from BundleTools import bundleTools3 as bt3
import PreProcessingTools as ppt

#Calculo de interseccion de clusters individuales con mallado cortical
def intersection (start, end, folder_path, bundles_path, I_flag, M_flag, R_flag, theta):

    t=time()

    for i in range(start,end+1):     #Iteración sobre todos los sujetos
        #---------- Sujeto a evaluar ----------#
        if i < 10:
            sub = '00' + str(i)    #Los primeros 9 sujetos tienen prefijo "00" (001-009)...

        else:
            sub = '0' + str(i)     #...y el resto tiene prefijo "0" (por ejemplo, 010-079).

        #---------- Setup (mallado y creación de carpetas) ----------#
        Lhemi_path = folder_path + '/lh.obj'
        # Ruta en donde se encuentra el mallado del hemisferio izquierdo (en formato .obj)
        Rhemi_path = folder_path + '/rh.obj'
        # Ruta en donde se encuentra el mallado del hemisferio derecho (en formato .obj)

        intersection= 'intersection/'
        membership= 'membership/'
        intersection_path = 'intersection_results/'+ intersection + sub + '/xd/'
        membership_path= 'intersection_results/'+ membership + sub + '/xd/'

        if not os.path.exists(membership_path):
            os.makedirs(membership_path)
            #Creación de carpeta para guardar las tuplas de conectividad (In,Fn) de cada cluster, en formato .membershipdata.

        if not os.path.exists(intersection_path):
            os.makedirs(intersection_path)
            #Creación de carpeta para guardar la información de intersección de cada cluster, en formato .intersectiondata.

        print('\nProcesando sujeto: ' + sub)

        # ---------- Selección de fascículos a procesar ----------#
        individual_path = bundles_path
        # Ruta en donde se encuentran los clusters del sujeto (en formato .bundles)
        #sp.call(['intersection/make']);  #Compilación del código C

        if R_flag==0:
            #---------- Intersección clusters con hemisferio izquierdo (Lhemi) ----------#
            t1 = time()    #Cálculo de tiempo de ejecución

            #---------- Llamada a código C++ ----------#
            print('\nLhemi direct')
            sp.call(['intersection/./interx', Lhemi_path, Rhemi_path, individual_path, '', intersection_path + 'L-L_direct/', '',
                     membership_path + 'L-L_direct/', '', str(I_flag), str(M_flag), str(R_flag)]) #Llamada al código C
            #Nota: Al llamar al código C en la línea anterior, los argumentos vacíos ('') se deben a que el código original se diseñó para trabajar con ambos hemisferios al mismo tiempo...
            #...sin embargo, esta implementación trabaja los hemisferios por separado, razón por la cual la mitad de los argumentos son vacíos.

            print('Tiempo para Lhemi Direct: ' + str(time()-t1) + '[s]')

            """
            #---------- Intersección clusters con hemisferio derecho (Rhemi) ----------#
            t2 = time()    #Cálculo de tiempo de ejecución

            #---------- Llamada a código C++ ----------#
            print('\nRhemi direct')
            sp.call(['intersection/./interx', Rhemi_path, '', individual_path, '', intersection_path + 'R-R_direct/', '',
                     membership_path + 'R-R_direct/', '', str(I_flag), str(M_flag), str(R_flag)])

            print('Tiempo para RHemi Direct: ' + str(time()-t2) + '[s]')
            """

        elif R_flag==1:
            #---------- Intersección clusters SENTIDO INVERSO con hemisferio izquierdo (Lhemi) ----------#
            t3 = time()    #Cálculo de tiempo de ejecución

            #---------- Llamada a código C++ ----------#
            print('\nLhemi inverse')
            sp.call(['intersection/./interx', Lhemi_path, Rhemi_path, individual_path, '', intersection_path + 'L-L_inverse/', '',
                     membership_path + 'L-L_inverse/', '', str(I_flag), str(M_flag), str(R_flag)])

            print('Tiempo para Lhemi Inverse: ' + str(time()-t3) + '[s]')

            """
            #---------- Intersección clusters SENTIDO INVERSO con hemisferio derecho (Rhemi) ----------#
            t4 = time()    #Cálculo de tiempo de ejecución

            #---------- Llamada a código C++ ----------#
            print('\nRhemi inverse')
            sp.call(['intersection/./interx', Rhemi_path, '', individual_path, '', intersection_path + 'R-R_inverse/', '',
                     membership_path + 'R-R_inverse/', '', str(I_flag), str(M_flag), str(R_flag)])

            print('Tiempo para Rhemi Inverse: ' + str(time()-t4) + '[s]')
            """
    print('Tiempo total: ' + str(time()-t)+'[s]' )

def intersection2 (start, end, folder_path, bundles_path, I_flag, M_flag, R_flag, theta):

    t=time()

    for i in range(start,end+1):     #Iteración sobre todos los sujetos
        #---------- Sujeto a evaluar ----------#
        if i < 10:
            sub = '00' + str(i)    #Los primeros 9 sujetos tienen prefijo "00" (001-009)...

        else:
            sub = '0' + str(i)     #...y el resto tiene prefijo "0" (por ejemplo, 010-079).

        #---------- Setup (mallado y creación de carpetas) ----------#
        Lhemi_path = folder_path + 'meshes_obj/' + sub + '/lh.obj'
        # Ruta en donde se encuentra el mallado del hemisferio izquierdo (en formato .obj)
        Rhemi_path = folder_path + 'meshes_obj/' + sub + '/rh.obj'
        # Ruta en donde se encuentra el mallado del hemisferio derecho (en formato .obj)
        intersection_path = folder_path + sub +'/intersection_all/intersection/'
        membership_path= folder_path +sub+ '/intersection_all/membership/'

        if not os.path.exists(membership_path):
            os.makedirs(membership_path)
            #Creación de carpeta para guardar las tuplas de conectividad (In,Fn) de cada cluster, en formato .membershipdata.

        if not os.path.exists(intersection_path):
            os.makedirs(intersection_path)
            #Creación de carpeta para guardar la información de intersección de cada cluster, en formato .intersectiondata.

        print('\nProcesando sujeto: ' + sub)

        # ---------- Selección de fascículos a procesar ----------#
        individual_path = "/home/joaquin/MT/parcelacion_john/Archi_y_codigos/ARCHI/subs/001/merged/transformedT1/"
        # Ruta en donde se encuentran los clusters del sujeto (en formato .bundles)
        #sp.call(['intersection/make']);  #Compilación del código C

        #LH
        sp.call(['intersection/./interx', Lhemi_path, '', individual_path, '', intersection_path + 'LH/', '',
                 membership_path + 'LH/', '', str(I_flag), str(M_flag), str(R_flag)])
        #RH
        #sp.call(['intersection/./interx', Rhemi_path, '', individual_path, '', intersection_path + 'RH/', '',
                 #membership_path + 'RH/', '', str(I_flag), str(M_flag), str(R_flag)])

    print('Tiempo total: ' + str(time()-t)+'[s]' )


#Clasificación de clusters (intra e inter-hemisferio)
def bundle_classification (start, end, bundles_path, theta):

    t= time()

    for i in range(start, end+1):     #Iteración sobre todos los sujetos

        #---------- Sujeto a evaluar ----------#
        if i < 10:
            sub = '00' + str(i)    #Los primeros 9 sujetos tienen prefijo "00" (001-009)...

        else:
            sub = '0' + str(i)     #...y el resto tiene prefijo "0" (por ejemplo, 010-079).

        print('\nProcesando sujeto ' + sub)

        # #---------- Setup (mallado y creación de carpetas) ----------#
        t1 = time()    #Cálculo de tiempo de ejecución

        individual_path = bundles_path + 'inter-clustered_QB_' + str(theta)+'/' + sub + '/individual_clusters_aligned/T1/' # Ruta en donde se encuentran los clusters del sujeto (en formato .bundles)
        bundle_classification='bundle_classification_QB_'+str(theta)+'/'
        intersection='intersection_QB_'+str(theta)+'/'
        classification_path = bundles_path + 'intersection/'+ bundle_classification + sub + '/' # Ruta en donde se guardarán los bundles clasificados según su tipo de intersección.
        cluster_names = [os.path.basename(os.path.splitext(x)[0]) for x in glob.glob(individual_path + '*.bundles')] # Lista con nombres de los clusters a evaluar

        membership='membership_QB_'+str(theta)+'/'
        Ldirect_path = bundles_path + 'intersection/'+ membership + sub + '/L-L_direct/' # Ruta de tupla de conectividad en el caso Ldirect.
        Linverse_path = bundles_path + 'intersection/'+ membership + sub + '/L-L_inverse/' # Ruta de tupla de conectividad en el caso Linverse.
        Rdirect_path = bundles_path + 'intersection/'+ membership + sub + '/R-R_direct/' # Ruta de tupla de conectividad en el caso Rdirect.
        Rinverse_path = bundles_path + 'intersection/'+ membership + sub + '/R-R_inverse/' # Ruta de tupla de conectividad en el caso Rinverse.

        # #---------- Creación de diccionario ----------# #
        memb_dict = {} # Diccionario que contiene la tupla de conectividad de cada cluster.

        for name in cluster_names:  # Iteración sobre cada cluster.
            Ldirect = Ldirect_path + name + '.membershipdata'
            Linverse = Linverse_path + name + '.membershipdata'
            Rdirect = Rdirect_path + name + '.membershipdata'
            Rinverse = Rinverse_path + name + '.membershipdata'

            with open(Ldirect, 'rb') as f1, open(Linverse, 'rb') as f2, open(Rdirect, 'rb') as f3, open(Rinverse, 'rb') as f4:  # Lectura de tupla de conectividad para los 4 casos posibles.
                Ldirect_flags = tuple(f1.read( 2 ))
                Linverse_flags = tuple(f2.read( 2 ))
                Rdirect_flags = tuple(f3.read( 2 ))
                Rinverse_flags = tuple(f4.read( 2 ))

            memb_dict[name] = {'Ldirect': Ldirect_flags, 'Linverse': Linverse_flags, 'Rdirect': Rdirect_flags, 'Rinverse': Rinverse_flags} # Se añaden las tuplas al diccionario.

        # #---------- Creación de directorios de clasificación ----------# #
        paths=['L-L/', 'R-R/', 'L-R/', 'R-L/', 'Discarded/', 'Ambiguous/']
        for i in range(len(paths)):
            if not os.path.exists(classification_path + paths[i]):
                os.makedirs(classification_path + paths[i])

        ambiguous = [] # Lista de clusters ambiguos.

        # # ---Bundle Clasification--- # #
        for cluster in memb_dict:   # Iteración sobre cada cluster.
            Ldirect_flags = memb_dict[cluster]['Ldirect']  # Tupla de conectividad Ldirect.
            Linverse_flags = memb_dict[cluster]['Linverse']  # Tupla de conectividad Linverse.
            Rdirect_flags = memb_dict[cluster]['Rdirect']  # Tupla de conectividad Rdirect.
            Rinverse_flags = memb_dict[cluster]['Rinverse']  # Tupla de conectividad Rinverse.

            count = sum(Ldirect_flags) + sum(Linverse_flags) + sum(Rdirect_flags) + sum(Rinverse_flags)    # Cantidad de flags iguales a 1 (se utiliza en el quinto caso).

            #--- Clasificación de clusters ---#
            if (Rdirect_flags == (1,1) or Rinverse_flags == (1,1)) and (Ldirect_flags == (0,0) and Linverse_flags == (0,0)):
                # Clusters que sólo conectan con el hemisferio derecho, en cualquier sentido.
                copyfile(individual_path + cluster + '.bundles', classification_path + 'R-R/' + cluster + '.bundles')
                copyfile(individual_path + cluster + '.bundlesdata', classification_path + 'R-R/' + cluster + '.bundlesdata')

            elif (Ldirect_flags == (1,1) or Linverse_flags == (1,1)) and (Rdirect_flags == (0,0) and Rinverse_flags == (0,0)):
                # Clusters que sólo conectan con el hemisferio izquierdo, en cualquier sentido.
                copyfile(individual_path + cluster + '.bundles', classification_path + 'L-L/' + cluster + '.bundles')
                copyfile(individual_path + cluster + '.bundlesdata', classification_path + 'L-L/' + cluster + '.bundlesdata')

            elif (Rdirect_flags == (1,0) and Linverse_flags == (1,0)) and (Rinverse_flags ==(0,0) and Ldirect_flags == (0,0)):
                copyfile(individual_path + cluster + '.bundles', classification_path + 'R-L/' + cluster + '.bundles')
                copyfile(individual_path + cluster + '.bundlesdata', classification_path + 'R-L/' + cluster + '.bundlesdata')

            elif (Ldirect_flags == (1,0) and Rinverse_flags == (1,0)) and (Linverse_flags ==(0,0) and Rdirect_flags == (0,0)):
                copyfile(individual_path + cluster + '.bundles', classification_path + 'L-R/' + cluster + '.bundles')
                copyfile(individual_path + cluster + '.bundlesdata', classification_path + 'L-R/' + cluster + '.bundlesdata')

            elif count <= 1:    # Si el perfil tiene un flag igual a 1 o menos, es porque no conecta, o probablemente sólo conecte en un extremo.
                copyfile(individual_path + cluster + '.bundles', classification_path + 'Discarded/' + cluster + '.bundles')
                copyfile(individual_path + cluster + '.bundlesdata', classification_path + 'Discarded/' + cluster + '.bundlesdata')

            else:   # Si no cumple cualquiera de los casos anteriores, se considera ambiguo y se necesita más procesamiento.
                copyfile(individual_path + cluster + '.bundles', classification_path + 'Ambiguous/' + cluster + '.bundles')
                copyfile(individual_path + cluster + '.bundlesdata', classification_path + 'Ambiguous/' + cluster + '.bundlesdata')
                ambiguous.append(cluster)

        for cluster in ambiguous:
            Ld = memb_dict[cluster]['Ldirect'] #Perfil de conectividad en caso Ldirect.
            Li = memb_dict[cluster]['Linverse'] #Perfil de conectividad en caso Linverse.
            Rd = memb_dict[cluster]['Rdirect'] #Perfil de conectividad en caso Rdirect.
            Ri = memb_dict[cluster]['Rinverse'] #Perfil de conectividad en caso Rinverse.

            # Lectura de datos de intersección en los cuatro casos posibles (triángulos intersectados, fibras que intersectan).
            Ld_int = bundles_path + 'intersection/'+ intersection + sub + '/L-L_direct/' + cluster + '.intersectiondata'
            Li_int = bundles_path + 'intersection/'+ intersection + sub + '/L-L_inverse/' + cluster + '.intersectiondata'
            Rd_int = bundles_path + 'intersection/'+ intersection + sub + '/R-R_direct/' + cluster + '.intersectiondata'
            Ri_int = bundles_path + 'intersection/'+ intersection + sub + '/R-R_inverse/' + cluster + '.intersectiondata'

            # Variables de conteo para cada caso (intra o inter, y L o R).
            intra_Ld = 0
            intra_Li = 0
            intra_Rd = 0
            intra_Ri = 0
            inter_LR = 0
            inter_RL = 0

            # Dependiendo del perfil de conectividad, obtenemos los datos de intersección de cada uno.
            if Ld == (1,1): # Caso: posiblemente intra hemisferio izquierdo.
                Ld_FnTri = bt.read_intersection(Ld_int)[1]
                intra_idx = np.where(Ld_FnTri != 1000000)  # No consideramos las fibras ínter-hemisferio en un cluster intra-hemisferio.
                intra_Ld = len(Ld_FnTri[intra_idx])    #Cantidad de fibras intra-hemisferio.

            if Li == (1,1): # Caso: posiblemente intra hemisferio izquierdo.
                Li_FnTri = bt.read_intersection(Li_int)[1]
                intra_idx = np.where(Li_FnTri != 1000000)  # No consideramos las fibras ínter-hemisferio en un cluster intra-hemisferio.
                intra_Li = len(Li_FnTri[intra_idx])    #Cantidad de fibras intra-hemisferio.

            if Rd == (1,1): # Caso: posiblemente intra hemisferio derecho.
                Rd_FnTri = bt.read_intersection(Rd_int)[1]
                intra_idx = np.where(Rd_FnTri != 1000000)  # No consideramos las fibras ínter-hemisferio en un cluster intra-hemisferio.
                intra_Rd = len(Rd_FnTri[intra_idx])    #Cantidad de fibras intra-hemisferio.

            if Ri == (1,1): # Caso: posiblemente intra hemisferio derecho.
                Ri_FnTri = bt.read_intersection(Ri_int)[1]
                intra_idx = np.where(Ri_FnTri != 1000000)  # No consideramos las fibras ínter-hemisferio en un cluster intra-hemisferio.
                intra_Ri = len(Ri_FnTri[intra_idx])    #Cantidad de fibras intra-hemisferio.

            if (Ld == (1,1) or Ld == (1,0)) and (Ri == (1,1) or Ri == (1,0)): # Caso: posiblemente ínter hemisferio de izquierda a derecha.
                Ld_FnTri = bt.read_intersection(Ld_int)[1] #Índices de triángulos en el extremo final.

                fib_direct = bt.read_intersection(Ld_int)[-1]  # Índices de fibras que intersectaron en el sentido directo (provenientes del hemisferio izquierdo).
                fib_inverse = bt.read_intersection(Ri_int)[-1]  # Índices de fibras que intersectaron en el sentido inverso (que llegan al hemisferio derecho).

                fib_direct_inter = fib_direct[np.where(Ld_FnTri == 1000000)]    # Índices de fibras que sólo intersectaron en su extremo inicial en sentido directo (potencialmente ínter-hemisferio).

                inter_LR = sum(np.isin(fib_inverse, fib_direct_inter))  # Cantidad de fibras en sentido inverso que coinciden con las intersectadas en sentido directo (es decir, fibras potencialmente inter izquierda a derecha).

            if (Rd == (1,1) or Rd == (1,0)) and (Li == (1,1) or Li == (1,0)): # Caso: posiblemente ínter hemisferio de derecha a izquierda.
                Rd_FnTri = bt.read_intersection(Rd_int)[1] #Índices de triángulos en el extremo final.

                fib_direct = bt.read_intersection(Rd_int)[-1]  # Índices de fibras que intersectaron en el sentido directo (provenientes del hemisferio derecho).
                fib_inverse = bt.read_intersection(Li_int)[-1]  # Índices de fibras que intersectaron en el sentido inverso (que llegan al hemisferio izquierdo).

                fib_direct_inter = fib_direct[np.where(Rd_FnTri == 1000000)]    # Índices de fibras que sólo intersectaron en su extremo inicial en sentido directo (potencialmente ínter-hemisferio).

                inter_RL = sum(np.isin(fib_inverse, fib_direct_inter))  # Cantidad de fibras en sentido inverso que coinciden con las intersectadas en sentido directo (es decir, fibras potencialmente inter derecha a izquierda).

            dom = max(intra_Ld,intra_Li,intra_Rd,intra_Ri,inter_LR,inter_RL)   # Máximo de cantidad de fibras en todos los casos considerados.
            dom_idx = np.argmax([intra_Ld,intra_Li,intra_Rd,intra_Ri,inter_LR,inter_RL])   # Índice del máximo (argmax).

            #No intersection
            if dom == 0:    # Si no hay máximo, ninguno es mayor, y por tanto no es claro el tipo de intersección => se descarta.
                move(classification_path + 'Ambiguous/' + cluster + '.bundles', classification_path + 'Discarded/' + cluster + '.bundles')
                move(classification_path + 'Ambiguous/' + cluster + '.bundlesdata', classification_path + 'Discarded/' + cluster + '.bundlesdata')

            #L_direct
            elif dom_idx == 0:  # Si argmax es el caso intra Ldirect, se clasifica como intra hemisferio izquierdo.
                move(classification_path + 'Ambiguous/' + cluster + '.bundles', classification_path + 'L-L/' + cluster + '.bundles')
                move(classification_path + 'Ambiguous/' + cluster + '.bundlesdata', classification_path + 'L-L/' + cluster + '.bundlesdata')
                #os.remove(Li_int)

            #L_inverse
            elif dom_idx == 1:  # Si argmax es el caso intra Linverse, se clasifica como intra hemisferio izquierdo.
                move(classification_path + 'Ambiguous/' + cluster + '.bundles', classification_path + 'L-L/' + cluster + '.bundles')
                move(classification_path + 'Ambiguous/' + cluster + '.bundlesdata', classification_path + 'L-L/' + cluster + '.bundlesdata')
                #os.remove(Ld_int)

            #R_direct
            elif dom_idx == 2:  # Si argmax es el caso intra Rdirect, se clasifica como intra hemisferio derecho.
                move(classification_path + 'Ambiguous/' + cluster + '.bundles', classification_path + 'R-R/' + cluster + '.bundles')
                move(classification_path + 'Ambiguous/' + cluster + '.bundlesdata', classification_path + 'R-R/' + cluster + '.bundlesdata')
                #os.remove(Ri_int)

            #R_inverse
            elif dom_idx == 3:  # Si argmax es el caso intra Rinverse, se clasifica como intra hemisferio derecho.
                move(classification_path + 'Ambiguous/' + cluster + '.bundles', classification_path + 'R-R/' + cluster + '.bundles')
                move(classification_path + 'Ambiguous/' + cluster + '.bundlesdata', classification_path + 'R-R/' + cluster + '.bundlesdata')
                #os.remove(Rd_int)

            #L-R
            elif dom_idx == 4:  # Si argmax es el caso inter LR, se clasifica como inter hemisferio izquierdo a derecho.
                move(classification_path + 'Ambiguous/' + cluster + '.bundles', classification_path + 'L-R/' + cluster + '.bundles')
                move(classification_path + 'Ambiguous/' + cluster + '.bundlesdata', classification_path + 'L-R/' + cluster + '.bundlesdata')

            #R-L
            elif dom_idx == 5:  # Si argmax es el caso inter RL, se clasifica como inter hemisferio derecho a izquierdo.
                move(classification_path + 'Ambiguous/' + cluster + '.bundles', classification_path + 'R-L/' + cluster + '.bundles')
                move(classification_path + 'Ambiguous/' + cluster + '.bundlesdata', classification_path + 'R-L/' + cluster + '.bundlesdata')

        print('\nTiempo de ejecución: ' + str(time()-t1) + '[s] para sujeto: '+str(sub))

    print('\nTiempo Total de ejecución: ' + str(time()-t) + '[s]')

#Calculo de intersecciones finales
def intersection_classification(start, end, bundles_path, theta):

    t= time()

    #Iteración sobre todos los sujetos
    for i in range (start, end+1):
        #---------- Sujeto a evaluar ----------#
        if i < 10:
            sub = '00' + str(i)    #Los primeros 9 sujetos tienen prefijo "00" (001-009)...

        else:
            sub = '0' + str(i)     #...y el resto tiene prefijo "0" (por ejemplo, 010-079).

        print('\nProcesando sujeto ' + sub)

        t1 = time()

        intersection= 'intersection_QB_'+str(theta)+'/'
        bundle_clasification= 'bundle_classification_QB_'+str(theta)+'/'

        #---------- Perfiles de conectividad ----------#
        Ldirect_path = bundles_path +'intersection/' + intersection + sub + '/L-L_direct/'
        Linverse_path = bundles_path +'intersection/' + intersection + sub + '/L-L_inverse/'
        Rdirect_path = bundles_path +'intersection/' + intersection + sub + '/R-R_direct/'
        Rinverse_path = bundles_path +'intersection/' + intersection + sub + '/R-R_inverse/'

        #---------- Creación de directorios para intersecciones finales ----------#
        paths=['/Final/L-L/', '/Final/L-R/', '/Final/R-L/', '/Final/R-R/', '/Final/Inter/']
        for i in range(len(paths)):
            if not os.path.exists(bundles_path +'intersection/' + intersection + sub + paths[i]):
                os.makedirs(bundles_path +'intersection/' + intersection + sub + paths[i])

        #---------- Clasificaciones iniciales de clusters (obtenidas con bundle_classification) ----------#
        LL_path = bundles_path +'intersection/' + bundle_clasification + sub + '/L-L/'
        LR_path = bundles_path +'intersection/' + bundle_clasification + sub + '/L-R/'
        RL_path = bundles_path +'intersection/' + bundle_clasification + sub + '/R-L/'
        RR_path = bundles_path +'intersection/' + bundle_clasification + sub + '/R-R/'

        #---------- Lectura de bundles en cada categoría ----------#
        LL_list = [os.path.basename(os.path.splitext(x)[0]) for x in glob.glob(LL_path + '*.bundles')] #Todos los clusters que son intra L hemi.
        LR_list = [os.path.basename(os.path.splitext(x)[0]) for x in glob.glob(LR_path + '*.bundles')] #Todos los clusters que son inter L to R hemi.
        RL_list = [os.path.basename(os.path.splitext(x)[0]) for x in glob.glob(RL_path + '*.bundles')] #Todos los clusters que son inter R to L hemi.
        RR_list = [os.path.basename(os.path.splitext(x)[0]) for x in glob.glob(RR_path + '*.bundles')] #Todos los clusters que son intra R hemi.

        #---------- Clusters con resultados extraños (contradictorios o ambiguos) ----------#
        LL_weird = []
        LR_weird = []
        RL_weird = []
        RR_weird = []

        #---------- Clasificación de los clusters intra L hemi ----------#
        #La idea aquí es filtrar las fibras que sólo intersectan en su extremo inicial, y mantener las fibras intra-hemisferio.
        for cluster in LL_list:
            Ld_infile = Ldirect_path + cluster + '.intersectiondata'  #Intersección en sentido directo.
            Li_infile = Linverse_path + cluster + '.intersectiondata'  #Intersección en sentido inverso.

            # En primera instancia trabajamos con intersección en sentido directo
            if os.path.exists(Ld_infile):
                InTri, FnTri, InPoints, FnPoints, fib_index = bt.read_intersection(Ld_infile)
                intra_idx = np.where(FnTri != 1000000)   #Índices de fibras intra-hemisferio (las fibras que sólo conectan en un extremo tienen asignado FnTri = 1000000).
                InTri = InTri[intra_idx]                 #Triángulos del extremo inicial correspondientes sólo a fibras intra-hemisferio.
                n_tri = np.uint32(len(InTri))            #Cantidad de triángulos intra.

                #En caso de que no hayan triángulos intra, intentamos lo mismo con el sentido inverso.
                if n_tri == 0:
                    InTri, FnTri, InPoints, FnPoints, fib_index = bt.read_intersection(Li_infile)
                    intra_idx = np.where(FnTri != 1000000)
                    InTri = InTri[intra_idx]
                    n_tri = np.uint32(len(InTri))

                #En caso de que sí hayan triángulos intra...
                else:
                    FnTri = FnTri[intra_idx]           #Triángulos del extremo inicial correspondientes sólo a fibras intra-hemisferio.
                    InPoints = InPoints[intra_idx]     #Puntos del extremo inicial correspondientes sólo a fibras intra-hemisferio.
                    FnPoints = FnPoints[intra_idx]     #Puntos del extremo final correspondientes sólo a fibras intra-hemisferio.
                    fib_index = fib_index[intra_idx]   #Índices de las fibras intra-hemisferio (relativas al total de fibras del cluster).

                f = open(Ldirect_path + '../Final/L-L/' + cluster + '.intersectiondata', 'wb+')    #Archivo de intersección final a escribir.

                #El archivo .intersectiondata guarda lo siguiente:
                f.write(n_tri)     #Número de triángulos intersectados.

                f.write(InTri)     #Triángulos del extremo inicial.
                f.write(FnTri)     #Triángulos del extremo final.

                f.write(InPoints)  #Puntos del extremo inicial.
                f.write(FnPoints)  #Puntos del extremo final.

                f.write(fib_index) #Índices de fibras intra-hemisferio.

                f.close()

            #Si no hay intersección en sentido directo, utilizamos la encontrada en sentido inverso.
            elif os.path.exists(Li_infile):
                InTri, FnTri, InPoints, FnPoints, fib_index = bt.read_intersection(Li_infile)
                intra_idx = np.where(FnTri != 1000000)
                InTri = InTri[intra_idx]
                n_tri = np.uint32(len(InTri))

                #Si no hay triángulos, entonces descartamos el cluster (aunque esta situación no debiese ocurrir).
                if n_tri == 0:
                    move(bundles_path +'intersection/' + bundle_clasification + sub + '/L-L/' + cluster + '.bundles',
                         bundles_path +'intersection/' + bundle_clasification + sub + '/Discarded/' + cluster + '.bundles')
                    move(bundles_path +'intersection/' + bundle_clasification + sub + '/L-L/' + cluster + '.bundlesdata',
                         bundles_path +'intersection/' + bundle_clasification + sub + '/Discarded/' + cluster + '.bundlesdata')
                    continue

                else:
                    FnTri = FnTri[intra_idx]
                    InPoints = InPoints[intra_idx]
                    FnPoints = FnPoints[intra_idx]
                    fib_index = fib_index[intra_idx]

                f = open(Ldirect_path + '../Final/L-L/' + cluster + '.intersectiondata', 'wb+')

                f.write(n_tri)

                f.write(FnTri)
                f.write(InTri)

                f.write(FnPoints)
                f.write(InPoints)

                f.write(fib_index)

                f.close()

        #---------- Clasificación de los clusters inter L to R hemi ----------#
        #La idea aquí es filtrar las fibras que intersectan en ambos extremos (intra), y mantener las fibras inter-hemisferio.
        for cluster in LR_list:

            #Tomamos el sentido directo, es decir, de L a R.
            infile = Ldirect_path + cluster + '.intersectiondata'

            #El sentido directo debiese existir, por lo que si no existe lo agregamos a la lista de clusters extraños para inspección.
            if not os.path.exists(infile):
                LR_weird.append(cluster)
                continue

            #Lectura de intersección en sentido directo.
            InTri_direct, FnTri_direct, InPoints_direct, FnPoints_direct, fib_index_direct = bt.read_intersection(infile)

            #Tomamos el sentido inverso, es decir, de R a L.
            infile = Rinverse_path + cluster + '.intersectiondata'

            #Si no existe lo agregamos a la lista de extraños.
            if not os.path.exists(infile):
                LR_weird.append(cluster)
                continue

            #Lectura de intersección en sentido inverso.
            InTri_reverse, FnTri_reverse, InPoints_reverse, FnPoints_reverse, fib_index_reverse = bt.read_intersection(infile)

            #Obtenemos los índices de las fibras que están tanto en sentido directo como en sentido inverso (que aparecen en ambos archivos .intersectiondata).
            fib_index_direct = fib_index_direct[np.in1d(fib_index_direct, fib_index_reverse)]
            fib_index_reverse = fib_index_reverse[np.in1d(fib_index_reverse, fib_index_direct)]

            #Obtenemos los triángulos del extremo inicial en ambos sentidos (donde los triángulos del extremo inicial en sentido inverso, serán los triángulos finales al considerar las fibras como ínter-hemisferio).
            InTri_direct = InTri_direct[fib_index_direct.argsort()]
            InTri_reverse = InTri_reverse[fib_index_reverse.argsort()]

            #Reordenamos los índices.
            fib_index_direct = np.sort(fib_index_direct)

            #Cantidad de triángulos obtenidos.
            n_tri = np.uint32(len(InTri_direct))

            #Si no quedan triángulos, descartamos el cluster.
            if n_tri == 0:
                move(bundles_path +'intersection/' + bundle_clasification + sub + '/L-R/' + cluster + '.bundles',
                     bundles_path +'intersection/' + bundle_clasification + sub + '/Discarded/' + cluster + '.bundles')
                move(bundles_path +'intersection/' + bundle_clasification + sub + '/L-R/' + cluster + '.bundlesdata',
                     bundles_path +'intersection/' + bundle_clasification + sub + '/Discarded/' + cluster + '.bundlesdata')
                continue

            #Guardado de archivo .intersectiondata final.
            f = open(Ldirect_path + '../Final/Inter/' + cluster + '.intersectiondata', 'wb+')

            f.write(n_tri)

            f.write(InTri_direct)
            f.write(InTri_reverse)

            f.write(InPoints_direct)
            f.write(InPoints_reverse)

            fib_index_direct = np.sort(fib_index_direct)
            f.write(fib_index_direct)

            f.close()

        #---------- Clasificación de los clusters inter R to L hemi ----------#
        #Se realiza el mismo procedimiento que para L to R, pero a la inversa.
        for cluster in RL_list:
            infile = Rdirect_path + cluster + '.intersectiondata'

            if not os.path.exists(infile):
                RL_weird.append(cluster)
                continue

            InTri_direct, FnTri_direct, InPoints_direct, FnPoints_direct, fib_index_direct = bt.read_intersection(infile)

            infile = Linverse_path + cluster + '.intersectiondata'

            if not os.path.exists(infile):
                RL_weird.append(cluster)
                continue

            InTri_reverse, FnTri_reverse, InPoints_reverse, FnPoints_reverse, fib_index_reverse = bt.read_intersection(infile)

            fib_index_direct = fib_index_direct[np.in1d(fib_index_direct, fib_index_reverse)]
            fib_index_reverse = fib_index_reverse[np.in1d(fib_index_reverse, fib_index_direct)]

            InTri_direct = InTri_direct[fib_index_direct.argsort()]
            InTri_reverse = InTri_reverse[fib_index_reverse.argsort()]

            fib_index_direct = np.sort(fib_index_direct)

            n_tri = np.uint32(len(InTri_direct))

            if n_tri == 0:
                move(bundles_path +'intersection/' + bundle_clasification + sub + '/R-L/' + cluster + '.bundles',
                     bundles_path +'intersection/' + bundle_clasification + sub + '/Discarded/' + cluster + '.bundles')
                move(bundles_path +'intersection/' + bundle_clasification + sub + '/R-L/' + cluster + '.bundlesdata',
                     bundles_path +'intersection/' + bundle_clasification + sub + '/Discarded/' + cluster + '.bundlesdata')
                continue;

            f = open(Ldirect_path + '../Final/Inter/' + cluster + '.intersectiondata', 'wb+')

            f.write(n_tri)

            f.write(InTri_reverse)
            f.write(InTri_direct)

            f.write(InPoints_reverse)
            f.write(InPoints_direct)

            fib_index_direct = np.sort(fib_index_direct)
            f.write(fib_index_direct)

            f.close()

        #---------- Clasificación de los clusters intra R hemi ----------#
        #Se realiza el mismo procedimiento que para intra L hemi.
        for cluster in RR_list:
            Rd_infile = Rdirect_path + cluster + '.intersectiondata'
            Ri_infile = Rinverse_path + cluster + '.intersectiondata'

            if os.path.exists(Rd_infile):
                InTri, FnTri, InPoints, FnPoints, fib_index = bt.read_intersection(Rd_infile)
                intra_idx = np.where(FnTri != 1000000)
                InTri = InTri[intra_idx]
                n_tri = np.uint32(len(InTri))

                if n_tri == 0:
                    print(cluster)
                    InTri, FnTri, InPoints, FnPoints, fib_index = bt.read_intersection(Ri_infile)
                    intra_idx = np.where(FnTri != 1000000)
                    InTri = InTri[intra_idx]
                    n_tri = np.uint32(len(InTri))

                else:
                    FnTri = FnTri[intra_idx]
                    InPoints = InPoints[intra_idx]
                    FnPoints = FnPoints[intra_idx]
                    fib_index = fib_index[intra_idx]

                f = open(Ldirect_path + '../Final/R-R/' + cluster + '.intersectiondata', 'wb+')

                f.write(n_tri)

                f.write(InTri)
                f.write(FnTri)

                f.write(InPoints)
                f.write(FnPoints)

                f.write(fib_index)

                f.close()

            elif os.path.exists(Ri_infile):
                InTri, FnTri, InPoints, FnPoints, fib_index = bt.read_intersection(Ri_infile)
                intra_idx = np.where(FnTri != 1000000)
                InTri = InTri[intra_idx]
                n_tri = np.uint32(len(InTri))

                if n_tri == 0:
                    move(bundles_path +'intersection/' + bundle_clasification + sub + '/R-R/' + cluster + '.bundles',
                         bundles_path +'intersection/' + bundle_clasification + sub + '/Discarded/' + cluster + '.bundles')
                    move(bundles_path +'intersection/' + bundle_clasification + sub + '/R-R/' + cluster + '.bundlesdata',
                         bundles_path +'intersection/' + bundle_clasification + sub + '/Discarded/' + cluster + '.bundlesdata')
                    continue;

                else:
                    FnTri = FnTri[intra_idx]
                    InPoints = InPoints[intra_idx]
                    FnPoints = FnPoints[intra_idx]
                    fib_index = fib_index[intra_idx]

                f = open(Ldirect_path + '../Final/R-R/' + cluster + '.intersectiondata', 'wb+')

                f.write(n_tri)

                f.write(FnTri)
                f.write(InTri)

                f.write(FnPoints)
                f.write(InPoints)

                f.write(fib_index)

                f.close()

        print('Tiempo de ejecución: ' + str(time()-t1) + '[s]')

    print('Tiempo Total: ' + str(time()-t)+'[s]')


#Path
# folders_path='Archi_y_codigos/ARCHI/subs/'
# interclusters_path=folders_path + 'filtered_centroids6/'
# #Parámetros
# theta_QB=21
# start=1
# end=79
# # i_flag= 1 #intersection flag
# # m_flag= 1 #membership flag
# # r_flag1= 0 #reverse flag (Si y no)
# # r_flag2= 1
# # #Preprocessingtools
# path_matrix_inter='/transformmatrices/trmTalToT1.trm'
# ppt.transform_individual_clusters(start, end, folders_path, path_matrix_inter, theta_QB)

# #Funciones
# intersection(start, end, folders_path, interclusters_path, i_flag, m_flag, r_flag1, theta_QB)
# intersection(start, end, folders_path, interclusters_path, i_flag, m_flag, r_flag2, theta_QB)
# bundle_classification(start, end, interclusters_path, theta_QB)
# intersection_classification(start, end, interclusters_path, theta_QB)
