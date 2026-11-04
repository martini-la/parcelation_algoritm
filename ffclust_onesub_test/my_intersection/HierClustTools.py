
import sys, os
import numpy as np
from time import time
import random
import nipy.algorithms.graph.graph as FG
import nipy.algorithms.clustering.hierarchical_clustering as HC
import bundleTools3 as BT3
#%%
"""Leer grafo y hacer el árbol"""

def get_tree(graphfile):

    #read graph file and create fff graph
    
    #open file
    inFile = open(graphfile, 'r')
    #read vertices number and edges number
    line = inFile.readline()
    l = line.split()
    V = int(l[0])
    E = int(l[1])
    
    #create 2D arrays for edges and weights
    edges = np.zeros([E, 2], 'i')
    weights = np.zeros([E, 1], 'f')
    print('edges: ', E, ' vertices: ', V)
    
    
    #read file lines (edges)
    count = 0
    while True:
    	line = inFile.readline()
    	if line == '': break
    	l = line.split()
    	edges[count, 0] = int(l[0])
    	edges[count, 1] = int(l[1])
    	weights[count] = float(l[2])
    	count = count + 1
    
    #close file
    inFile.close()
    
    #create fff graph (affinity graph to be clusterized)
    G = FG.WeightedGraph(V, edges, weights)
    
    # G.show(X=None, ax=None)
    
    
    #print graph values
    G.V #vertices num
    G.E #edges num
    ed = G.get_edges() #edges
    we = G.get_weights() #weights
    
    #clusterize graph using hierarchical clustering
    #t: Weigthed Forest
    t = HC.average_link_graph(G)
    t.plot()
    
    return t

def read_mesh_obj( infile ):
    vertex = []
    polygons = []
    len_vertex = 0
    len_polygons = 0
    
    with open(infile, 'r') as f:

        for l in f:
            e = l.split()
            
            if len(e) == 0:
                continue
            else:
                ID = e[0]
            
            if ID == 'v':
                vertex.append(e[1:])
                len_vertex += 1
                ID = 'n'
                
            elif ID == 'f':
                polygons.append([int(e[1].split('//')[0])-1, int(e[2].split('//')[0])-1, int(e[3].split('//')[0])-1])
                len_polygons += 1
                ID = 'n'
                
    vertex = np.asarray(vertex, dtype = np.float32)
    polygons = np.asarray(polygons, dtype = int)
    
    return vertex, polygons


def savewforest_txt(forestfile,f):
    lines=[]
    lines.append(str(len(f.parents))+'\n')
    for i in range(len(f.parents)):
        if f.children[i]==[]:
            child=[-1,-1]
        else:
            child=f.children[i]
        lines.append( str(f.parents[i])+' '+str(child[0])+' '+str(child[1])+'\n')
    ar=open(forestfile,'wt')
    ar.writelines(lines);ar.close()


"""asignar nuevos clusters a cada vertice"""
def new_clusters_x_vertex(labels, clusteres_x_dist_max,labels_arr):
    # reformar la lista de new_clusters
    new_clusters_v2=[0]*len(labels)
    for new,ant in clusteres_x_dist_max.items():
        for anterior in ant:
            new_clusters_v2[anterior]=new
            
    # asignar nuevos clusters
    count=0
    for p in range(len(labels)):
        if count==0:
            new_labels_arr_lh=np.where(labels_arr!=labels[p],labels_arr,new_clusters_v2[p])
            count=1
        else:
            new_labels_arr_lh=np.where(new_labels_arr_lh!=labels[p],new_labels_arr_lh,new_clusters_v2[p])

    return new_clusters_v2, new_labels_arr_lh

# """formar las nuevas parcelas"""
# def new_parcels_c_mem(new_clusters, polygons, new_labels_arr):
#     dicti={}
#     dicti_ID_tri={}
#     for s in np.unique(new_clusters):
#         dicti[str(s)]=np.zeros((0,3))
#         dicti_ID_tri[str(s)]=[]
    
#     lista_ID_tri_negros=[]
#     Lpolygons_etiq=[]
#     for h,i in enumerate(polygons):
#         label_a=new_labels_arr[i[0]]
#         label_b=new_labels_arr[i[1]]
#         label_c=new_labels_arr[i[2]]  
#         if label_a==label_b & label_a==label_c & label_a!=-1:    #solo guardo en el dicc triang con los tres vertices iguales
#             dicti[str(label_a)]=np.append(dicti[str(label_a)],[i],axis= 0) #indice de los vertices de cada triangulo
#             dicti_ID_tri[str(label_a)].append(h)          # indices de cada triangulos
#             Lpolygons_etiq.append(str(label_a))          #etiqueta de cada triangulo en lista
#         else:
#             Lpolygons_etiq.append(str('-1'))
#             lista_ID_tri_negros.append(h)
#             continue

#     return dicti_ID_tri, lista_ID_tri_negros


"""formar las nuevas parcelas"""

def new_parcels(new_clusters, polygons,new_labels_arr):
    dicti={}
    dicti_ID_tri={}
    for s in np.unique(new_clusters):
        dicti[str(s)]=np.zeros((0,3))
        dicti_ID_tri[str(s)]=[]
    
        
    
    lista_ID_tri_negros=[]
    Lpolygons_etiq=[]
    for h,i in enumerate(polygons):
        label_a=new_labels_arr[i[0]]
        label_b=new_labels_arr[i[1]]
        label_c=new_labels_arr[i[2]]  
        if label_a==label_b & label_a!=-1:
            dicti[str(label_a)]=np.append(dicti[str(label_a)],[i],axis= 0) #indice de los vertices de cada triangulo
            dicti_ID_tri[str(label_a)].append(h)          # indices de cada triangulos
            Lpolygons_etiq.append(str(label_a))          #etiqueta de cada triangulo en lista
        elif label_b==label_c & label_b!=-1:
            dicti[str(label_b)]=np.append(dicti[str(label_b)],[i],axis= 0)
            dicti_ID_tri[str(label_b)].append(h)
            Lpolygons_etiq.append(str(label_b))
        elif label_a==label_c & label_a!=-1:
            dicti[str(label_a)]=np.append(dicti[str(label_a)],[i],axis= 0)
            dicti_ID_tri[str(label_a)].append(h)
            Lpolygons_etiq.append(str(label_a))
        elif label_a!=label_c & label_a!=label_b & label_b!=label_c:
            if h!=0:            
                lab=Lpolygons_etiq[h-1]
                if lab!='-1':                
                    dicti[lab]=np.append(dicti[lab],[i],axis= 0)
                    dicti_ID_tri[lab].append(h)
                Lpolygons_etiq.append(lab)
                if lab =='-1':
                    lista_ID_tri_negros.append(h)
        else:
            Lpolygons_etiq.append(str('-1'))
            lista_ID_tri_negros.append(h)
            continue
    return dicti_ID_tri, lista_ID_tri_negros


    
"""guardar nuevas parcelas"""
def save_new_parcels(dicti_ID_tri,parcellation_path):
    # triangulos=list(dicti_ID_tri.values())
    
    #Crear carpeta Parcellation
    if not os.path.exists(parcellation_path):
        os.mkdir(parcellation_path)
        print('Carpeta Parcellation creada')
    else: 
        print('La carpeta Parcellation ya existe')
    
    #Crear carpeta segun numero de clusters
    write_path = parcellation_path + f'/{len(dicti_ID_tri)}_clusters'
    if not os.path.exists(write_path):
        os.mkdir(write_path)
        print(f'Carpeta {len(dicti_ID_tri)}_clusters creada')
    else: 
        print(f'La carpeta {len(dicti_ID_tri)}_clusters ya existe')
        
    for m,n in dicti_ID_tri.items():
        BT3.write_parcels(write_path +f'/{m}.parcelsdata', n)
        
        
"""cargar nuevas parcelas"""
def load_parcels (parcels_path): #hard cc final

    parcels = dict()
        
    for Dir in os.listdir(parcels_path):
        parcels[Dir.split('.')[0]] = BT3.read_parcels(parcels_path + "/" + Dir)
            
    return parcels

"""visualizar parcelas independientes"""
def visualize_independent_parcels(mesh_lh_path, L_sp):
    
    import visualizationTools as vt
      
    #Parcelas finales a graficar
    final_parcels = set()

    for k in L_sp.keys():
        final_parcels.add(k)

    fp = list(final_parcels)
    
    import yph_claulib as yc
    # import math  
    
    N_clusters = len(fp)
    paleta=yc.paleta_de_colores(N_clusters)
       
    #Directorios de los mallados corticales
    Lhemi_path = mesh_lh_path

    #Lectura de mallados
    Lvertex, Lpolygons = BT3.read_mesh_obj(Lhemi_path)
    
    Lhemi = vt.Polygon(Lvertex, Lpolygons);
    
    Lhemi.setOpacity(0.1);

    #Creación del render a visualizar
    render = vt.Render();
    
    #Se renderizan los mallados
    render.AddActor(Lhemi);
    
    color_dict={}
    #Para cada parcela del hemisferio ...
    for k, v in L_sp.items():
        #Se selecciona un color de la paleta
        color= paleta[fp.index(k)]
        color_dict[k]=color
        print(color)       
        if len(v) == 0:
            continue
        #Se renderizan los triángulos y polígonos.
        sp_tri = vt.Polygon(Lvertex, Lpolygons[v]);
        sp_tri.setColor((color[0], color[1], color[2]));
        sp_tri.setOpacity(1.0)
        render.AddActor(sp_tri);
       
       
    render.Start();
    del render
    return color_dict

"""visualizar parcelacion basic"""
def visualize_parcellation_basic(mesh_lh_path, L_sp, lista_ID_tri_negros, previous_color_dict={}):
    
    import visualizationTools as vt
   
    
    #Parcelas finales a graficar
    final_parcels = set()

    for k in L_sp.keys():
        final_parcels.add(k)

    fp = list(final_parcels)
    
    import yph_claulib as yc
    # import math  
    
    N_clusters = len(fp)
    paleta=yc.paleta_de_colores(N_clusters*2)
       
    #Directorios de los mallados corticales
    Lhemi_path = mesh_lh_path

    #Lectura de mallados
    Lvertex, Lpolygons = BT3.read_mesh_obj(Lhemi_path)
    
    Lhemi = vt.Polygon(Lvertex, Lpolygons);
    
    Lhemi.setOpacity(0.1);

    #Creación del render a visualizar
    render = vt.Render();
    
    #Se renderizan los mallados
    render.AddActor(Lhemi);
    
    color_repetidos=[]
    for p in list(previous_color_dict.keys()):
        if p in fp:
            color_repetidos.append(previous_color_dict[p])
           
                      
    color_dict={}
    a=0
    #Para cada parcela del hemisferio ...
    for k, v in L_sp.items():
       #Se selecciona un color de la paleta
       if k in list(previous_color_dict.keys()) :
           color = previous_color_dict[k]
       else:
           
           color= paleta[a]
           
           while color in color_repetidos or color in list(color_dict.values()) :
               color=paleta[a+1]
               a+=1
               
               
       color_dict[k]=color
       print(color)       
       if len(v) == 0:
          continue
        #Se renderizan los triángulos y polígonos.
       sp_tri = vt.Polygon(Lvertex, Lpolygons[v]);
       sp_tri.setColor((color[0], color[1], color[2]));
       sp_tri.setOpacity(1.0)
       render.AddActor(sp_tri);
    
        
   
      
    # agregar región negra
    color_negro=(0,0,0)
    negro_tri = vt.Polygon(Lvertex, Lpolygons[lista_ID_tri_negros]);
    negro_tri.setColor(color_negro)
    negro_tri.setOpacity(1.0)
    render.AddActor(negro_tri);
    
       
    render.Start();
    del render
    return color_dict
    
"""visualizar parcelacion herencia"""
def visualize_parcellation_herencia(mesh_lh_path, roi, parcels_dict,previous_color_dict={}):
    
    import bundleTools3 as BT3
    import visualizationTools as vt
   
    #Parcelas a graficar
    final_parcels = set()

    for k in roi:
        final_parcels.add(k)

    fp = list(final_parcels)
    
    #colores
    import yph_claulib as yc
        
    N_clusters = len(fp)
    paleta=yc.paleta_de_colores(N_clusters)
   
    
    #Directorios de los mallados corticales
    Lhemi_path = mesh_lh_path

    #Lectura de mallados
    Lvertex, Lpolygons = BT3.read_mesh_obj(Lhemi_path)
    
    Lhemi = vt.Polygon(Lvertex, Lpolygons);
    
    Lhemi.setOpacity(1);

    #Creación del render a visualizar
    render = vt.Render();
    
    #Se renderizan los mallados
    render.AddActor(Lhemi);
    
    color_dict={}
    
    color_repetidos=[]
    for p in list(previous_color_dict.keys()):
        if p in fp:
            color_repetidos.append(previous_color_dict[p])
        
    
    
    a=0
    for p in fp:        
        triangulos=parcels_dict[str(p)]
        if len(triangulos) == 0:
            continue
        sp_tri = vt.Polygon(Lvertex, Lpolygons[triangulos]);
        if p in list(previous_color_dict.keys()) :
            color = previous_color_dict[p]
        else:
            
            color= paleta[a]
            
            while color in color_repetidos or color in list(color_dict.values()) :
                color=paleta[a+1]
                a+=1
        
        color_dict[p]=color
        
        sp_tri.setColor((color[0], color[1], color[2]));
        sp_tri.setOpacity(1.0)
        render.AddActor(sp_tri);
        
    # agregar región negra
    # color_negro=(0,0,0)
    # negro_tri = vt.Polygon(Lvertex, Lpolygons[lista_ID_tri_negros]);
    # negro_tri.setColor(color_negro)
    # negro_tri.setOpacity(1.0)
    # render.AddActor(negro_tri);
    
    
    render.Start();
    del render
    return color_dict

#%%
# """visualizar parcels level 1..........intento de dividir las parcelas en tonalidades de un mismo color en el subnivel"""
# color_dict={}
# color_dict['60']='red'
# color_dict['63']='green'
# color_dict['64']='brown'
# import visualizationTools as vt
# import vtk_helpers
# from vtk_helpers import *

# def visualize_parcellation_c_memoria(mesh_path, parcels_dict,color_dict, herencia, paleta):
    
   
#     #Parcelas finales a graficar
#     final_parcels = set()
    
#     for k in parcels_dict.keys():
#         final_parcels.add(k)
        
#     fp = list(final_parcels)
    
#     #Directorios de los mallados corticales
#     Lhemi_path = mesh_path

#     #Lectura de mallados
#     Lvertex, Lpolygons = BT3.read_mesh_obj(Lhemi_path)
    
#     Lhemi = vt.Polygon(Lvertex, Lpolygons);
    
#     Lhemi.setOpacity(0.1);

#     #Creación del render a visualizar
#     render = vt.Render();
    
#     #Se renderizan los mallados
#     render.AddActor(Lhemi);
    
#     #Para cada parcela del hemisferio ...
#     # for k, v in parcels_dict.items():
#     #     #Se selecciona un color de la paleta
#     #     for ke, va in herencia:
#     #         if k in va:
#     #             color_base=color_dict[ke]
#     #     colores=        
#     #     color=color_dict
#     #     color= paleta[fp.index(k)]
#     #     color_dict[k]=color
#     #     print(color)       
#     #     if len(v) == 0:
#     #         continue
#     #     #Se renderizan los triángulos y polígonos.
#     #     sp_tri = vt.Polygon(Lvertex, Lpolygons[v]);
#     #     sp_tri.setColor((color[0], color[1], color[2]));
#     #     sp_tri.setOpacity(1.0)
#     #     render.AddActor(sp_tri);
    
    
#     for k,v in herencia.items():
#         color_base=color_dict[str(k)]
#         colores = paleta[color_base]
#         for ind,va in enumerate(v):
#             color=colores[ind]
#             sp_tri = vt.Polygon(Lvertex, Lpolygons[parcels_dict[str(va)]]);
#             sp_tri.setColor((color[0], color[1], color[2]));
#             sp_tri.setOpacity(1.0)
#             render.AddActor(sp_tri);
            
#     # agregar región negra
#     color_negro=(0,0,0)
#     negro_tri = vt.Polygon(Lvertex, Lpolygons[lista_ID_tri_negros_1]);
#     negro_tri.setColor(color_negro)
#     negro_tri.setOpacity(1.0)
#     render.AddActor(negro_tri);
    
#     #agregar centroides
#     # vertex=Lvertex[int(centroides[0])]
#     # vertex=np.reshape(vertex, (1,3))
#     # for c in range(1,len(centroides)):
#     #     v=np.reshape(Lvertex[int(centroides[c])],(1,3))
#     #     vertex=np.append(vertex,v, axis=0)
        
   
#     # point = vt.Points(vertex)
#     # point.setColor((1,1,1))
#     # render.AddActor(point)
    
#     render.Start();
#     del render
#     return 
    
# visualize_parcellation_c_memoria("../data/meshes/001/lh.obj",parcels_dict_1,color_dict, herencia_subnivel_1, paleta)  #seed = 10

#%%
# """visualizar parcelacion con memoria"""...........intento de mostrar las subparcelas en tonalidades de un mismo color
# import visualizationTools as vt
# def visualize_parcellation_c_memoria(mesh_path, L_sp,color_dict_0, lista_ID_tri_negros):  
    
     
#     #Parcelas finales a graficar
#     final_parcels = set()

#     for k in L_sp.keys():
#         final_parcels.add(k)


#     fp = list(final_parcels)
    
    
#     #Directorios de los mallados corticales
#     Lhemi_path = mesh_path

#     #Lectura de mallados
#     Lvertex, Lpolygons = BT3.read_mesh_obj(Lhemi_path)
    
#     Lhemi = vt.Polygon(Lvertex, Lpolygons);
    
#     Lhemi.setOpacity(0.1);

#     #Creación del render a visualizar
#     render = vt.Render();
    
#     #Se renderizan los mallados
#     render.AddActor(Lhemi);

    
#     red=[[1., 0.8352941, 0.8352941],[1., 0.6666667, 0.6666667], [1., 0.5019608, 0.5019608], [1., 0.33333334, 0.33333334],[1. , 0.16470589, 0.16470589], [1., 0., 0.], [0.83137256, 0., 0.], [0.6666667, 0., 0.],[0.5019608, 0., 0.]]
#     green = [[0.01960784, 0.20392157, 0.],[0.02745098, 0.4, 0.],[0.03921569, 0.51372549, 0.],[0.04705882, 0.63921569, 0.],[0.16862745, 0.69411765, 0.10980392],[0.47843137, 0.81568627, 0.38823529],[0.6745098 , 0.90980392, 0.57254902], [0.8 , 0.96470588, 0.68235294]]
    
#     color_dict={}
#     #Para cada parcela del hemisferio ...
#     for k, v in L_sp.items():
#         #Se selecciona un color de la paleta
#         # color= paleta[fp.index(k)]
#         verde=green[fp.index(k)]
#         color_dict[k]=color
#         print(color)       
#         if len(v) == 0:
#             continue
#         #Se renderizan los triángulos y polígonos.
#         sp_tri = vt.Polygon(Lvertex, Lpolygons[v]);
#         # sp_tri.setColor((color[0], color[1], color[2]));
#         sp_tri.setColor((verde[0], verde[1], verde[2]));
#         sp_tri.setOpacity(1.0)
#         render.AddActor(sp_tri);
    
#     # agregar región negra
#     color_negro=(0,0,0)
#     negro_tri = vt.Polygon(Lvertex, Lpolygons[lista_ID_tri_negros]);
#     negro_tri.setColor(color_negro)
#     negro_tri.setOpacity(1.0)
#     render.AddActor(negro_tri);
    
#     #agregar centroides
#     # vertex=Lvertex[int(centroides[0])]
#     # vertex=np.reshape(vertex, (1,3))
#     # for c in range(1,len(centroides)):
#     #     v=np.reshape(Lvertex[int(centroides[c])],(1,3))
#     #     vertex=np.append(vertex,v, axis=0)
        
   
#     # point = vt.Points(vertex)
#     # point.setColor((1,1,1))
#     # render.AddActor(point)
    
#     render.Start();
#     del render
#     return color_dict










