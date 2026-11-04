import os
import bundleTools3 as bt3
import time
import numpy as np
import random 
import nibabel as nib

from dipy.segment.metric import AveragePointwiseEuclideanMetric
from dipy.segment.clustering import QuickBundles
from dipy.io.pickles import save_pickle

"""
Guardamos en un único archivo los centroides y los clusters de todos los sujetos, con el fin 
de realizar el resto de las operaciones inter cluster. 
"""

class save_all_centroids:
    def __init__(self, input_dir, output_dir, parameters):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.parameters = parameters
        self.save_data()

    def save_data(self):

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        subs = os.listdir(self.input_dir)
        avance = 0

        all_centroids = []    #Arreglo que contiene todos los centroides
        all_data = []         #Arreglo que contiene todas las fibras
        ids = []              #Arreglo que contiene todos los identificadores (a qué sujeto pertenece cada cluster)
        
        t0 = time.time()

        for sub in subs: 

            # Lectura de clusters por sujeto
            clusters = np.array(bt3.read_bundle_severalbundles(os.path.join(self.input_dir, sub, "finalClusters.bundles"))[0], dtype=object)
            
            # Filtramos los clusters más pequeños.
            centroids, data = self.get_centroids(clusters, sub, self.parameters["remove_little_ones"], self.parameters["size"], self.parameters["nearest"])

            # Almacenamos la data en la lista y creamos un identificador que se repetirá tantas veces como 
            # largo tenga el cluster filtrado. De esta manera se podrá indexar sobre las fibras almacenadas. 

            
            all_centroids.append(centroids)
            all_data.append(data)
            ids = [sub]*len(data)
            avance = avance + 1
            print("avance: ", avance / len(subs) * 100, "\tProcesando sujeto ", sub)

        print('Tiempo Total filtrado de centroides: ' + str(time.time() - t0) + ' [s]')    
        print('Saving all_data...')

        np.save(os.path.join(self.output_dir, "all_centroids"), np.concatenate(all_centroids))
        np.save(os.path.join(self.output_dir, "all_data"), np.concatenate(all_data))
        np.save(os.path.join(self.output_dir, "ids"), np.array(ids))
        print('Tiempo Total guardado de centroides: ' + str(time.time() - t0) + ' [s]')
        
        #os._exit(00) #Restart Kernel for freeing memory
        
        print('Saved...')

    def get_centroids(self, clusters, sub, remove_little_ones = False, size=False, nearest=False):
        centroids = []
        original_data = []

        og_centroids = np.array(bt3.read_bundle_severalbundles(os.path.join(self.input_dir, sub, "centroids.bundles"))[0], dtype = object)
    
        if(remove_little_ones):
            removed1 = 0
            removed2 = 0
            for i,cluster in enumerate(clusters):
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
            print(str(removed1)+' and '+str(removed2)+' removed from '+str(len(clusters))+'. '+str(len(centroids))+' remaining centroids.'+str(len(np.concatenate(clusters)))+' remaining fibers')
        else:
            for cluster in clusters:
                centroids.append(np.mean(cluster, axis = 0))
                original_data.append(np.array(cluster))
        return np.array(centroids, dtype=object), np.array(original_data, dtype = object)
    

"""
Lodea los centroides del archivo con todos los sujetos. 
"""

#Carga de centroides filtrados
def load_filtered_centroids (centroids_file):
    #Carga de identificadores => elemento i indica a qué sujeto pertenece el centroide i.
    # ids = np.load(centroids_path + 'ids.npy')
    
    #Carga de centroides intra-sujeto y conversión a ArraySequence (AS, formato utilizado por DIPY/QuickBundles)
    intra_centroids = np.load(centroids_file)
    intra_centroids_AS = nib.streamlines.ArraySequence(intra_centroids)
    
    # return ids, intra_centroids_AS
    return intra_centroids_AS

"""
Aplica clustering intersujeto. 
"""

#Ejecucion de algoritmo de interclustering QuickBundles
def interclustering (theta_threshold, centroids_data, output_dir):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #Como las fibras ya están muestreadas con 21 puntos, se evita el resampleo interno de QB...
    metric = AveragePointwiseEuclideanMetric()
    
    #Umbral de distancia de QB en mm (notar que debe ser un float)
    theta = float(theta_threshold)

    #Creación de método QB
    qb = QuickBundles(threshold = theta, metric = metric)
    
    #Clustering
    t0 = time.time()
    inter_clusters = qb.cluster(centroids_data)
    
    print('Tiempo total interclustering: ' + str(time.time() - t0) + ' [s]')
    
    #Guardado de clustering en pkl
    save_pickle(output_dir + 'centroids_QB_' + str(theta_threshold) + 'mm.pkl', inter_clusters)
    
    return inter_clusters

"""
Realizamos un guardado de los clusters y centroides intersujeto filtrados según la cantidad de sujetos que componen dicho cluster.
Consideramos como clusters representativos aquellos que contienen centroides de al menos un 88% de la población. 
Los clusters restantes serán descartados debido a que no se consideran representativos de la población de sujetos. 
"""

