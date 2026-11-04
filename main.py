from BundleTools import bundleTools as bt
from BundleTools import bundleTools3 as bt3
import os
from PreProcessingTools import filter_allsubs, apply_ffclust_allsubs, quantify_tractograms, elbow_method
# from InterClusterTools import save_all_centroids, load_filtered_centroids, interclustering

def main():

    path = os.getcwd()
    original_subs = "/home/martin/tesis/original_tractograms"
    num_intento = 3

    """
    Aplicación del filtrado de sujetos
    """
    # Filtrado de sujetos
    filter_thresh = 60
    filtered_allsubs_output_dir = os.path.join(path, "intento_" + str(num_intento), "intrasujeto", "filtered_tractograms")
    # filter_allsubs(original_subs, filtered_allsubs_output_dir, filter_thresh)

    """
    Cuantificación de fibras del filtrado
    """
    # quantify_tractograms(filtered_allsubs_output_dir, "60mmfiltered_MNI_21p.bundles", path)

    """
    Aplicación del elbow method para encontrar puntos preliminares óptimos de ffclust.
    """
    # elbow_method(path, num_intento, filter_thresh, subs_limit = 3)
    """
    Aplicación del clustering intrasujeto usando FFClust
    """

    filtered_tractograms_dir = os.path.join(path,"intento_"+ str(num_intento), "intrasujeto", "filtered_tractograms")
    ffclust_clusters_dir = os.path.join(path,"intento_"+ str(num_intento), "intrasujeto", "ffclust_clusters")

    parameters_ffclust = {
        "points": "0,3,10,17,20",
        "ks": "150,125,125,125,150",
        "thr-seg": 6,
        "thr-join": 6,
        "filter_thresh": filter_thresh
    }
    # blacklist = []
    apply_ffclust_allsubs(filtered_tractograms_dir, ffclust_clusters_dir, parameters_ffclust)

    """
    Guardamos todos los centroides de los clusters 
    """

    allsubsdata_dir = os.path.join(path,"intento_"+str(num_intento), "intersujeto", "allsubsdata")
    # Si no se desea filtrar los clusters, dejar todos los parámetros como False.
    parameters_save_all_centroids = {
        "remove_little_ones": 10, # tamaño mínimo que debe tener un cluster (10 fibras)
        "size": 1.5, # Fibras de menos de size [mm] no se consideran.
        "nearest": 40 # 
    }
    # save_all_centroids(ffclust_clusters_dir, allsubsdata_dir, parameters_save_all_centroids)

    """
    Lodeamos los centroides del archivo numpy que almacenó todos los clusters, con el fin de aplicar
    quickbundles como clustering intersujeto. 
    """
    
    # allsubs_centroids_file = os.path.join(allsubsdata_dir, "all_centroids.npy")    
    # allsubs_centroids = load_filtered_centroids(allsubs_centroids_file)

    """
    Aplicamos clustering intersujeto. 
    """
    
    # theta_QB = 21
    # QB_output_dir = os.path.join(path, "intersujeto", "QB_output")
    # interclustering(theta_QB, allsubs_centroids, QB_output_dir)
    


if __name__ == "__main__":
    main()

