import numpy as np
import os
from bundleTools import read_bundle
from sklearn.cluster import MiniBatchKMeans
import matplotlib.pyplot as plt

# Función para calcular la suma de cuadrados intra-cluster
def calculate_inertia(points, clusters):
    kmeans = MiniBatchKMeans(n_clusters=clusters, random_state=42)
    kmeans.fit(points)
    return kmeans.inertia_

# Establecemos la ruta del directorio que contiene las tractografías originales
path = os.getcwd()
original_tractograms = os.path.join(path,"intento_2","intrasujeto","filtered_tractograms")

# Ordenamos los sujetos y escojemos primeros 5. Luego establecemos los puntos de interés. 
subs = sorted(os.listdir(original_tractograms))
points_of_interest = [0,4,10,16,20]
cluster_numbers = [25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600]

for sub in subs:
    # Importamos el archivo de tractografía 
    print("Trabajando sujeto ",sub)
    bundle_file = os.path.join(original_tractograms, sub,"40mmfiltered_MNI_21p.bundles")
    try:
        tractogram = read_bundle(bundle_file)
        tractogram_len = len(tractogram)
    except Exception as e:
        print(f"Error al leer el archivo {bundle_file}: {e}")
        continue

    # Diccionario para almacenar los puntos de interés
    points_dict = {point_index: [] for point_index in points_of_interest}

    print("Guardando puntos de interés.")
    for fiber in tractogram: 
        # Iteramos entre los valores de points of interest, los cuales pueden ser 
        # [0,4,10,16,20]
        for point in points_of_interest:
            # Almacenamos los puntos de interés en el diccionario
            points_dict[point].append(fiber[point])

    # Convertir listas en arrays numpy para facilitar el cálculo    
    for point_index in points_of_interest:
        points_dict[point_index] = np.array(points_dict[point_index])
    print("Puntos de interés guardados...")

    # Graficamos el método del codo y guardamos la gráfica como archivo PNG
    print("Iniciando creación de plots.")
    plt.figure()
    # Calculamos la inercia para los números de clusters especificados para cada punto de interés
    for point_index in points_of_interest:
        print("point: ", point_index)
        inertias = []
        for k in cluster_numbers:
            inertia = calculate_inertia(points_dict[point_index], k)
            inertias.append(inertia)
            print("\tn_clusters: ",k,"\tinertia: ", inertia)
        plt.plot(cluster_numbers, inertias, label=f'Punto de interés {point_index}')
    
    # Crear el directorio si no existe
    output_dir = os.path.join(path, "intento_2", "intrasujeto", "elbow_method")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Configurar las propiedades del gráfico
    plt.xlabel('Número de Clusters')
    plt.ylabel('Inercia')
    plt.title(f'Método del Codo - Sujeto {sub}')
    plt.xticks(cluster_numbers, rotation=90)
    plt.legend()
    plt.grid(True)

    # Guardar y mostrar el gráfico
    output_file = os.path.join(output_dir, f"elbow_method_{sub}.png")
    plt.savefig(output_file)
    plt.close()  # Cerrar la figura para liberar memoria
    print(f"Gráfica guardada como {output_file}")
        




