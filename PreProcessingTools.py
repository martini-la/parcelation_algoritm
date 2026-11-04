from BundleTools import bundleTools as bt
from BundleTools import bundleTools3 as bt3
import os
import time
import subprocess
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from sklearn.cluster import MiniBatchKMeans

"""
Filtrado de tractogramas dado un umbral
"""
def filter_allsubs(input_dir, output_dir, min_size):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    subs = os.listdir(input_dir)
    avance = 0
    for sub in subs:
        input_bundle = os.path.join(input_dir, sub, "3Msift_t_MNI_21p.bundles")

        output_sub_dir = os.path.join(output_dir, sub)
        if not os.path.exists(output_sub_dir):
            os.makedirs(output_sub_dir)

        output_bundle = os.path.join(output_sub_dir, str(min_size) + "mmfiltered_MNI_21p.bundles")
        command = ["python", "filtering/filtering.py", input_bundle, output_bundle, str(min_size)]
        subprocess.run(command)
        avance += 1
        print('Tractografia filtrada: ', sub, "\tAvance (%): ", avance / len(subs) * 100)

"""
Conteo de fibras en archivos .bundles
"""

class quantify_tractograms():
    def __init__(self, bundles_dir, bundle_name, path):

        self.path = path
        subs = os.listdir(bundles_dir)

        # Contamos las fibras abriendo el .bundle en un diccionario
        fiber_counts = {}
        for sub in subs:
            bundle_file = os.path.join(bundles_dir, sub, bundle_name)      
            num_fibers = int(self.count_fibers_in_bundle(bundle_file).rstrip(','))
            fiber_counts[sub] = num_fibers

        # Mostrar los resultados
        suma = 0
        num_subjects = len(fiber_counts)
        for sub, count in fiber_counts.items():
            print(f"Sujeto: {sub} - Número de fibras: {count}")
            suma += int(count)

        # Calcular y mostrar el promedio
        promedio = suma / num_subjects
        print(f"\n\tPromedio del número de fibras por sujeto tras filtrado: {promedio}")
        print(f"\n\tfiber_counts.values: {fiber_counts.values()}")

        # Calcular la varianza
        suma_cuadrados_diferencias = sum((count - promedio) ** 2 for count in fiber_counts.values())
        varianza = suma_cuadrados_diferencias / (num_subjects - 1)
        desviacion_estandar = math.sqrt(varianza)
        print(f"\n\tVarianza del número de fibras por sujeto tras filtrado: {varianza}")
        print(f"\n\tDesviación estandar del número de fibras por sujeto tras filtrado: {desviacion_estandar}")
        
        # Graficar la distribución de los datos
        self.plot_distribution(fiber_counts, save_dir=os.path.dirname(os.path.dirname(bundles_dir)))

    def count_fibers_in_bundle(self, file_path):
        """Función para contar el número de fibras en un archivo .bundle."""
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if 'curves_count' in line:
                        # Suponiendo que la línea tiene el formato 'curves_count: valor'
                        _, value = line.split(':')
                        return value.strip()
        except Exception as e:
            print(f"Error al leer el archivo {file_path}: {e}")
            return None

    def plot_distribution(self, fiber_counts, save_dir):
        """Función para graficar la distribución de los datos."""
        counts = list(fiber_counts.values())

        plt.figure(figsize=(12, 6))

        # Histograma
        plt.subplot(1, 2, 1)
        plt.hist(counts, bins=30, edgecolor='black')
        plt.title('Distribución del número de fibras')
        plt.xlabel('Número de fibras')
        plt.ylabel('Frecuencia')

        # Gráfico de caja
        plt.subplot(1, 2, 2)
        plt.boxplot(counts, vert=False)
        plt.title('Gráfico de caja del número de fibras')
        plt.xlabel('Número de fibras')

        plt.tight_layout()
        
        # Guardar el gráfico como archivo de imagen
        plt.savefig(os.path.join(self.path, save_dir, "filtered_fibers_quantity_distribution.png"))
        print("El gráfico de distribución se ha guardado como filtered_fibers_quantity_distribution.png")



"""
Elbow method
"""

# Función para calcular la suma de cuadrados intra-cluster
def calculate_inertia(points, clusters):
    kmeans = MiniBatchKMeans(n_clusters=clusters, random_state=42)
    kmeans.fit(points)
    return kmeans.inertia_

def elbow_method(path, num_intento, filter_thresh, subs_limit = 100):
    filtered_tractograms = os.path.join(path,"intento_" + str(num_intento),"intrasujeto","filtered_tractograms")

    # Ordenamos los sujetos y escojemos primeros 5. Luego establecemos los puntos de interés. 
    subs = sorted(os.listdir(filtered_tractograms))
    points_of_interest = [0,4,10,16,20]
    cluster_numbers = [25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600]

    count_subs = 0
    for sub in subs:
        # Importamos el archivo de tractografía 
        print("Trabajando sujeto ",sub)
        bundle_file = os.path.join(filtered_tractograms, sub, str(filter_thresh) + "mmfiltered_MNI_21p.bundles")
        try:
            tractogram = bt.read_bundle(bundle_file)
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
        output_dir = os.path.join(path, "intento_" + str(num_intento), "intrasujeto", "elbow_method")
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
        count_subs += 1
        if count_subs > subs_limit:
            break

"""
Aplicación de ffclust
"""

def apply_ffclust_allsubs(input_dir, output_dir, parameters, blacklist = []):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    subs = os.listdir(input_dir)
    avance = 0

    # Compilar el código C con gcc
    subprocess.run(['gcc', '-fPIC', '-shared', '-O3', '-o', 'FFClust-master/segmentation_clust_v1.2/segmentation.so',
                    'FFClust-master/segmentation_clust_v1.2/segmentation.c', '-fopenmp', '-ffast-math'])

    for sub in subs:
        print("trabajando sujeto ", sub)
        if sub in blacklist:
            print(f"Sujeto {sub} está en la lista negra. Saltando...")
            continue
        
        input_bundle = os.path.join(input_dir, sub, str(parameters["filter_thresh"])+"mmfiltered_MNI_21p.bundles")
        output_sub_dir = os.path.join(output_dir, sub)
        if not os.path.exists(output_sub_dir):
            os.makedirs(output_sub_dir)

        # Construir el comando con todos los parámetros
        command = [
            "python3", "FFClust-master/main.py",
            "--infile", input_bundle,
            "--outdir", output_sub_dir,
            "--ks"
        ] + parameters["ks"].split(',') 
        
        try:
            subprocess.run(command, stdout=None, stderr=None, text=True, check=True)
            print(f"Clustering exitoso para {sub}.")
        except subprocess.CalledProcessError as e:
            print(f"Error en la ejecución del comando para {sub}: {e}")
            print(f"Código de retorno: {e.returncode}")

        avance += 1
        print(f'Avance (%): {avance / len(subs) * 100}', "\tSujeto: ", sub)


