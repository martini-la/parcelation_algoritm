import os

def count_fibers_in_bundle(file_path):
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

def main():

    # Obtenemos el path de la carpeta target
    path = os.getcwd()
    filtered_tractograms = os.path.join(path, "intento_3", "intrasujeto","filtered_tractograms")
    subs = os.listdir(filtered_tractograms)

    # Contamos las fibras abriendo el .bundle en un diccionario
    fiber_counts = {}
    for sub in subs:
        bundle_file = os.path.join(filtered_tractograms, sub, "20mmfiltered_MNI_21p.bundles")      
        num_fibers = int(count_fibers_in_bundle(bundle_file).rstrip(','))
        fiber_counts[sub] = num_fibers

    # Mostrar los resultados
    suma = 0
    num_subjects = len(fiber_counts)
    for sub, count in fiber_counts.items():
        print(f"Sujeto: {sub} - Número de fibras: {count}")
        suma = suma + int(count)

    # Calcular y mostrar el promedio
    promedio = suma / num_subjects
    print(f"\n\tPromedio del número de fibras por sujeto tras filtrado: {promedio}")
    print(f"\n\tfiber_counts.values: {fiber_counts.values()}")

    # Calcular la varianza
    suma_cuadrados_diferencias = sum((count - promedio) ** 2 for _, count in fiber_counts.items())
    varianza = suma_cuadrados_diferencias / (num_subjects - 1)
    print(f"\n\tVarianza del número de fibras por sujeto tras filtrado: {varianza}")

if __name__ == "__main__":
    main()
