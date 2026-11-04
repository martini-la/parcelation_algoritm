import os
import subprocess

def run_segmentation_clust():
    try:
        # Ejecutamos el código de segmentación
        result = subprocess.run(['gcc', '-fPIC', '-shared', '-O3', '-o', 'FFClust-master/segmentation_clust_v1.2/segmentation.so',
                                'FFClust-master/segmentation_clust_v1.2/segmentation.c', '-fopenmp', '-ffast-math'],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Verificar el código de salida
        if result.returncode == 0:
            print("El algoritmo de segmentación se ejecutó correctamente.")
        else:
            print("Hubo un problema al ejecutar el algoritmo de segmentación.")
            print("Error:", result.stderr.decode('utf-8'))

    except subprocess.CalledProcessError as e:
        print("Error al ejecutar el comando:")
        print(e)
        print("Output del comando (stdout):", e.stdout.decode('utf-8'))
        print("Error del comando (stderr):", e.stderr.decode('utf-8'))

def main():
    """
    Aplicación preliminar de ffclust
    """
    path = os.getcwd()
    filter_thresh = 60
    input_file = os.path.join(path, "112516", f"{filter_thresh}mmfiltered_MNI_21p.bundles")

    output_dir = "results"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    run_segmentation_clust()

    parameters_ffclust = {
        "points": "0,3,10,17,20",
        "ks": "150,125,125,125,150",
        "thr-seg": 6,
        "thr-join": 6,
        "filter_thresh": filter_thresh
    }

    # Construir el comando con todos los parámetros
    sub = "112516"
    command = [
        "python3", "FFClust-master/main.py",
        "--infile", input_file,
        "--outdir", output_dir,
        "--ks"
    ] + parameters_ffclust["ks"].split(',') 

    try:
        subprocess.run(command, stdout=None, stderr=None, text=True, check=True)
        print(f"Clustering exitoso para {sub}.")
    except subprocess.CalledProcessError as e:
        print(f"Error en la ejecución del comando para {sub}: {e}")
        print(f"Código de retorno: {e.returncode}")

if __name__ == "__main__":
    main()
