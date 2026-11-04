Fiber Filtering
======================
Este código permite filtrar las fibras de un sujeto en base a su largo en mm, eliminando aquellas que no superan un umbral definido.
El largo de una fibra se aproxima como:

<img src="https://latex.codecogs.com/svg.latex?L(f)=(N_p-1)\cdot(p_2&space;-&space;p_1)" title="L(f)=(N_p-1)\cdot(p_2 - p_1)" />

donde <img src="https://latex.codecogs.com/svg.latex?p_1,&space;p_2" title="p_1, p_2" /> son el primer y segundo punto de la fibra respectivamente, y <img src="https://latex.codecogs.com/svg.latex?N_p" title="N_p" /> es el número de puntos (en este caso <img src="https://latex.codecogs.com/svg.latex?N_p=21" title="N_p=21" />).

Cabe mencionar que este código es rudimentario y no es necesariamente la mejor implementación.

## Dependencias de Código
Para utilizar el código, es necesario instalar las siguientes librerías:
- Numpy: https://numpy.org/

### OPCIÓN 1: Instalación de dependencias via pip3 en Ubuntu
```
pip3 install numpy
```

### OPCIÓN 2: Instalación de dependencias via apt en Ubuntu
```
sudo apt install python3-numpy
```

## Datos de ejemplo
En el siguiente enlace hay un sujeto de la base de datos ARCHI, con sus fibras remuestreadas a 21 puntos por fibra.
https://drive.google.com/drive/folders/1-qYE4iCXVQHoxkwSgcqW1V2ExnsSvk4D?usp=sharing

## Sintaxis
```
python filtering.py input_bundle.bundles output_bundle.bundles min_size_in_mm
```

- **input_bundle**: Archivo de tractografía de entrada (debe estar remuestreada a 21 puntos) en formato .bundles.
- **output_bundle**: Archivo de tractografía filtrada de salida, en formato .bundles y .bundlesdata.
- **min_size_in_mm**: Tamaño mínimo de fibra en milímetros; cualquier fibra más corta será eliminada.

## Ejemplo de uso
```
python filtering.py example_sub_resampled.bundles example_sub_resampled_filtered.bundles 40
```

Nota: A pesar de que la llamada al código sólo pida el archivo .bundles, el archivo .bundlesdata correspondiente debe estar presente en la misma carpeta que el archivo .bundles.
