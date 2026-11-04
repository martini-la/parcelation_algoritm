Ejecutar archivo "trm.py" para aplicar una matriz de transformación afín al archivo de tractografía.

Uso (desde Terminal):

python trm.py fibras_input.bundles fibras_output.bundles matriz_transformacion.trm

- fibras_input.bundles: Archivo tractografía de entrada, a la cual se le aplicará la matriz.
- fibras_output.bundles: Archivo en el cual se guardará la tractografía transformada. 
- matriz_transformacion.trm: Archivo que contiene la matriz de transformación afín a aplicar a la tractografía.

Nota 1: No quitar los archivos "bundleTools.py" y "bundleTools.pyc".
Nota 2: El código está escrito y probado en Python 2.