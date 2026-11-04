import intersection_utils as iu
import os

#Paths
meshes_path= 'meshes/'
bundles_path = "set_tractografia/"  
results_path = "results/"

print(os.listdir(meshes_path))
print(os.listdir(bundles_path))
print(os.listdir(results_path))


#Funciones
#=========

#Ejecutar intersección. 
#arg1: Directorio con mallados lh y rh .obj
#arg2: Directorio con archivos .bundles a intersectar.
#arg3: Directorio de resultados
#Return: Escribe datos de intersección en archivo binario .intersectiondata en el siguiente formato: InTri, FnTri, InPoints, FnPoints, fib_index, fiber_class
iu.intersection(meshes_path, bundles_path, results_path)

# #Leer datos de intersección.
# #arg1: Archivo de intersección .intersectiondata
# #Return: Triángulos iniciales, triángulos finales, puntos iniciales, puntos finales, indice de la fibra que intersecta esos triángulos, mallados intersectados (1: LH-LH, 2: RH-RH, 3: LH-RH, 4: RH-LH)
# InTri, FnTri, InPoints, FnPoints, fib_index, fiber_class = iu.read_intersection(results_path + 'centroids.intersectiondata')

# print(InTri.shape)
# # #Seleccionar los índices de cada clase 
# # #arg1: fiber_class obtenido de read_intersection()
# # #Return: Listas que contiene los índices de cada clase.
# fibs_lh, fibs_rh, fibs_inter_l_r, fibs_inter_r_l = iu.return_fiberclass(fiber_class)


# fibers_file = bundles_path + 'centroids.bundles'
# path_classfile = 'results/classfiles/'
# iu.write_fiberclass_files(fibs_lh, fibers_file, fib_index, path_classfile, 'fibs_lh')
# iu.write_fiberclass_files(fibs_rh, fibers_file, fib_index, path_classfile, 'fibs_rh')    
# iu.write_fiberclass_files(fibs_inter_l_r, fibers_file, fib_index, path_classfile, 'fibs_l_r')
# iu.write_fiberclass_files(fibs_inter_r_l, fibers_file, fib_index, path_classfile, 'fibs_r_l')
# parcels_gii_lh = 'parcels_gii/lh.r.aparc.annot.gii'
# parcels_gii_rh = 'parcels_gii/rh.r.aparc.annot.gii' 
# iu.connectome_gifti(results_path + 'centroids.intersectiondata', meshes_path, parcels_gii_lh, parcels_gii_rh, results_path, 'conmatrix')


#IN PROGRESS...
#parcels_pd_lh = 'parcels_parcelsdata/LH/'
#parcels_pd_rh = 'parcels_parcelsdata/RH/'
#iu.load_parcels(parcels_pd_lh)
#iu.connectome_parcelsdata(results_path + 'tractography-streamline-regularized-deterministic_001_21p.intersectiondata', meshes_path, parcels_pd_lh, parcels_pd_rh, results_path, 'conmatrix_pd')

