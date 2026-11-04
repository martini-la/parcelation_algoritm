import numpy as np
import vtk
print(vtk.__version__)
import os
import bundleTools3 as bt3

def read_bundle_severalbundles(file_path):
    with open(file_path, 'rb') as f:
        bundle_names = []
        fibers = []
        
        while True:
            # Leer el nombre del bundle
            name_length = int.from_bytes(f.read(4), 'little')
            if not name_length:
                break
            name = f.read(name_length).decode('utf-8')
            bundle_names.append(name)
            
            # Leer el número de puntos en el bundle
            p = int.from_bytes(f.read(4), 'little')
            print(f"Bundle name: {name}, number of points: {p}")  # Mensaje de depuración
            
            if p < 0:
                raise ValueError(f"Invalid number of points: {p}")
            
            if p == 0:
                continue
            
            # Leer los vértices
            vertex = np.frombuffer(f.read(p * 3 * 4), 'f').reshape(-1, 3)
            fibers.append(vertex)
        
        return fibers, bundle_names

# Función para convertir fibras a vtkPolyData
def fibers_to_vtk_polydata(fibers):
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()

    point_id = 0
    for fiber in fibers:
        line = vtk.vtkPolyLine()
        line.GetPointIds().SetNumberOfIds(len(fiber))
        for i, point in enumerate(fiber):
            points.InsertNextPoint(point)
            line.GetPointIds().SetId(i, point_id)
            point_id += 1
        lines.InsertNextCell(line)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    return polydata

# Función para visualizar vtkPolyData
def visualize_polydata(polydata):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.1, 0.2, 0.3)  # Fondo azul oscuro

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    render_window.Render()
    interactor.Start()

path = os.getcwd()

bundle_file = os.path.join(path, "visualization_test", "100206", "centroids.bundles")
hie_file = os.path.join(path, "visualization_test", "QB_21.hie")
fibers, bundle_names = read_bundle_severalbundles(bundle_file)

# Convertir las fibras a vtkPolyData
# Aquí asumimos que quieres visualizar todas las fibras
all_fibers = [fiber for bundle in fibers for fiber in bundle]
polydata = fibers_to_vtk_polydata(all_fibers)

# Visualizar las fibras
visualize_polydata(polydata)