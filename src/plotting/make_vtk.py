import sys

import vtk as vtk
import numpy as np

def save_vtk_unstructured_grid(path, points, cells, point_data):
    """
    Modified from: https://stackoverflow.com/questions/65683933/write-a-vtk-file-from-adaptive-mesh-refinement-with-python

    Create a vtk grid file containing quads from a list of points,
    a list of cells and additional point data.
    The list of cells references the points inside the point list via the row index.

    N Points: [[x_0, y_0, z_0],
               [x_1, y_1, z_1],
               ...
               [x_(n-1), y_(n-1), z_(n-1)]]

    M Cells: [[i_00, i_01, i_02, i_03],
              [i_10, i_11, i_12, i_13],
              ...
              [i_(m-1)0, i_(m-1)1, i_(m-1)2, i_(m-1)3]]

    E.g.:
    Cell: p0 x------x p1    =>      Cell indices inside the cell array:
             |      |               [0, 1, 2, 3]
             |      |
             |      |
          p2 x------x p3

    :param path: Save path as string
    :param points: Nx3 numpy array of point coordinates
    :param cells: Mx4 numpy array of point indices for a mesh consisting of quads.
    :param point_data: Nx1 numpy array of containing data at the point coordinates.
    """
    vtk_points = vtk.vtkPoints()
    vtk_cells = vtk.vtkCellArray()

    # insert points
    for p in points:
        vtk_points.InsertNextPoint(p[0], p[1], 0)

    # insert the quad cells
    for idx in cells:
        # create a new quad cell
        quad = vtk.vtkQuad()
        quad.GetPointIds().SetId(0, int(idx[0]))
        quad.GetPointIds().SetId(1, int(idx[1]))
        quad.GetPointIds().SetId(2, int(idx[2]))
        quad.GetPointIds().SetId(3, int(idx[3]))
        vtk_cells.InsertNextCell(quad)

    # # create the point data
    # data = vtk.vtkDoubleArray()
    # data.SetNumberOfComponents(1)
    # data.SetNumberOfValues(len(point_data))
    # for i, d in enumerate(point_data):
    #     print(i, d)
    #     data.SetTuple(i, d)

    # create the grid from the points and cells
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(vtk_points)
    grid.SetCells(vtk.VTK_QUAD, vtk_cells)
    #grid.GetPointData().SetScalars(data)

    # write the grid
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(path)
    writer.SetInputData(grid)
    writer.Write()

def main():
    path = "boundary.vtu"
    points_path = sys.argv[1]
    points = np.loadtxt(points_path)
    cell_path = sys.argv[2]
    cells = np.loadtxt(cell_path)
    point_data = np.zeros(np.size(points))
    save_vtk_unstructured_grid(path, points, cells, point_data)

if __name__ == "__main__":
    main()