import pyvista as pv
import numpy as np

def create_quad4_square(side_length=1.0, n=4, z_translation=0.0):
    """
    Create a PyVista mesh of a square made of quad4 elements.

    Parameters
    ----------
    side_length : float, optional
        Length of the square sides (default is 1.0).
    n : int, optional
        Number of quad elements along each side (default is 4).
    z_translation : float, optional
        Translation of the square in the z-direction (default is 0.0).

    Returns
    -------
    mesh : pyvista.UnstructuredGrid
        A PyVista mesh representing a structured square grid of quad4 elements.
    """
    # Generate grid points
    x = np.linspace(-side_length/2, side_length/2, n+1)
    y = np.linspace(-side_length/2, side_length/2, n+1)
    xx, yy = np.meshgrid(x, y)
    zz = np.full_like(xx, z_translation)

    # Flatten to a (N, 3) array of point coordinates
    points = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]

    # Build quad connectivity
    cells = []
    for j in range(n):
        for i in range(n):
            # node indexing pattern
            p0 = j*(n+1) + i
            p1 = p0 + 1
            p2 = p0 + (n+1) + 1
            p3 = p0 + (n+1)
            cells.extend([4, p0, p1, p2, p3])  # 4 nodes per quad

    # Convert to numpy arrays
    cells = np.array(cells)
    celltypes = np.full(n*n, 9, dtype=np.uint8)  # 9 = VTK_QUAD

    # Create mesh
    mesh = pv.UnstructuredGrid(cells, celltypes, points)

    return mesh


# Example usage
if __name__ == "__main__":
    mesh = create_quad4_square(side_length=66.0, n=100, z_translation=122)
    plotter = pv.Plotter()
    plotter.add_mesh(mesh, show_edges=True, color='lightblue')
    plotter.show()
