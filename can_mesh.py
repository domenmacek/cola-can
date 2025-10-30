import numpy as np
import pyvista as pv
import pygmsh

# ===========================================
# PARAMETERS
# ===========================================
R_body = 33.0
R_cap = 27.0
height_total = 122
n_theta = 192
n_body = 100
n_loft = 10
loft_height = 9
cap_slope = 0.00
cap_mesh_size = 5.  # target triangle size inside cap (smaller -> finer)
VTK_QUAD = 9
VTK_TRIANGLE = 5

# ===========================================
# BODY + LOFT (unchanged)
# ===========================================
pts = []
cells = []
cell_types = []

z_top = height_total / 2
z_bottom = -height_total / 2
z_top_loft_start = z_top - loft_height
z_bottom_loft_end = z_bottom + loft_height
theta_vals = np.linspace(0, 2*np.pi, n_theta, endpoint=False)

# body points
z_vals = np.linspace(z_bottom_loft_end, z_top_loft_start, n_body + 1)
for z in z_vals:
    for th in theta_vals:
        pts.append([R_body*np.cos(th), R_body*np.sin(th), z])

def idx_body(j, i): return j * n_theta + (i % n_theta)
for j in range(n_body):
    for i in range(n_theta):
        a = idx_body(j,i)
        b = idx_body(j,i+1)
        c = idx_body(j+1,i+1)
        d = idx_body(j+1,i)
        cells.extend([4,a,b,c,d])
        cell_types.append(VTK_QUAD)

# loft
def add_loft(z_start, z_end, R_start, R_end, upward=True):
    base = len(pts)
    for j in range(n_loft + 1):
        t = j / n_loft
        z = (1 - t) * z_start + t * z_end
        R = (1 - t) * R_start + t * R_end
        for th in theta_vals:
            pts.append([R*np.cos(th), R*np.sin(th), z])
    for j in range(n_loft):
        for i in range(n_theta):
            a = base + j*n_theta + i
            b = base + j*n_theta + (i+1)%n_theta
            c = base + (j+1)*n_theta + (i+1)%n_theta
            d = base + (j+1)*n_theta + i
            order = [a,b,c,d] if upward else [d,c,b,a]
            cells.extend([4]+order)
            cell_types.append(VTK_QUAD)
    return base + n_loft*n_theta  # index of outer ring

top_ring_start = add_loft(z_top_loft_start, z_top, R_body, R_cap, upward=True)
bottom_ring_start = add_loft(z_bottom_loft_end, z_bottom, R_body, R_cap, upward=False)

# ===========================================
# UNSTRUCTURED TRIANGULAR CAPS USING PYGMSH (robust)
# ===========================================
def generate_cap_with_gmsh(boundary_pts, mesh_size=0.05):
    """
    Generate an unstructured triangular mesh inside the polygon defined by boundary_pts.
    boundary_pts : (N,3) array of 3D coordinates lying on the cap plane (z constant).
    Returns: mesh_points (M,3), triangle_cells (K,3) (indices into mesh_points)
    """
    # polygon points must be 3D; pygmsh.add_polygon accepts list-of-lists
    polygon_pts = boundary_pts.tolist()

    # Use pygmsh.geo.Geometry() context; add_polygon is supported across versions
    with pygmsh.geo.Geometry() as geom:
        # add_polygon will create the boundary loop and a plane surface inside
        poly = geom.add_polygon(polygon_pts, mesh_size=mesh_size)
        mesh = geom.generate_mesh()
    points = mesh.points  # (M,3)
    # mesh.cells is dict-like in newer meshio; use get_cells_type
    try:
        tri_cells = mesh.get_cells_type("triangle")
    except Exception:
        # backward compatible fallback
        # mesh.cells might be a list of (type, data)
        tri_cells = None
        for cell_block in mesh.cells:
            if cell_block.type == "triangle":
                tri_cells = cell_block.data
                break
        if tri_cells is None:
            raise RuntimeError("No triangle cells returned by Gmsh mesh")
    return points, tri_cells

def add_unstructured_cap(pts, cells, cell_types, ring_start, n_theta, z, mesh_size=0.05):
    """
    Create an unstructured triangular cap whose boundary nodes are exactly
    the loft outer ring nodes (no duplicates).
    """
    import numpy as np
    ring_pts = np.array([pts[ring_start + i] for i in range(n_theta)])
    ring_pts[:, 2] = z  # flatten exactly on z-plane

    # === Generate gmsh mesh ===
    gmsh_pts, gmsh_tri = generate_cap_with_gmsh(ring_pts, mesh_size=mesh_size)
    gmsh_pts = np.array(gmsh_pts)

    # === Identify boundary nodes (the ones on the polygon edge) ===
    # We'll find which Gmsh nodes are nearly on the boundary polygon by radius.
    # Compute distance from center:
    r_cap = np.mean(np.sqrt(ring_pts[:, 0]**2 + ring_pts[:, 1]**2))
    gmsh_r = np.sqrt(gmsh_pts[:, 0]**2 + gmsh_pts[:, 1]**2)
    boundary_mask = np.abs(gmsh_r - r_cap) < 1e-6

    # === Map boundary gmsh nodes to existing ring nodes ===
    # Use angular matching (safer than XY rounding)
    gmsh_theta = np.arctan2(gmsh_pts[:, 1], gmsh_pts[:, 0])
    ring_theta = np.arctan2(ring_pts[:, 1], ring_pts[:, 0])

    # wrap to [0,2pi)
    gmsh_theta = np.mod(gmsh_theta, 2*np.pi)
    ring_theta = np.mod(ring_theta, 2*np.pi)

    boundary_map = {}
    for gi, is_boundary in enumerate(boundary_mask):
        if not is_boundary:
            continue
        # find nearest ring angle
        idx = np.argmin(np.abs((ring_theta - gmsh_theta[gi] + np.pi) % (2*np.pi) - np.pi))
        boundary_map[gi] = ring_start + idx

    # === Add new interior points ===
    base = len(pts)
    new_indices = {}
    for gi in range(len(gmsh_pts)):
        if gi in boundary_map:
            continue
        new_indices[gi] = base + len(new_indices)
        pts.append(gmsh_pts[gi].tolist())

    # === Add triangles ===
    for tri in gmsh_tri:
        global_ids = []
        for gi in tri:
            if gi in boundary_map:
                global_ids.append(boundary_map[gi])
            else:
                global_ids.append(new_indices[gi])
        cells.extend([3] + global_ids)
        cell_types.append(VTK_TRIANGLE)


# def add_unstructured_cap(pts, cells, cell_types, ring_start, n_theta, z, mesh_size=0.05):
#     """
#     Create an unstructured triangular cap whose boundary nodes are exactly the loft outer ring.
#     - ring_start: index in pts where the n_theta outer-ring nodes begin
#     - n_theta: number of nodes in that ring
#     - z: z coordinate of cap
#     - mesh_size: target element size for gmsh
#     """
#     import numpy as np

#     # Collect outer ring points (use their exact XY coords to ensure conformity)
#     ring_pts = np.array([pts[ring_start + i] for i in range(n_theta)])
#     # Ensure the ring points are planar at z
#     ring_pts[:, 2] = z

#     # Generate Gmsh triangle mesh inside polygon defined by ring_pts
#     gmsh_pts, gmsh_tri = generate_cap_with_gmsh(ring_pts, mesh_size=mesh_size)

#     # Append gmsh points to global pts (they may include the boundary points again)
#     base = len(pts)
#     pts.extend(gmsh_pts.tolist())

#     # Create a mapping from gmsh boundary nodes to existing ring node indices to avoid duplicates
#     # Find points in gmsh_pts that are on the boundary (close to radius R_cap)
#     # For robustness, match by exact XY coordinate: we know gmsh used the boundary coordinates we passed,
#     # so any gmsh point with identical XY should map to the corresponding ring node.
#     gmsh_pts_arr = np.array(gmsh_pts)
#     boundary_map = {}  # gmsh_index -> existing global index
#     tol = 1e-9
#     # Build a dict mapping (rounded x,y) -> ring global index for quick lookup
#     ring_xy_map = {}
#     for i in range(n_theta):
#         x, y = round(ring_pts[i,0], 10), round(ring_pts[i,1], 10)
#         ring_xy_map[(x,y)] = ring_start + i

#     for gi, p in enumerate(gmsh_pts_arr):
#         key = (round(p[0],10), round(p[1],10))
#         if key in ring_xy_map:
#             boundary_map[gi] = ring_xy_map[key]
#     # Append triangle cells; use mapped global indices if available, else base+gi
#     for tri in gmsh_tri:
#         global_ids = []
#         for gi in tri:
#             if gi in boundary_map:
#                 global_ids.append(boundary_map[gi])
#             else:
#                 global_ids.append(base + int(gi))
#         cells.extend([3] + global_ids)
#         cell_types.append(VTK_TRIANGLE)

# Add both caps (top and bottom)
add_unstructured_cap(pts, cells, cell_types, top_ring_start, n_theta, z_top, mesh_size=cap_mesh_size)
add_unstructured_cap(pts, cells, cell_types, bottom_ring_start, n_theta, z_bottom, mesh_size=cap_mesh_size)

# ===========================================
# BUILD GRID
# ===========================================
    
pts = np.array(pts)
cells = np.array(cells, dtype=np.int64)
cell_types = np.array(cell_types, dtype=np.uint8)
grid = pv.UnstructuredGrid(cells, cell_types, pts)

# ===========================================
# DISPLAY
# ===========================================
p = pv.Plotter(window_size=(900,700))
p.add_mesh(grid, color=(0.85,0.05,0.05), smooth_shading=True, show_edges=True)
p.set_background("white")
p.add_axes()
p.camera_position = 'yz'
p.show(title="Coca-Cola Can – Unstructured Triangular Caps (Gmsh)")
