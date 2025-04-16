# cantilever_contact.py - Tutorial on using FEniCSx for cantilever beam contact problem
"""
This tutorial demonstrates how to solve a hyperelastic cantilever beam problem
with a rigid circular indenter using FEniCSx. The beam is fixed at the left end
and the rigid circle is pressed onto the beam from the top.
"""

from dolfinx import default_scalar_type, fem, log, mesh, plot
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.fem.petsc import NonlinearProblem
from mpi4py import MPI
import numpy as np
import pyvista
import ufl

# -------------------------------
# Define the beam geometry and mesh
# -------------------------------
# Parameters for the beam
length = 10.0     # Length of the beam
height = 1.0      # Height of the beam
nx, ny = 100, 10  # Number of elements in each direction

# Create a rectangular mesh for the beam
domain = mesh.create_rectangle(
    MPI.COMM_WORLD,
    [[0.0, 0.0], [length, height]],
    [nx, ny],
    cell_type=mesh.CellType.triangle
)

# -------------------------------
# Define a function space over the domain
# -------------------------------
# Vector function space for displacement field (2 components in 2D)
V = fem.functionspace(domain, ("Lagrange", 2, (domain.geometry.dim,)))
v = ufl.TestFunction(V)          # Test function
u = fem.Function(V, name="Displacement")  # Solution function

# -------------------------------
# Material properties for hyperelasticity
# -------------------------------
# Young's modulus and Poisson's ratio
E = 1e5
nu = 0.3
# Lame parameters
mu = fem.Constant(domain, E / (2 * (1 + nu)))  # Shear modulus
lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))  # First Lame parameter

# -------------------------------
# Hyperelastic Constitutive Model (Neo-Hookean)
# -------------------------------
d = domain.geometry.dim
I = ufl.Identity(d)                    # Identity tensor
F = ufl.variable(I + ufl.grad(u))      # Deformation gradient
J = ufl.det(F)                         # Volume ratio
C = F.T * F                            # Right Cauchy-Green tensor
I1 = ufl.tr(C)                         # First invariant of C

# Strain energy for compressible Neo-Hookean material
psi = (mu / 2) * (I1 - 3) - mu * ufl.ln(J) + (lmbda / 2) * ufl.ln(J)**2

# First Piola-Kirchhoff stress tensor
P = ufl.diff(psi, F)

# -------------------------------
# Define the rigid circle indenter
# -------------------------------
# Parameters for the circle
circle_radius = 0.5
circle_center_x = length / 2.0  # Place it in the middle of the beam
circle_center_y = height + circle_radius * 0.5  # Initial position slightly above the beam

# Create a controllable indenter displacement
indent_depth = fem.Constant(domain, 0.0)  # Will be updated in the time loop

# -------------------------------
# Mark boundaries for boundary conditions
# -------------------------------
def left_boundary(x): return np.isclose(x[0], 0.0)
def top_boundary(x): return np.isclose(x[1], height)

# Find facets on the boundaries
fdim = domain.topology.dim - 1  # Dimension of facets (1 for 2D mesh)
left_facets = mesh.locate_entities_boundary(domain, fdim, left_boundary)
top_facets = mesh.locate_entities_boundary(domain, fdim, top_boundary)

# Mark facets with different tags
marked_facets = np.hstack([left_facets, top_facets])
marked_values = np.hstack([
    np.full_like(left_facets, 1),  # Left boundary tag: 1
    np.full_like(top_facets, 2)    # Top boundary tag: 2
])
sorted_facets = np.argsort(marked_facets)
facet_tag = mesh.meshtags(domain, fdim, marked_facets[sorted_facets], marked_values[sorted_facets])

# -------------------------------
# Set up boundary conditions
# -------------------------------
# Fixed displacement at left boundary (cantilever)
u_bc = np.array((0.0, 0.0), dtype=default_scalar_type)  # Zero displacement
left_dofs = fem.locate_dofs_topological(V, facet_tag.dim, facet_tag.find(1))
bcs = [fem.dirichletbc(u_bc, left_dofs, V)]

# -------------------------------
# Define contact formulation
# -------------------------------
# Penalty method parameters
penalty_coefficient = fem.Constant(domain, 1e6)  # Penalty stiffness

# Contact function using a smoothed penalty approach
def contact_penalty(u):
    """
    Calculate contact forces between beam and circular indenter using penalty method
    
    Args:
        u: Displacement field
    
    Returns:
        Vector field of contact forces
    """
    # Get current coordinates (reference + displacement)
    x = ufl.SpatialCoordinate(domain)
    current_pos = x + u
    
    # Calculate position of indenter center (with indentation)
    indenter_x = circle_center_x
    indenter_y = circle_center_y - indent_depth.value
    
    # Calculate distance vector from current point to indenter center
    dx = current_pos[0] - indenter_x
    dy = current_pos[1] - indenter_y
    
    # Calculate distance to indenter surface (negative means penetration)
    r = ufl.sqrt(dx**2 + dy**2)
    distance = r - circle_radius
    
    # Normal direction (pointing outward from circle)
    normal_x = dx / (r + 1e-14)  # avoid division by zero
    normal_y = dy / (r + 1e-14)
    normal = ufl.as_vector([normal_x, normal_y])
    
    # Smoothed penetration function
    # negative value indicates penetration
    penetration = ufl.min_value(distance, 0.0)
    
    # Penalty force is proportional to penetration depth
    # Direction is along the normal, pointing outward from the circle
    penalty_force = penalty_coefficient * penetration * normal
    
    return penalty_force

# -------------------------------
# Define forces and weak form
# -------------------------------
# Body force (gravity or other distributed load)
B = fem.Constant(domain, default_scalar_type((0.0, -0.1)))  # Small gravity force

# Integration measures
metadata = {"quadrature_degree": 4}
dx = ufl.Measure("dx", domain=domain, metadata=metadata)
ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tag, metadata=metadata)

# Weak form: internal elastic forces + contact forces - external forces = 0
F_internal = ufl.inner(ufl.grad(v), P) * dx           # Internal elastic forces
F_external = ufl.inner(v, B) * dx                     # External body forces
F_contact = ufl.inner(v, contact_penalty(u)) * dx     # Contact forces

# Total weak form
F_form = F_internal - F_external + F_contact

# -------------------------------
# Set up nonlinear solver
# -------------------------------
problem = NonlinearProblem(F_form, u, bcs)
solver = NewtonSolver(domain.comm, problem)
solver.atol = 1e-8
solver.rtol = 1e-8
solver.convergence_criterion = "incremental"
solver.max_it = 25  # Increased max iterations for contact problems

# -------------------------------
# Visualization setup
# -------------------------------
pyvista.start_xvfb()
plotter = pyvista.Plotter(window_size=(800, 600), off_screen=True, notebook=False)
plotter.open_gif("cantilever_contact.gif", fps=3)

# Create the visualization mesh
topology, cells, geometry = plot.vtk_mesh(u.function_space)
function_grid = pyvista.UnstructuredGrid(topology, cells, geometry)

# Prepare displacement field for visualization (3D vectors for PyVista)
values = np.zeros((geometry.shape[0], 3))
function_grid["u"] = values
function_grid.set_active_vectors("u")

# Create warped mesh for deformation visualization
warped = function_grid.warp_by_vector("u", factor=1.0)
warped.points[:, 2] = 0.0  # Ensure all points stay in the z=0 plane
beam_actor = plotter.add_mesh(warped, show_edges=True, color="lightblue")

# Create a scalar field for displacement magnitude visualization
Vs = fem.functionspace(domain, ("Lagrange", 2))
magnitude = fem.Function(Vs)
us = fem.Expression(ufl.sqrt(sum([u[i]**2 for i in range(len(u))])), 
                    Vs.element.interpolation_points())

# Add the rigid circle indenter to the plot
circle = pyvista.Circle(radius=circle_radius, resolution=50)
circle.translate([circle_center_x, circle_center_y, 0.0], inplace=True)
circle_actor = plotter.add_mesh(circle, color="red")

# -------------------------------
# Time-stepping loop to gradually apply indentation
# -------------------------------
log.set_log_level(log.LogLevel.INFO)
max_indent = 1.0   # Maximum indentation depth
num_steps = 10     # Number of load steps

for step in range(1, num_steps + 1):
    # Update indentation depth
    indent_depth.value = step * max_indent / num_steps
    print(f"\nStep {step}/{num_steps} - Indentation depth: {indent_depth.value}")
    
    # Update circle position for visualization
    circle.translate([0.0, -max_indent / num_steps, 0.0], inplace=True)
    
    # Solve the nonlinear problem
    try:
        num_its, converged = solver.solve(u)
        assert converged, f"Newton solver did not converge at step {step}"
        u.x.scatter_forward()
        print(f"  Newton iterations: {num_its}")
    except Exception as e:
        print(f"  Solver failed: {e}")
        break
    
    # Update visualization
    function_grid["u"][:, :len(u)] = u.x.array.reshape(geometry.shape[0], len(u))
    magnitude.interpolate(us)
    function_grid["mag"] = magnitude.x.array
    
    warped_n = function_grid.warp_by_vector("u", factor=1.0)
    warped.points[:, :] = warped_n.points
    warped.set_active_scalars("mag")
    
    plotter.update_scalar_bar_range([0, np.max(magnitude.x.array)])
    if step == 1:
        plotter.view_xy()
        plotter.reset_camera()
        # Zoom out a bit to see both beam and indenter
        plotter.camera.Zoom(0.8)
    
    plotter.write_frame()

plotter.close()
print("\nSimulation completed successfully.")