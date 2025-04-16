# Cantilever Beam with Rigid Contact Tutorial

This tutorial demonstrates how to solve a hyperelastic cantilever beam problem with a rigid circular indenter using FEniCSx. The code simulates a 2D cantilever beam that deforms as a rigid circle is pressed onto it from above.

## Problem Description

The simulation models a 2D cantilever beam made of a hyperelastic material (Neo-Hookean) with the following setup:
- The beam is fixed at the left end (cantilever boundary condition)
- A rigid circular indenter is pressed onto the beam from above
- A penalty contact formulation is used to model the interaction between the beam and the indenter
- The indentation depth increases gradually over time

## Installation

To run this code, you need to have FEniCSx and its dependencies installed:

1. Create a conda environment:
```bash
module load miniconda
mamba create -n fenicsx-env
mamba activate fenicsx-env
mamba install -c conda-forge fenics-dolfinx mpich pyvista
pip install imageio
pip install gmsh
pip install PyYAML
```

2. Activate the environment:
```bash
mamba activate fenicsx-env
```

## Running the Simulation

1. Make sure you have activated the FEniCSx environment:
```bash
mamba activate fenicsx-env
```

2. Run the simulation:
```bash
python cantilever_contact.py
```

3. The results will be saved as an animated GIF file named `cantilever_contact.gif` in the same directory.

## Customizing the Simulation

You can modify the following parameters in the code to customize the simulation:

### Geometry Parameters
- `length` and `height`: Dimensions of the cantilever beam
- `nx` and `ny`: Mesh refinement in x and y directions

### Material Parameters
- `E`: Young's modulus
- `nu`: Poisson's ratio

### Indenter Parameters
- `circle_radius`: Radius of the rigid circular indenter
- `circle_center_x` and `circle_center_y`: Initial position of the indenter
- `max_indent`: Maximum indentation depth
- `num_steps`: Number of incremental load steps

### Contact Parameters
- `penalty_coefficient`: Stiffness of the penalty contact formulation

## Understanding the Code

The code is organized into the following main sections:

1. **Mesh and Function Space Setup**: Defines the beam geometry and creates a finite element function space for the displacement field.

2. **Material Model**: Implements a Neo-Hookean hyperelastic constitutive model using the Unified Form Language (UFL).

3. **Boundary Conditions**: Sets up the cantilever boundary condition (fixed at the left end).

4. **Contact Formulation**: Implements a penalty contact model between the beam and the rigid circular indenter.

5. **Weak Form**: Defines the weak form of the problem, including elastic forces, contact forces, and external forces.

6. **Solver**: Sets up and configures the nonlinear solver for the problem.

7. **Visualization**: Prepares PyVista visualization of the deformation and the indenter.

8. **Time-stepping Loop**: Gradually increases the indentation depth and solves the nonlinear problem at each step.

## Notes

- The penalty contact method allows small penetrations between the beam and the indenter. Increasing the penalty parameter reduces penetration but may make the problem more difficult to solve.
- For large deformations or contact problems with friction, more advanced contact formulations may be needed.
- If the solver fails to converge, try:
  - Increasing the number of steps (`num_steps`)
  - Adjusting the penalty coefficient
  - Modifying the solver parameters