# FEniCSx Finite Element Analysis Tutorials

A collection of Jupyter notebooks demonstrating structural mechanics simulations using the FEniCSx finite element framework.

## Description

This repository provides educational examples of nonlinear finite element analysis for structural mechanics problems. The tutorials range from basic beam deflection to advanced convergence analysis, showcasing the capabilities of FEniCSx for solving complex mechanical problems.

## Installation

To run these notebooks, you need to have FEniCSx and its dependencies installed:

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

## Running the Notebooks

1. Make sure you have activated the FEniCSx environment:
```bash
mamba activate fenicsx-env
```

2. Launch Jupyter Lab:
```bash
jupyter lab
```

3. Open and run the notebooks

## Notebook Descriptions

### cantilever_circular_distributed_load.ipynb
- **Description:** 3D beam deflection simulation with a Neo-Hookean hyperelastic material model
- **Features:** Fixed beam at one end with circular distributed force at the other end
- **Output:** Animated visualization of beam deflection (3d_beam_deflection.gif)

### tutorial_1_analytical.ipynb
- **Description:** Validation of beam deflection against analytical solutions
- **Features:** Comparison between FEniCSx numerical results and engineering beam theory
- **Output:** Comparative plots and error analysis

### tutorial_2_mesh_refinement.ipynb
- **Description:** Mesh refinement studies for a cantilever beam problem
- **Features:** Both h-refinement (mesh density) and p-refinement (polynomial order)
- **Output:** Convergence plots and mesh quality metrics

### tutorial_3_failure_conditions.ipynb
- **Description:** Analysis of nonlinear solver convergence behavior
- **Features:** Investigation of convergence failure under different load magnitudes
- **Output:** Visualization of convergence boundaries and solver behavior

## Key Features

- Nonlinear hyperelastic material modeling
- Large deformation mechanics
- Incremental loading techniques
- Convergence analysis
- Mesh refinement studies
- Visualization of results
- Validation against analytical solutions

## Notes

- Each notebook includes detailed explanations of the problem setup, solution approach, and results
- For large models, computation time may be significant
- The PyVista visualization creates an animated GIF showing the deflection over time

These tutorials provide a comprehensive introduction to nonlinear finite element analysis using FEniCSx, from basic validation to advanced convergence studies.