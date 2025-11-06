# Computational Fluid Dynamics  
## Comprehensive Finite Volume Solver

![MATLAB](https://img.shields.io/badge/Language-MATLAB-blue)
![FVM](https://img.shields.io/badge/Numerical%20Method-Finite%20Volume-orange)
![CFD](https://img.shields.io/badge/Topic-Computational%20Fluid%20Dynamics-brightgreen)
![University](https://img.shields.io/badge/AUB-American%20University%20of%20Beirut-red)

> **Author:** Charbel Ziade  
> **Instructor:** Prof. Fadl Moukalled  
> **Course:** Computational Fluid Dynamics  
> **Program:** Master of Engineering in Mechanical Engineering (Energy Track)  
> **Semester:** Spring 2024  

---

## Project Overview

This repository contains a step-by-step development of a **general 2D CFD solver** using the **Finite Volume Method (FVM)** in MATLAB.  
It follows the progression of the course through four structured assignments, each introducing a new physical or numerical complexity:

| Assignment | Focus | Key Concept |
|-------------|--------|-------------|
| **1** | Diffusion Term Discretization | Steady 2-D heat conduction, orthogonal & non-orthogonal grids |
| **2** | Convection Term Discretization | Upwind & SMART schemes for advection transport |
| **3** | Unsteady Diffusion | Time integration with Backward Euler & Crank–Nicolson |
| **4** | Incompressible Flow Solver | SIMPLE algorithm for coupled pressure–velocity–temperature fields |

Each folder includes its **own README** explaining the governing equations, discretization method, numerical schemes, and results.

---

## Numerical Framework

All solvers are implemented using the **Finite Volume Method**, ensuring conservation at the control-volume level.  

### General Equation
\[
\frac{\partial (\rho \phi)}{\partial t} + \nabla \cdot (\rho \vec{v}\phi) = \nabla \cdot (\Gamma \nabla \phi) + S_\phi
\]

Where φ represents the transported scalar (temperature, velocity, etc.).  

---

## Assignment Summaries

### [Assignment 1 – Diffusion Term Discretization](README_DiffusionTerm.md)
- Develops a **2-D steady diffusion solver**.  
- Validates with the **orthogonal benchmark case** (Moukalled et al., 2016, Example 5.2 p. 227).  
- Extends to **non-orthogonal structured meshes** using **Transfinite Interpolation** (TFI).  
- Solved using **Gauss–Seidel** and **TDMA** iterative methods.  


📂 Folder: `Diffusion Term/`  
📄 Detailed README: [README_DiffusionTerm.md](README_DiffusionTerm.md)

---

### [Assignment 2 – Convection Term Discretization](README_Advection.md)
- Solves **advection** with a **known velocity field**.  
- Implements:
  - **First-order Upwind Scheme** (stable but diffusive)  
  - **Modified SMART Scheme** (accurate & bounded via deferred correction)  
- Benchmarked on the **transverse step-profile convection test**.  
- Compares both schemes along the domain centerline.

📂 Folder: `Advection Term/`  
📄 Detailed README: [README_Advection.md](README_Advection.md)

---

### [Assignment 3 – Transient (Unsteady) Heat Conduction](README_TransientTerm.md)
- Adds **time-dependence** to the diffusion equation.  
- Implements two time integration schemes:
  - **Backward Euler (1st order implicit)**  
  - **Crank–Nicolson (2nd order accurate)**  
- Solves the **2-D transient conduction problem** with adiabatic and Dirichlet boundaries.  
- Validated against analytical solutions at t = 10, 30, 60 s.

Folder: `Transient Term/`  
Detailed README: [README_TransientTerm.md](README_TransientTerm.md)

---

### [Assignment 4 – Incompressible Flow (SIMPLE Algorithm)](README_IncompFlow.md)
- Final solver combining **momentum, continuity, and energy** equations.  
- Uses **SIMPLE algorithm** with **Rhie–Chow interpolation** to couple pressure and velocity.  
- Employs **SMART convection** and **central differencing diffusion**.  
- Calculates full hydrodynamic and thermal fields for air.  
- Demonstrates convergence and pressure-velocity consistency.

Folder: `Incompressible Flow Equation/`  
Detailed README: [README_IncompFlow.md](README_IncompFlow.md)

---

## Numerical Techniques Used

| Category | Techniques Implemented |
|-----------|------------------------|
| **Spatial Discretization** | Central Differencing, Upwind, SMART (Modified SMART), Deferred Correction |
| **Temporal Discretization** | Backward Euler, Crank–Nicolson |
| **Pressure–Velocity Coupling** | SIMPLE Algorithm (with Rhie–Chow correction) |
| **Solvers** | Gauss–Seidel, TDMA (line-by-line), Sparse Matrix Assembly |
| **Mesh Handling** | Transfinite Interpolation for Non-Orthogonal Grids |
| **Convergence Control** | Under-Relaxation, Residual Monitoring|

---

## Repository Structure

```
Computational Fluid Dynamics/
│
├── Diffusion Term/
│   ├── DiffusionTerm.m
│   ├── Orthogonal_Diffusion.m
│   ├── README_DiffusionTerm.md
│
├── Advection Term/
│   ├── Upwind_GaussSeidel.m
│   ├── Upwind_TDMA.m
│   ├── SMART_GaussSeidel.m
│   ├── ModifiedSMART_GaussSeidel.m
│   ├── README_Advection.md
│
├── Transient Term/
│   ├── Transient_Analytical.m
│   ├── Prototype_2x2.m
│   ├── CrankNicholson.m
│   ├── BACKWARD_EULER.m
│   ├── README_TransientTerm.md
│
├── Assignment4_IncompressibleFlow/
│   ├── Incomp_Flow.m
│   ├── README_IncompFlow.md
│
└── README.md   ← (this file)
```

---

## References

- Moukalled, F., Mangani, L., & Darwish, M. (2016).  
  *The Finite Volume Method in Computational Fluid Dynamics.* Springer.  


---

## Author

**Charbel Ziade**  
Master of Engineering in Mechanical Engineering (Energy Track) – American University of Beirut (AUB)  
[email: caz09@mail.aub.edu](mailto:caz09@mail.aub.edu)  
*([Linkedin: Charbel Ziade](https://www.linkedin.com/in/charbel-ziade-b908141a4/))*  
