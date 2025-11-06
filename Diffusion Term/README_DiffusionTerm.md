# 🧮 Assignment 1 – 2D Diffusion Term Discretization (Finite Volume Method)

![MATLAB](https://img.shields.io/badge/Language-MATLAB-blue)
![FVM](https://img.shields.io/badge/Method-Finite%20Volume%20Method-orange)
![CFD](https://img.shields.io/badge/Field-Computational%20Fluid%20Dynamics-brightgreen)

> **Course:** Computational Fluid Dynamics  
> **Instructor:** Prof. Fadl Moukalled  
> **Author:** Charbel Ziade  
> **Institution:** American University of Beirut (AUB)  
> **Semester:** Spring 2024  

---

## Problem Overview

This assignment develops and verifies a **finite-volume solver** for 2-D diffusion-type problems.  
Two cases are addressed:

1. A **benchmark orthogonal problem** reproduced from *Moukalled et al.*, 2016 (page 227) to validate the solver on a uniform grid.  
2. A **general non-orthogonal diffusion problem** with variable Γ and S.

---

## 🧩 Part 1 – Orthogonal Diffusion Benchmark (*Moukalled et al.*, 2016, Example 5.2, p. 227*)

### Mathematical Model
\[
-∇·(Γ∇T)=0
\]

**Boundary Conditions**

| Boundary | Type | Value |
|-----------|-------|-------|
| West | Dirichlet | T = 100 °C |
| East | Dirichlet | T = 0 °C |
| North | Neumann (adiabatic) | ∂T/∂y = 0 |
| South | Neumann (adiabatic) | ∂T/∂y = 0 |

**Domain & Parameters**

| Parameter | Symbol | Value |
|------------|---------|-------|
| Width | Lx | 1 m |
| Height | Ly | 1 m |
| Diffusion coefficient | Γ | 1 W/m·K |
| Grid | 20 × 20 uniform cells |

---

### Discretization and Numerical Method

- Finite Volume Method with **central differencing** for diffusion.  
- **Uniform orthogonal grid** → no non-orthogonal correction needed.  
- Linear algebraic form:
  \[
  a_C T_C = a_E T_E + a_W T_W + a_N T_N + a_S T_S + b_C
  \]
- Sparse matrix solved iteratively using Gauss–Seidel and TDMA (line-by-line).

---

### Results and Validation

- Numerical solution shows a **linear temperature gradient** from left to right.  
- Comparison with analytical solution:
  \[
  T(x)=100(1-x/L_x)
  \]
- Mean absolute error < 1 % for 20 × 20 grid.  
- Confirms correct coefficient assembly and boundary implementation.

> 📄 Implemented in `Orthogonal_Diffusion.m`

---

## 🧭 Part 2 – General Non-Orthogonal Diffusion Problem

### Governing Equation
\[
-∇·(Γ∇T)=S
\]
with variable coefficients:
\[
S=T \frac{(2x − 0.2y)}{400},\quad
Γ=T \frac{(x^2 + e^{0.25y})}{400}
\]

### Mesh and Geometry

- **Transfinite Interpolation (TFI)** for non-orthogonal structured mesh.  
- Computed cell centroids, face centroids, areas, and volumes.  
- Normal vector S split into E (aligned) and T (transverse) components.  
- Geometric ratios:
  \[
  g_\text{diffusion}=|E|/Δ_{PC},\quad g_f=V_P/(V_P+V_C)
  \]
- Visualization → `Diffusion_Mesh.fig`

---

### Numerical Implementation

- Non-linear problem solved iteratively; Γ(T) and S(T) updated each iteration.  
- Sparse matrix A of size [(nx − 1)(ny − 1)]² assembled.  
- Solvers: Gauss–Seidel and TDMA.  
- Under-relaxation (α = 1 for this case).  
- Grid refinement 20×20 → 30×30 → 40×40 → grid-independent solution at 40×40.

---

### Results Summary

- Smooth temperature distribution consistent with variable Γ and S.  
- Stable and monotonic convergence observed.  
- Confirms robustness of the non-orthogonal formulation.

---

## 🧾 Code Structure

```
Assignment1_Diffusion/
│
├── orthogonal_case_227.m      # textbook orthogonal benchmark (Moukalled p. 227)
├── diffusion_solver.m         # non-orthogonal diffusion solver
├── mesh_generator.m           # transfinite interpolation mesh
├── boundary_conditions.m      # Dirichlet/Neumann setup
├── gauss_seidel.m             # iterative solver
├── tdma_solver.m              # 1-D TDMA solver
├── Diffusion_Mesh.fig         # sample mesh plot
└── README_Assignment1.md      # this file
```

---

## 🧠 References

- Moukalled, F., Mangani, L., & Darwish, M. (2016).  
  *The Finite Volume Method in Computational Fluid Dynamics.* Springer. (Ch. 5, Example 5.2, p. 227)  
- Lecture Notes – MECH 663: Diffusion Term Discretization  
- MATLAB Documentation – Sparse Matrices and Vectorization  

---

## 🖋️ Author

**Charbel Ziade**  
MSc Mechanical Engineering (Energy Track) – AUB  
📧 [charbelziade2000@hotmail.com](mailto:charbelziade2000@hotmail.com)  
🌐 *(LinkedIn link to be added)*  
