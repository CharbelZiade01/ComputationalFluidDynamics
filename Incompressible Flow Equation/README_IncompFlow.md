# Assignment 4 – 2D Incompressible Flow Solver (SIMPLE Algorithm)

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

In this final CFD assignment, the solver developed in previous tasks (diffusion, convection, and transient heat conduction) is extended to simulate **two-dimensional incompressible flow** using the **SIMPLE algorithm** (*Semi-Implicit Method for Pressure-Linked Equations*).  

The goal is to obtain the **hydrodynamic (u, v, p)** and **thermal (T)** fields of an air flow within a rectangular domain, using the **SMART scheme** for convection.  

### Fluid Properties
| Property | Symbol | Value |
|-----------|---------|--------|
| Density | ρ | 0.8 kg/m³ |
| Dynamic Viscosity | μ | 5 × 10⁻⁵ Pa·s |
| Specific Heat | cₚ | 1.03 kJ/kg·K |
| Thermal Conductivity | k | 0.036 W/m·K |

---

## Governing Equations

The system of equations for steady, laminar, incompressible flow consists of:

### **Continuity:**
\[
\nabla \cdot \vec{v} = 0
\]

### **Momentum (Navier–Stokes):**
\[
\rho (\vec{v}\cdot\nabla)\vec{v} = -\nabla p + \mu \nabla^2 \vec{v}
\]

### **Energy:**
\[
\rho c_p (\vec{v}\cdot\nabla T) = k \nabla^2 T
\]

---

## Numerical Method: SIMPLE Algorithm

Because the **continuity equation does not explicitly contain p**, a **pressure-correction** approach is required to couple velocity and pressure.

### Algorithm Steps
1. Initialize u, v, p, and mass flow rate (ṁ) with guessed fields.  
2. **Solve momentum equations** (for u and v) → obtain intermediate velocities u*, v*.  
3. **Apply Rhie–Chow interpolation** to compute face velocities and mass-flow rates → avoid checkerboard pressure.  
4. **Assemble and solve the pressure-correction equation** for p′.  
5. **Update fields:**
   \[
   u = u^* + d_u (p'_{west} - p'_{east}), \quad
   v = v^* + d_v (p'_{south} - p'_{north})
   \]
   \[
   p = p + \alpha_p p'
   \]
6. **Recalculate mass-flow rates and repeat** until convergence.  
7. **Solve energy equation** using updated velocity field (advection + diffusion).  

### Rhie–Chow Interpolation
Prevents pressure–velocity decoupling on collocated grids:
\[
\vec{v}_f = \bar{\vec{v}} - \frac{1}{a_P}(\nabla p - \overline{\nabla p})
\]

### Pressure Correction Equation
Derived from continuity + momentum:
\[
a_P p'_P = a_E p'_E + a_W p'_W + a_N p'_N + a_S p'_S + b_P
\]
where  
\[
b_P = \sum_f \dot{m}_f^*
\]
and under-relaxation is applied to p and velocity updates for stability.  

---

## Implementation Structure

Each equation (momentum, pressure, and energy) is solved sequentially using **finite volume discretization** and **SMART scheme** for convection.  
Under-relaxation factors are used to enhance stability and ensure smooth convergence.

### **Solver Loop**
```
Initialize (u, v, p, T)
│
├─► Solve Momentum (u*, v*)
│
├─► Rhie–Chow interpolation → compute face velocities
│
├─► Assemble & solve Pressure Correction (p')
│
├─► Correct (u, v, p) fields
│
├─► Solve Energy equation (advection + diffusion)
│
└─► Repeat until residuals < 10⁻⁶
```

---

## Discretization Summary

| Term | Scheme | Notes |
|------|---------|-------|
| Convection | SMART (upwind-biased high-resolution) | Sharp gradients, no overshoot |
| Diffusion | Central Differencing | Second-order accurate |
| Pressure–Velocity Coupling | SIMPLE + Rhie–Chow | Collocated grid correction |
| Under-Relaxation | Velocity = 0.7, Pressure = 0.3 | Ensures convergence |

---

## Results Summary

- The SIMPLE solver converged to stable velocity and pressure fields.  
- The flow pattern satisfied both **momentum** and **continuity** equations.  
- **Temperature field** distribution was influenced by both advection and diffusion.  
- Rhie–Chow interpolation successfully removed pressure oscillations.  
- Residual plots showed smooth decay → numerical stability achieved.  

---

## File Structure

```
Incompressible Flow Equation/
│   ├── Incomp_Flow.m
│   ├── README_IncompFlow.md
```

---

## References

- Moukalled, F., Mangani, L., & Darwish, M. (2016). *The Finite Volume Method in Computational Fluid Dynamics.* Springer. [Ch. 6 – Incompressible Flow and SIMPLE Algorithm](https://link.springer.com/chapter/10.1007/978-3-319-16874-6_5#Sec21)  
 

---

## Author

**Charbel Ziade**  
Master of Engineering in Mechanical Engineering (Energy Track) – American University of Beirut (AUB)  
[email: caz09@mail.aub.edu](mailto:caz09@mail.aub.edu)  
*([Linkedin: Charbel Ziade](https://www.linkedin.com/in/charbel-ziade-b908141a4/))*  

