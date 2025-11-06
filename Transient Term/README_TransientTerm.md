# ⏱️ Assignment 3 – 2D Transient Heat Conduction (Unsteady Diffusion)

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

In this assignment, the solver includes the **transient (unsteady)** term.  
The problem models **2D unsteady heat conduction** governed by:

\[
\frac{\partial T}{\partial t} = \alpha \left( \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} \right)
\]

where  
- **T(x, y, t)**: temperature field  
- **α**: thermal diffusivity  

Two time-integration schemes are implemented:
- **First-order implicit (Backward Euler)**  
- **Second-order Crank–Nicolson**  

---

## Domain & Setup

- **Geometry:** rectangular plate of size *a = 0.1 m, b = 0.2 m*  
- **Grid:** 20 × 40 control volumes  
- **Time step:** Δt = 1 s  
- **Simulation duration:** 0 ≤ t ≤ 60 s  
- **Comparison times:** t = 10, 30, and 60 s  
- **Analytical solution available for validation**

---

## Boundary Conditions

- **North & South:** Dirichlet (T = 0 K)  
- **East & West:** Neumann (zero heat flux)  
- **Initial condition:** T = 1000 K throughout the domain  

---

## Steady-State Reference

Before adding the transient term, the steady-state diffusion equation (from Assignment 1) is solved to verify setup correctness.  
The steady solution yields **T = 0 K everywhere**, as expected for these boundary conditions.

---

## Fully Implicit Scheme – Backward Euler

The **Backward Euler method** is a first-order accurate, unconditionally stable scheme.  
The discretized form for each control volume becomes:

\[
a_C T_C^{n+1} = a_E T_E^{n+1} + a_W T_W^{n+1} + a_N T_N^{n+1} + a_S T_S^{n+1} + b_C + \rho c_p \frac{V_P}{\Delta t} T_C^n
\]

Implementation details:
- Extend temperature to 3D array: *T(i, j, n)*  
- Loop over time steps updating coefficients and source terms  
- Store and plot temperature at selected times

### Example results (2×2 prototype test):
| Time (s) | Avg T (K) |
|-----------|-----------|
| 1         | 999.23    |
| 30        | 977.95    |
| 60        | 955.65    |

The temperature decays gradually as heat diffuses to the cold boundaries.

---

## Second-Order Scheme – Crank–Nicolson

The **Crank–Nicolson scheme** provides **second-order accuracy** in both time and space.  
It averages the diffusion term between current and previous time steps:

\[
\frac{T^{n+1} - T^n}{\Delta t} = \frac{1}{2} \left[ \alpha \nabla^2 T^{n+1} + \alpha \nabla^2 T^n \right]
\]

Implemented using a **two-step approach**:
1. Divide the time step Δt by 2.  
2. Apply correction using intermediate results from previous steps.

The scheme is slightly more computationally expensive but yields smoother and more accurate transient profiles compared to Backward Euler.

---

## Results and Discussion

- Both Backward Euler and Crank–Nicolson yield consistent cooling behavior.  
- The **temperature field decreases with time**, approaching steady state (T = 0 K).  
- Crank–Nicolson shows slightly higher accuracy and reduced numerical diffusion.  
- Analytical and numerical solutions agree closely at t = 10, 30, and 60 s.

---

## Code Structure

```
Transient Term/
   ├── Transient_Analytical.m
   ├── Prototype_2x2.m
   ├── CrankNicholson.m
   ├── BACKWARD_EULER.m
   ├── README_TransientTerm.md

```

---

## References

- Moukalled, F., Mangani, L., & Darwish, M. (2016). *The Finite Volume Method in Computational Fluid Dynamics*. Springer. [Ch. 5 – Transient Diffusion](https://link.springer.com/chapter/10.1007/978-3-319-16874-6_5#Sec21)  


---

## Author

**Charbel Ziade**  
Master of Engineering in Mechanical Engineering (Energy Track) – American University of Beirut (AUB)  
[email: caz09@mail.aub.edu](mailto:caz09@mail.aub.edu)  
*([Linkedin: Charbel Ziade](https://www.linkedin.com/in/charbel-ziade-b908141a4/))*  
