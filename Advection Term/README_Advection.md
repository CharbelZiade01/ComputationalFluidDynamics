# 🌬️ Assignment 2 – 2D Convection Term Discretization (Finite Volume Method)

![MATLAB](https://img.shields.io/badge/Language-MATLAB-blue)
![FVM](https://img.shields.io/badge/Method-Finite%20Volume%20Method-orange)
![CFD](https://img.shields.io/badge/Field-Computational%20Fluid%20Dynamics-brightgreen)

> **Course:** Computational Fluid Dynamics  
> **Instructor:** Prof. Fadl Moukalled  
> **Author:** Charbel Ziade  
> **Institution:** American University of Beirut (AUB)  
> **Semester:** Spring 2024  

---

##  Problem Overview

In this assignment, the solver handles the **convection (advection) term** assuming a **known velocity field**.  

The governing transport equation for a scalar variable φ is:

\[
\nabla \cdot (\rho \vec{v} \phi) = 0
\]

which in integral form becomes:

\[
\sum_f (\rho \vec{v} \phi)_f \cdot \vec{S}_f = 0
\]

The main challenge is that the value of φ at the **face centroids** must be interpolated from cell-center values using suitable **discretization schemes**.

---

## Objective

- Implement **first-order upwind** and **higher-order modified SMART** convection schemes  
- Compare the two in terms of accuracy and stability  
- Validate using the **pure convection of a transverse step profile** benchmark (square domain, 30 × 25 elements)

---

## Upwind Scheme

The **upwind scheme** determines the face value of φ_f based on the flow direction:

\[
\phi_f =
\begin{cases}
\phi_\text{upstream}, & \text{if } \dot{m}_f > 0 \\
\phi_\text{downstream}, & \text{otherwise}
\end{cases}
\]

This yields the discretized equation:
\[
a_C \phi_C + a_E \phi_E + a_W \phi_W + a_N \phi_N + a_S \phi_S = b_C
\]

where an additional contribution appears in the neighbor coefficients:
\[
- \max(-\dot{m}_f, 0)
\]

While the upwind method is **stable**, it introduces **numerical diffusion**, leading to smearing of sharp gradients.

---

## SMART and Modified SMART Schemes

To improve accuracy, the **SMART (Sharp and Monotonic Algorithm for Realistic Transport)** family of schemes is used.

These are based on the **Normalized Variable Formulation (NVF)**:

\[
\tilde{\phi}_C = \frac{\phi_C - \phi_U}{\phi_D - \phi_U}
\]

where **U**, **C**, and **D** denote the upstream, current, and downstream centroids respectively.

### Standard SMART Scheme
\[
\tilde{\phi}_f =
\begin{cases}
\frac{3}{4}\tilde{\phi}_C + \frac{3}{8}, & 0 \le \tilde{\phi}_C \le \frac{5}{6} \\
1, & \frac{5}{6} \le \tilde{\phi}_C \le 1 \\
\tilde{\phi}_C, & \text{otherwise}
\end{cases}
\]

### Modified SMART Scheme
\[
\tilde{\phi}_f =
\begin{cases}
3\tilde{\phi}_C, & 0 \le \tilde{\phi}_C \le \frac{1}{6} \\
\frac{3}{4}\tilde{\phi}_C + \frac{3}{8}, & \frac{1}{6} \le \tilde{\phi}_C \le \frac{7}{10} \\
\frac{1}{3}\tilde{\phi}_C + \frac{2}{3}, & \frac{7}{10} \le \tilde{\phi}_C \le 1 \\
\tilde{\phi}_C, & \text{otherwise}
\end{cases}
\]

The **modified version** offers a better balance between **accuracy** and **stability**, limiting overshoots while sharply resolving gradients.

---

## Deferred Correction Implementation

The **deferred correction** technique combines the stability of upwind with the accuracy of SMART by adjusting the source term:

\[
b_C = b_C - \sum_f \dot{m}_f \cdot (\phi_f^{SMART} - \phi_f^{UPWIND})
\]

Implementation details:
- Loops are performed over **faces**, not cells.  
- Separate loops handle internal east and north faces.  
- Boundary faces are treated individually.

---

## Boundary Conditions

- **West & South:** Inlet (known φ) → contributes to a_C and b_C  
- **East & North:** Outlet (zero-flux) → contributes only to b_C, usually zero  

---

## Results

- Simulations performed on a **30 × 25 mesh**.  
- **Upwind**: stable but overly diffused results.  
- **Modified SMART**: accurately preserves the step profile with minimal oscillation.  
- Comparison along the vertical centerline validates expected trends versus analytical results.

> *Visualization of scalar field available in* `Advection.fig`.

---

## Code Structure

```
Advection Term/
   ├── Upwind_GaussSeidel.m
   ├── Upwind_TDMA.m
   ├── SMART_GaussSeidel.m
   ├── ModifiedSMART_GaussSeidel.m
   ├── README_Advection.md
```

---

## References

- Moukalled, F., Mangani, L., & Darwish, M. (2016). *The Finite Volume Method in Computational Fluid Dynamics*. Springer.(https://link.springer.com/chapter/10.1007/978-3-319-16874-6_5#Sec21)  


---

## Author

**Charbel Ziade**  
Master of Engineering in Mechanical Engineering (Energy Track) – American University of Beirut (AUB)  
[email: caz09@mail.aub.edu](mailto:caz09@mail.aub.edu)  
*([Linkedin: Charbel Ziade](https://www.linkedin.com/in/charbel-ziade-b908141a4/))*  
