---
title: 'Webster equation: pH formulation and discretization'
---


The Webster's equation (or horn equation) describes the plane-wave propagation in acoustic ducts of space-dependent section. A time-varying geometry is also considered in this work, mainly for use in articulations of the vocal tract.

As this component is built mainly to be used in the vocal tract, visco-thermal losses are ignored. Combined with the plane wave assumption, the validity range of the model should be ~0-6000 Hz for the vocal tract.

## Notations

| Constant       | Notation | Dimension  | Typical value (in air) |
| -------------- | -------- | ---------- | ---------------------- |
| Rest density   | $\rho_0$ | $M.L^{-3}$ | $1.2\ kg.m^{-3}$       |
| Speed of sound | $c_0$    | $L.T^{-1}$ | $340\ m.s^{-1}$        |
| Tube length    | $l_0$    | $L$        | $17\ cm$ (vocal tract) |

| Time-space dependent variable | Notation | Dimension         | S.I.unit    |
| ----------------------------- | -------- | ----------------- | ----------- |
| Density                       | $\rho$   | $M.L^{-3}$        | $kg.m^{-3}$ |
| Pressure                      | $P$      | $M.L^{-1}.T^{-2}$ | $Pa$        |
| Velocity (axial)              | $v$      | $L.T^{-1}$        | $m.s^{-1}$  |
| Cross-section area            | $A$      | $L^2$             | $m^2$       |
| Perimeter                     | $\Gamma$ | $L$               | $m$         |

## PDE formulation

The system considered directly follow from the reduction of Euler equations to 1D by considering homogenous fields on cross-sections and normal velocity continuity at the tube closed boundaries. The kinetic energy associated with the non-axial velocity is neglected and the system partially linearized.

The state reduces to $x = [v, \rho, A]^\intercal$, associated with Hamiltonian
$$
    \mathcal H (x) = \frac{1}{2} \left(\rho_0 \vert v \vert^2_A + \frac{c_0^2}{\rho_0} \vert \rho \vert_A^2\right),
$$
with norms meant over the length of the tube and weighted by the cross section area $A$.

Dynamics are
$$
\begin{equation}
    \begin{bmatrix}
        \partial_t v \\\
        \partial_t \rho \\\
        \partial_t A
    \end{bmatrix}
    =
    \begin{bmatrix}
        0 &-\partial_x\left(\frac{\bullet}{A}\right) & 0 \\\
        -\frac{1}{A} \partial_x & 0 & 0 \\\
        0 & 0 & 0
    \end{bmatrix}
    \underbrace{
    \begin{bmatrix}
        A \rho_0 v \\\
        A \frac{c_0^2}{\rho_0} \rho \\\
        \frac{1}{2} \left(\rho_0 v^2 + \frac{c_0^2 \rho^2}{\rho_0}\right)
    \end{bmatrix}
    }_{\nabla \mathcal H}
    +
    \begin{bmatrix}
        0 \\\
        - \frac{\rho_0}{A} \\\
        1
    \end{bmatrix}
    \underbrace{\partial_t A}_u,
    \label{base_webster}
\end{equation}
$$
with cross-section area variations as a distributed input to the system. 
The effective pressure applied to the walls is the power-conjugated output retrieved by skew-symmetry:

$$
\underbrace{P}_y =  
-
\begin{bmatrix}
    0 &
    - \frac{\rho_0}{A} &
    1
\end{bmatrix}
\underbrace{
\begin{bmatrix}
    A \rho_0 v \\\
    A \frac{c_0^2}{\rho_0} \rho \\\
    \frac{1}{2} \left(\rho_0 v^2 + \frac{c_0^2 \rho^2}{\rho_0}\right)
\end{bmatrix}
}_{\nabla \mathcal H}.
$$

Note that the pressure here includes an unexpected term related to the last element of the Hamiltonian gradient. This follows from the linearization process: it is a second order perturbation to the received pressure by the walls, non-physical but structurally obtained to preserve the power balance. 

Boundary conditions at $x=0$ and $x=l_0$ are considered to be either of the from of a volume flow or pressure control (with the other as the power conjugated output).



### Separation of small area variations

The above formulation results in linearity of the internal dynamics but still makes appear a nonlinear term corresponding to $\delta_A \mathcal H$ (a.k.a. the spurious pressure term from before). In view of numerical simulation, this is not relevant while $A$ is treated as a parameter externally forced with no dynamics associated to it. For the vocal tract, it is indeed practical to consider $A$ as such and directly set its value to corresponds to articulations. However, vibrations of the surrounding tissues around the vocal tract are known to have a significant impact on the effective damping at low frequencies. To represent this effect, the acoustic propagation equation must be coupled to a mechanical model of the tissues. The coupled system would then be internally nonlinear and require additional care in the solver design. In order to prevent that, the area is decomposed as
$$
 A(x, t) = A_0(x, t) + \widetilde{A}(x, t),
$$
where $A_0$ is a controlled geometry variation representing articulations and $\widetilde{A}$ is a small amplitude area variation due to wall vibrations (he wall mechanical model being define independently). After replacement of $A$ in \eqref{base_webster}, the dynamics are further linearized with respect to $\widetilde{A}$, yielding
$$
\begin{equation}
    \begin{bmatrix}
        \partial_t v \\\
        \partial_t \rho \\\
        \partial_t A_0
    \end{bmatrix}
    =
    \begin{bmatrix}
        0 &-\partial_x\left(\frac{\bullet}{A_0}\right) & 0 \\\
        -\frac{1}{A_0} \partial_x & 0 & 0 \\\
        0 & 0 & 0
    \end{bmatrix}
    \underbrace{
    \begin{bmatrix}
        A_0 \rho_0 v \\\
        A_0 \frac{c_0^2}{\rho_0} \rho \\\
        \frac{1}{2} \left(\rho_0 v^2 + \frac{c_0^2 \rho^2}{\rho_0}\right)
    \end{bmatrix}
    }_{\nabla \mathcal H}
    +
    \begin{bmatrix}
        0 & 0\\\
        - \frac{\rho_0}{A_0} & - \frac{\rho_0}{A_0}\\\
        1 & 0
    \end{bmatrix}
    \underbrace{
    \begin{bmatrix}
        \partial_t A_0 \\\
        \partial_t \widetilde{A}
    \end{bmatrix}}_u,
    \label{base_webster_lin}
\end{equation}
$$
associated to state $x = [v, \rho, A_0]^\intercal$ and Hamiltonian 

$$
    \mathcal H (x) = \frac{1}{2} \left(\rho_0 \vert v \vert^2_{A_0} + \frac{c_0^2}{\rho_0} \vert \rho \vert_{A_0}^2\right).
$$

<details><summary>Power balance (EDP)</summary>

The collapsed section's content goes here

</details>

## Spatial discretization (finite differences)

Finite differences on staggered grids are used for spatial discretization (see e.g. Trenchant 2018[@Tre18]). The presentation is made here for volume flow input on both end of the tube. 

$\rho$, $A_0$ and $\widetilde A$ are discretized on a primal grid with constant spatial step $h$ and $N$ edges. $v$ is discretized on the dual (interleaved) grid with $N-1$ edges. The finite difference matrix 
$$
		\mathbf D^- = \frac{1}{h}
		\begin{bmatrix}
			1 & 0 & 0  \\\
			-1 & 1 & 0 \\\
			0 & \ddots & \ddots \\\
			0 & 0 &  -1\\
		\end{bmatrix} \in \mathbb R^{N, N-1},
$$
builds on the dual grid the derivative of a quantity defined on the primal grid. $\mathbf D^+ = - (\mathbf D^-)^\intercal$ does the opposite (derivative from dual to primal). Considering piecewise constant values of the field variables on segments of their respective grids, a discrete equivalent to \eqref{base_webster_lin} is obtained as (in the following equation, $\boldsymbol{A_p}$ and $\boldsymbol{A_p}$ are vectors but are used as diagonal matrices in multiplicative expressions to simplify notations)

$$
\begin{align}
    \begin{bmatrix}
        \dot{\boldsymbol{v}} \\\
        \dot{\boldsymbol{\rho}} \\\
        \dot{\boldsymbol{A_p}}
    \end{bmatrix}
    =&
    \begin{bmatrix}
        0 &- \mathbf D^+ \boldsymbol{A_p}^{-1} & 0 \\\
        -\boldsymbol{A_p}^{-1}  \mathbf D^- & 0 & 0 \\\
        0 & 0 & 0
    \end{bmatrix}
    \underbrace{
    \begin{bmatrix}
        \boldsymbol{A_d} \rho_0 \boldsymbol{v} \\\
        \boldsymbol{A_p} \frac{c_0^2}{\rho_0} \boldsymbol{\rho} \\\
        \frac{1}{2} \left(\rho_0 \boldsymbol{v}^2 + \frac{c_0^2 \boldsymbol{\rho}^2}{\rho_0}\right)
    \end{bmatrix}
    }_{\nabla H}\nonumber\\\
    &+
    \begin{bmatrix}
        0 & 0 & 0 & 0\\\
        - \rho_0 h \boldsymbol{A_p}^{-1} & - \rho_0 h \boldsymbol{A_p}^{-1} & \rho_0 \boldsymbol{A_p}^{-1}\vert_{i=0}& -\rho_0 \boldsymbol{A_p}^{-1}\vert_{i=N-1}\\\
        1 & 0 & 0 & 0
    \end{bmatrix}
    \underbrace{
    \begin{bmatrix}
        \partial_t \boldsymbol{A_0} \\\
        \partial_t \widetilde{\boldsymbol{A}}\\\
        Q_{l} \\\
        Q_{r}
    \end{bmatrix}}_{\boldsymbol{u}},
    \label{semi_discrete_webster}
\end{align}
$$

with state $\mathbf x = [\boldsymbol{v}, \boldsymbol{\rho}, \boldsymbol{A}_p]^\intercal$ and Hamiltonian

$$
    H (\mathbf x) = \frac{1}{2} \left(\rho_0 \vert \boldsymbol{v} \vert^2_{\boldsymbol{A_d}} + \frac{c_0^2}{\rho_0} \vert \boldsymbol{\rho} \vert_{\boldsymbol{A_p}}^2\right).
$$

$\boldsymbol{A_p}$ and $\boldsymbol{A_d}$ are linked by an averaging operator to be specified later, as one choice is better with the chosen time discretization scheme. The notation $\rho_0 \boldsymbol{A_p}^{-1}\vert_{i=0}$ indicates that all values expect the one at index 0 are null (multiplication with a discrete Dirac function).

The conjugated output vector, including both the distributed pressure applied on the walls and the boundaries pressure writes
$$
    \boldsymbol{y} = 
    \begin{bmatrix}
        \boldsymbol{P_w} \\\
        \widetilde{\boldsymbol{P_w}} \\\
        - P_l \\\
        P_r
    \end{bmatrix}
    =
    \begin{bmatrix}
        0 &  \rho_0 h \boldsymbol{A_p}^{-1} & -1 \\\ 
        0 &  \rho_0 h \boldsymbol{A_p}^{-1} & 0 \\\ 
        0 & - \rho_0 \boldsymbol{A_p}^{-1}\vert_{i=0} & 0  \\\ 
        0 & \rho_0 \boldsymbol{A_p}^{-1}\vert_{i=N-1} &0 
    \end{bmatrix}
    \begin{bmatrix}
        \boldsymbol{A_d} \rho_0 \boldsymbol{v} \\\
        \boldsymbol{A_p} \frac{c_0^2}{\rho_0} \boldsymbol{\rho} \\\
        \frac{1}{2} \left(\rho_0 \boldsymbol{v}^2 + \frac{c_0^2 \boldsymbol{\rho}^2}{\rho_0}\right)
    \end{bmatrix}.
$$


<details><summary>Power balance (ODE continuous time)</summary>

The collapsed section's content goes here

</details>

## Time discretization (Störmer-Verlet)

Time discretization is made using a Störmer-Verlet scheme. It writes:
$$
\begin{align}
    \delta_{t+} \boldsymbol{v}^{n-\frac{1}{2}} &= - \mathbf D^+ \frac{c_0^2}{\rho_0} \boldsymbol{\rho}^n \\\
    \delta_{t+} \boldsymbol{\rho}^n &= - \boldsymbol{A_p}^{-1} \boldsymbol{D}^{-} \boldsymbol{A_d} \rho_0 \boldsymbol{v^{n+\frac{1}{2}}} - \rho_0 \boldsymbol{A_p}^{-1} \left(h \partial_t \boldsymbol{A_0} + h \partial_t \widetilde{\boldsymbol{A}} - Q_l \vert_{i=0} + Q_r\vert_{i=N-1}\right)^{n+\frac{1}{2}} \\\
    \delta_{t+} \boldsymbol{A_0}^{n-\frac{1}{2}} &= (\partial_t \boldsymbol{A_0})^{n+\frac{1}{2}}
\end{align}
$$

<details><summary>Power balance (ODE discrete time)</summary>

The collapsed section's content goes here

</details>