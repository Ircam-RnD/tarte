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

The tube length is considered here to be constant. This hypothesis can be somewhat relaxed for the vocal tract by modulating the length in a non-power preserving way, but still stable under slow variations. The axial variable is named $\xi$.

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
        0 &-\partial_{\xi}\left(\frac{\bullet}{A}\right) & 0 \\\
        -\frac{1}{A} \partial_{\xi} & 0 & 0 \\\
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

Boundary conditions at $\xi=0$ and $\xi=l_0$ are considered to be either of the from of a volume flow or pressure control (with the other as the power conjugated output).



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
        0 &-\partial_{\xi}\left(\frac{\bullet}{A_0}\right) & 0 \\\
        -\frac{1}{A_0} \partial_{\xi} & 0 & 0 \\\
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

The power balance is obtained directly by using the formula

$$
\begin{align*}
    \frac{d}{dt} \mathcal H &= \int_{\xi=0}^{\xi = l_0} \nabla \mathcal H^\intercal \partial_t \boldsymbol{x} \; d\xi, \\\
    &= \int_{\xi=0}^{\xi = l_0} - \partial_{\xi} \left( A_0 v c_0^2 \rho\right) - c_0^2 \rho (\partial_t A_0 + \partial_t \widetilde{A})+ 
    \frac{1}{2} \left(\rho_0 v^2 + \frac{c_0^2 \rho^2}{\rho_0}\right) \partial_t A_0 \; d\xi,\\\
    &= -[QP]_{\xi = 0}^{\xi=l_0} + 
    \int_{\xi=0}^{\xi = l_0} - c_0^2 \rho (\partial_t A_0 + \partial_t \widetilde{A})+ 
    \frac{1}{2} \left(\rho_0 v^2 + \frac{c_0^2 \rho^2}{\rho_0}\right) \partial_t A_0 \; d\xi.

\end{align*}
$$

As the interconnection operators is formally skew-symmetric, no dissipation terms are present. The energy exchanged with the walls, due to variations of $A_0$ and $\widetilde A$ includes the integral of the product of area variations and pressure along the tube. Is also includes the spurious second order pressure term. Energy exchanged at the tube boundaries writes as the product of volume flow and pressure at $\xi = 0$ and $\xi = l_0$.

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
    \frac{1}{h}
    \begin{bmatrix}
        0 &- \mathbf D^+ \boldsymbol{A_p}^{-1} & 0 \\\
        -\boldsymbol{A_p}^{-1}  \mathbf D^- & 0 & 0 \\\
        0 & 0 & 0
    \end{bmatrix}
    \underbrace{
    \begin{bmatrix}
        h\boldsymbol{A_d} \rho_0 \boldsymbol{v} \\\
        h\boldsymbol{A_p} \frac{c_0^2}{\rho_0} \boldsymbol{\rho} \\\
        h\frac{1}{2} \left(\rho_0 \boldsymbol{v}^2 + \frac{c_0^2 \boldsymbol{\rho}^2}{\rho_0}\right)
    \end{bmatrix}
    }_{\nabla H}\nonumber\\\
    &+
    \begin{bmatrix}
        0 & 0 & 0 & 0\\\
        - \rho_0  \boldsymbol{A_p}^{-1} & - \rho_0  \boldsymbol{A_p}^{-1} & \frac{1}{h} \rho_0 \boldsymbol{A_p}^{-1}\vert_{i=0}& - \frac{1}{h}\rho_0 \boldsymbol{A_p}^{-1}\vert_{i=N-1}\\\
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
    H (\mathbf x) = \frac{1}{2} h \left(\rho_0 \vert \boldsymbol{v} \vert^2_{\boldsymbol{A_d}} + \frac{c_0^2}{\rho_0} \vert \boldsymbol{\rho} \vert_{\boldsymbol{A_p}}^2\right).
$$

$\boldsymbol{A_p}$ and $\boldsymbol{A_d}$ are linked by an averaging operator to be specified later, as one choice is better with the chosen time discretization scheme. The notation $\rho_0 \boldsymbol{A_p}^{-1}\vert_{i=0}$ indicates that all values expect the one at index 0 are null (multiplication with a discrete Dirac function).

The conjugated output vector, including both the linear force density applied on the walls and the boundaries pressure writes
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
        0 &  \rho_0 \boldsymbol{A_p}^{-1} & -1 \\\ 
        0 &  \rho_0 \boldsymbol{A_p}^{-1} & 0 \\\ 
        0 & - \frac{1}{h} \rho_0 \boldsymbol{A_p}^{-1}\vert_{i=0} & 0  \\\ 
        0 & \frac{1}{h} \rho_0 \boldsymbol{A_p}^{-1}\vert_{i=N-1} &0 
    \end{bmatrix}
    \begin{bmatrix}
        h \boldsymbol{A_d} \rho_0 \boldsymbol{v} \\\
        h \boldsymbol{A_p} \frac{c_0^2}{\rho_0} \boldsymbol{\rho} \\\
        h \frac{1}{2} \left(\rho_0 \boldsymbol{v}^2 + \frac{c_0^2 \boldsymbol{\rho}^2}{\rho_0}\right)
    \end{bmatrix}.
$$


<details><summary>Power balance (ODE continuous time)</summary>

In the finite dimensional case, the power balance is easier to compute. Direct application of port-Hamiltonian theory yields:
$$
\frac{d}{dt} H  + \boldsymbol{y}^\intercal \boldsymbol{u} = 0.
$$
Elements of $\boldymbol u$ and $\boldymbol y$ are direct counterparts to the PDE case.

</details>

## Time discretization (Störmer-Verlet)

Time discretization is made using a Störmer-Verlet scheme. It writes:
$$
\begin{align}
    \delta_{t+} \boldsymbol{v}^{n-\frac{1}{2}} &= - \mathbf D^+ \frac{c_0^2}{\rho_0} \boldsymbol{\rho}^n \label{eq:update_v} \\\
    \delta_{t+} \boldsymbol{\rho}^n &= - \boldsymbol{A_p}^{-1} \boldsymbol{D}^{-} \boldsymbol{A_d} \rho_0 \boldsymbol{v^{n+\frac{1}{2}}} \nonumber \\\
    & \quad  - \rho_0 \boldsymbol{A_p}^{-1} \left(\partial_t \boldsymbol{A_0} + \partial_t \widetilde{\boldsymbol{A}} - \frac{1}{h}Q_l \vert_{i=0} + \frac{1}{h}Q_r\vert_{i=N-1}\right)^{n+\frac{1}{2}} \label{eq:update_rho}\\\
    \delta_{t+} \boldsymbol{A_0}^{n-\frac{1}{2}} &= (\partial_t \boldsymbol{A_0})^{n+\frac{1}{2}}
\end{align}
$$

<details><summary>Power balance (ODE discrete time)</summary>

The fully discrete power balance is slightly more challenging to get to. Indeed, the chosen scheme does not preserve a discrete equivalent to the continuous time energy, but rather a modified version of it given by:

$$
\begin{align}
    E^{n} &= K^{n} + P^{n} 
    -
    \frac{1}{4} \frac{h \rho_0}{2}\vert \vert \boldsymbol{v}^{n+\frac{1}{2}} - \boldsymbol{v}^{n-\frac{1}{2}}\vert \vert^2_{\boldsymbol{A_d}}, \nonumber \\\
    &\overset{\eqref{eq:update_v}}{=} K^{n} + P^{n} 
    -
    \frac{dt^2}{4} \frac{h c_0^4}{2 \rho_0 } \vert \vert D^+(\boldsymbol{\rho}^{n })\vert \vert^2_{\boldsymbol{A_d}}, \label{eq:energy_sv}
\end{align}
$$

with 

$$
\begin{align*}
    K^{n} = \frac{1}{2} h \rho_0 \vert\vert \mu_{t+}\boldsymbol{v}^{n-\frac{1}{2}}\vert\vert^2_{\boldsymbol{A_d}},
		%
		\quad P^{n} = \frac{1}{2} h \frac{c_0^2}{\rho_0}\vert\vert \boldsymbol{\rho}^{n}\vert\vert^2_{\boldsymbol{A_p}}.
\end{align*}
$$

In addition to that, the exchanged power writes 
$$
 P_{\rm ext}^{n+\frac{1}{2}} = h c_0^2 \mu_{t+} \boldsymbol{\rho}^n
 \left(\partial_t \boldsymbol{A_0} + \partial_t \widetilde{\boldsymbol{A}} - \frac{1}{h}Q_l \vert_{i=0} + \frac{1}{h}Q_r\vert_{i=N-1}\right)^{n+\frac{1}{2}},
$$

such that the discrete power balance is

$$
    E^{n+1} - E^n + dt P_{\rm ext}^{n+\frac{1}{2}} = 0.
$$

<details><summary>Proof</summary>
Take the scalar product of \eqref{eq:update_rho} and 

$$h\boldsymbol{A}_p\frac{c_0^2}{\rho_0} \mu_{t+} \boldsymbol{\rho}^n \approx (dt \; \partial_\boldsymbol{\rho} H)^{n+\frac{1}{2}}$$

resulting in 

$$
\begin{align*}
    P^{n+1} - P^{n} &= - hdt c_0^2 \left(\mu_{t+} \boldsymbol{\rho}^{n}, D^- \boldsymbol{A_d} \boldsymbol{v}_x^{n+\frac{1}{2}} \right)
    - dt P_{ext}^{n+\frac{1}{2}}
    , \\
    &= h dt c_0^2 \left(D^+ \mu_{t+} \boldsymbol{\rho}^{n}, \boldsymbol{A_d}  \boldsymbol{v}_x^{n+\frac{1}{2}}\right) - dt P_{ext}^{n+\frac{1}{2}},\\
    & = - \frac{h \rho_0}{2} \left(\boldsymbol{v}_x^{n+\frac{3}{2}} - \boldsymbol{v}_x^{n-\frac{1}{2}}, \boldsymbol{A_d}  \boldsymbol{v}_x^{n+\frac{1}{2}}\right) - dt P_{ext}^{n+\frac{1}{2}},\\
    &= - K^{n+1} + K^{n-1} + \frac{1}{4} \frac{h \rho_0}{2} \left(\vert\vert\boldsymbol{v}_x^{n+\frac{3}{2}} - \boldsymbol{v}_x^{n+\frac{1}{2}}\vert\vert_{\boldsymbol{A_d}} -\vert\vert\boldsymbol{v}_x^{n+\frac{1}{2}} - \boldsymbol{v}_x^{n-\frac{1}{2}}\vert\vert_{\boldsymbol{A_d}} \right) - dt P_{ext}^{n+\frac{1}{2}}
\end{align*}
$$

</details>

For the algorithm to be stable, the preserved energy must be positive definite. In order to ensure that, a bound on the third term of \eqref{eq:energy_sv} has to be derived:

$$
\begin{align*}
    \vert \vert D^+(\boldsymbol{\rho}^{n })\vert \vert^2_{\boldsymbol{A_d}} &= \frac{1}{h^2} \vert \vert \boldsymbol{\rho}^{n }_{0...N-1} - \boldsymbol{\rho}^{n }_{1...N}\vert \vert^2_{\boldsymbol{A_d}} \\\
    & \leq \frac{2}{h^2} \left(\vert \vert \boldsymbol{\rho}^{n }_{0...N-1} \vert\vert_{\boldsymbol{A_d}}^2 +  \vert \vert \boldsymbol{\rho}^{n }_{1...N} \vert\vert_{\boldsymbol{A_d}}^2\right)\\\
    & \leq \frac{4}{h^2} \left(\vert \vert \boldsymbol{\rho}^{n } \vert\vert_{\mu_{\xi_-}\boldsymbol{A_d}}^2\right)
\end{align*}
$$

The stability condition reduces to
$$
\begin{align*}
    P^n - \frac{dt^2}{4} \frac{h c_0^4}{2 \rho_0 } \frac{4}{h^2} \left(\vert \vert \boldsymbol{\rho}^{n } \vert\vert_{\mu_{\xi_-}\boldsymbol{A_d}}^2\right) &\geq 0, \\\
    \vert\vert \boldsymbol{\rho}^{n}\vert\vert^2_{\boldsymbol{A_p}} - \frac{dt^2}{h^2} c_0^2 \vert \vert \boldsymbol{\rho}^{n } \vert\vert_{\mu_{\xi_-}\boldsymbol{A_d}}^2 &\geq 0, \\\
    dt^2 c_0^2 \; \frac{\vert \vert \boldsymbol{\rho}^{n } \vert\vert_{\mu_{\xi_-}\boldsymbol{A_d}}^2}{\vert\vert \boldsymbol{\rho}^{n}\vert\vert^2_{\boldsymbol{A_p}}} &\leq h^2, \\\
    dt^2 c_0^2 \rm{max}(\frac{\mu_{\xi_-}\boldsymbol{A_d}}{\boldsymbol{A_p}}) \leq h^2,
\end{align*}
$$
meaning that under the suitable choice 

$$\boldsymbol A_p = \mu_{\xi_-} \boldsymbol A_d,$$

a standard-like CFL condition is recovered as

$$
    h^2 \geq dt^2 c_0^2.
$$

At the boundaries, 

$$
 \boldsymbol{A_p},
$$

can be of any value at least half the value of the boundaries of

$$
    \boldsymbol{A_d}.
$$


</details>