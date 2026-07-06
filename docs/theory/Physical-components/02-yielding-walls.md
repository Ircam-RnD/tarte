---
title: Distributed surface impedance (yielding walls)
---

In the vocal tract, surrounding tissues are not perfectly rigid. Due to the vibro-acoustics coupling with the fluid wave propagation, they may vibrate during phonation.Elasticity in the tissues then act as a significant source of dissipation in the vocal tract, with direct effect on the frequency and quality factor of its first resonant frequency (see e.g. Birkholz 2022 [@birkholz2022acoustic]).

## Surface impedance

A common approach to model the dynamics of the vocal tract tissues consists of assuming a linear and local behavior, which is represented in the frequency domain by a mechanical impedance per-unit-area $Z_{w}(\omega) = \frac{\sigma(\omega)}{v(\omega)}$ (with $\sigma(\omega)$ the normal stress at the boundary and $v(\omega)$ the velocity). In Ishizaka 1975 [@ishizaka1975direct], the authors present measurements of neck and cheek tissues impedances. They show that in the low frequency domain, the impedance is locally well represented by a (per-unit-area) mass-spring-damper system, thus writing in the Fourier domain
$$
  Z_{w}(\omega) = R_{w} + j \omega M_{w} + \frac{K_{w}}{j \omega},
$$
where per-unit-area parameters are the mass $M_{w}$ (kg$\cdot$m$^{-2}$), stiffness $K_{w}$ (N$\cdot$m$^{-3}$) and dissipation coefficient $R_{w}$ (kg$\cdot$s$^{-1}$m$^{-2}$). The distributed state $\boldsymbol{x} = [q_w, p_w \equiv M_w \dot q_w]$ composed of the tissue radial deformation $q_w$ and associated (per-area) momentum $p_w$ is defined on the fluid-structure interface $\partial \Omega_s$.
<!-- \begin{equation}
  \partial\Omega_s = \{\xi = (x, y, z) \in \partial S(x, t)  \}.
  \label{eq:fs_interface}
\end{equation} -->
It allows to compute the energy stored by the vocal tract walls as
$$
  H(q_w, p_w) = \int_{\partial \Omega_s} \frac{1}{2} \left( \frac{p_w^2}{M_w} + K_w q_w^2\right) dS,
$$
with associated dynamics
$$
  \begin{bmatrix}
    \partial_t q_w \\\
    \partial_t p_w \\\
    y = - \dot q_w
  \end{bmatrix}
  =
  \begin{bmatrix}
    0 & 1 & 0 \\\
    -1 & -R_w & 1 \\\
    0 & -1 & 0
  \end{bmatrix}
  \begin{bmatrix}
    K_w q_w \\\
    \frac{p_w}{M_w}\\\
    u = P_f
  \end{bmatrix},
$$
where $u= P_f$ is the acoustic pressure applied by the fluid to the wall, and $ y = - \dot q_w $ is the resulting velocity of the wall. 

## 1D reduction
The reduction of the two-dimensional representation (as $\partial \Omega_s$ is a surface) to a 1D representation is performed by further assuming that:

1. Both the physical parameters ($M_w, K_w, R_w$) and the radial displacement $q_w$ are homogeneous along the cross-sections perimeters (i.e. functions of $x$ only).
2. Cross-section area time-variations are small: $A(x, t) = A_0(x) + \widetilde{A}(x, t)$.

In this case, the energy rewrites 
$$
\begin{align}
  H(q_w, p_w) &= \int_{0}^{l_0}\int_{\partial S(x)} \frac{1}{2} \left( \frac{p_w^2}{M_w} + K_w q_w^2\right) d\Gamma(x)\, dx,\nonumber\\\
  &= \int_{0}^{l_0} \frac{\Gamma_0}{2} \left( \frac{p_w^2}{M_w} + K_w q_w^2\right) \, dx,
\end{align}
$$
where $\Gamma_0(x)$ is the cross-section perimeter. The associated governing equations are
$$
\begin{equation}
  \begin{bmatrix}
    \partial_t q_w \\\
    \partial_t p_w \\\
    y_w = - \Gamma_0 \dot q_w
  \end{bmatrix}
  =
  \begin{bmatrix}
    0 & \frac{1}{\Gamma_0} & 0 \\\
    -\frac{1}{\Gamma_0} & -\frac{R_w}{\Gamma_0} & 1 \\\
    0 & -1 & 0
  \end{bmatrix}
  \begin{bmatrix}
    \Gamma_0 K_w q_w \\\
    \Gamma_0 \frac{p_w}{M_w} \\\
    u_w = P_f
  \end{bmatrix},
  \label{eq:wall_impedance_pH}
\end{equation}
$$
where the only difference with the previous formulation is that the output of the system is now the 1D distributed (small) area variation $y = - \Gamma_0 \dot q_w = - \partial_t \widetilde A$.  

## Spatial discretization

The discrete equivalent to \eqref{eq:wall_impedance_pH} is obtained in a straightforward manner by making a finite difference approximation of the distributed variables $q_w$ and $p_w$ on the same grid with step-size $l_{d0}$ and replacing them by their vector equivalent $\boldsymbol q_w$ and $\boldsymbol p_w$. 
Similarly, $\boldsymbol \Gamma_0$, $\boldsymbol K_w$ and $\boldsymbol M_w$ are diagonal matrices discrete equivalent to $\Gamma_0$, $K_w$ and $M_w$.
The pHs writes

$$
\begin{equation}
  \begin{bmatrix}
    \frac{d}{dt} \boldsymbol q_w \\\
    \frac{d}{dt} \boldsymbol p_w \\\
    \boldsymbol y_w = - \boldsymbol \Gamma_0 \dot{\boldsymbol q}_w
  \end{bmatrix}
  =
  \frac{1}{l_{d0}}
  \begin{bmatrix}
    0 & \boldsymbol \Gamma_0^{-1} & 0 \\\
    \boldsymbol \Gamma_0^{-1} & -\boldsymbol R_w\boldsymbol \Gamma_0^{-1} & \mathbb I \\\
    0 & -\mathbb I & 0
  \end{bmatrix}
  \begin{bmatrix}
    l_{d0} \boldsymbol \Gamma_0 \boldsymbol K_w \boldsymbol q_w \\\
    l_{d0} \boldsymbol \Gamma_0 \dot{\boldsymbol q}_w \\\
    \boldsymbol u_w = l_{d0} (\boldsymbol P_f)
  \end{bmatrix},
  \label{eq:wall_impedance_discrete}
\end{equation}
$$

and corresponds to the discrete Hamiltonian 

$$
\begin{equation}
  H(\boldsymbol x = [\boldsymbol q_w, \boldsymbol p_w]^\intercal) = 
  \frac{l_{d0}}{2} \left(
    \boldsymbol q_w^\intercal \boldsymbol \Gamma_0 \boldsymbol K_w  \boldsymbol q_w
    + 
    \boldsymbol p_w^\intercal \boldsymbol \Gamma_0 \boldsymbol M_w^{-1} \boldsymbol p_w
  \right).
  \label{eq:wall_Ham}
\end{equation}
$$


The estimation of parameters $M_w, K_w$ and $R_w$ is not straightforward. Indeed, in-vivo measurements are hard to perform (as an example, in Ishizaka 1975, measurements are made on tissues from the outside of the neck). Moreover, parameter values are expected to be spatially dependent to correspond to e.g. different tissues of the tongue and palate. This motivates the use of optimization procedures to estimate the wall parameter values, as in e.g. Fleischer 2015[@fleischer2015formant]. 
