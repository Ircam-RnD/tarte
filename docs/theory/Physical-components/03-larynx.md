---
title: 'Larynx (vocal folds + flow) models'
---

Models of the larynx include two components and their coupling:

- The vocal folds mechanical model,
- The fluid flow model.

These two components are responsible for self-oscillation, and ultimately of sound production (at least for `voiced' sounds). 

The vocal folds are a complex assembly of tissues with various mechanical properties, which can to a certain extent be described by constitutive laws. Their dynamics can therefore be derived from the fundamental laws of mechanics, resulting in three dimensional PDEs models. From these, discretization can be made using high-fidelity methods like the finite element method, or using lumped parameter modeling. The two approaches are quite distinct in their objectives and computational complexity. In this work, low dimensional models are used to keep compute time low. As a rule of thumb, the real-time barrier should become a problem for vocal folds models with ~100 degrees of freedom with the current simulation method and code. The low dimensional models are known to be able to reproduce several qualitative properties of voice production. Large-scale parameter studies are also easier to perform than with larger models, due to the reduced computational cost. However, one needs to be extra cautious when interpreting simulation results as validation is generally hard, especially when comparing with in-vivo data.

Fluid flow between the vocal folds is also complex to represent accurately due to the time-varying domain and to the presence of turbulences. In this work, the classical reduced-order way of representing the glottal flow is used. It consists in a one-dimensional model, assuming incompressible, quasi-stationary flow that separates from the folds to form a jet at the exit of the larynx. 

## General model class

### Vocal folds

The dynamics are assumed to be written in terms of a state $\boldsymbol x_{\rm f} = [\boldsymbol p_{\rm f}, \boldsymbol q_{\rm f}]^\intercal$, including vectors of generalized coordinates $\boldsymbol q_{\rm f}$ and momenta $\boldsymbol p_{\rm f}$. The Hamiltonian (total energy) writes as a function of the state

$$
\begin{equation}
  H (\boldsymbol x_{\rm f}) = \frac{1}{2} \boldsymbol p_{\rm f}^\intercal \boldsymbol M^{-1} \boldsymbol p_{\rm f} + \underbrace{\frac{1}{2} \boldsymbol q_{\rm f}^\intercal \boldsymbol K \boldsymbol q_{\rm f}  + E_{\rm nl}(\boldsymbol q_{\rm f})}_{V(\boldsymbol q_{\rm f})},
  \label{eq:Hamiltonian}
\end{equation}
$$

with diagonal mass matrix $\boldsymbol M$, and symmetric stiffness matrix $\boldsymbol K$, both assumed positive definite. The potential energy $V(\boldsymbol q_{\rm f})$ is the sum of a quadratic term and of a general positive function of the generalized coordinates $E_{\rm nl}(\boldsymbol q_{\rm f})$, with $E_{\rm nl}(\boldsymbol 0) = 0$.

The vocal folds dynamics, including dissipation, evolve according to

$$
\begin{equation}
  \underbrace{
    \begin{bmatrix}
      \dot{\boldsymbol q_{\rm f}} \\
      \dot{\boldsymbol p_{\rm f}} \\
      \boldsymbol y
    \end{bmatrix}}_{\boldsymbol f}
  =
  \begin{bmatrix}
    \boldsymbol 0              & \boldsymbol I          & \boldsymbol 0           \\
    -\boldsymbol I & - \boldsymbol R_{0}      & \boldsymbol I           \\
    \boldsymbol 0              & -\boldsymbol I & \boldsymbol{0}
  \end{bmatrix}
  \underbrace{
    \begin{bmatrix}
      \nabla V           \\
      \boldsymbol M^{-1} \boldsymbol p_{\rm f} \\
      \boldsymbol u
    \end{bmatrix}}_{\boldsymbol e},
  \label{eq:dynamics_vf}
\end{equation}
$$

where $\boldsymbol u = \pm \boldsymbol{F}_{\rm f}$ is an external forcing and the associated output $\boldsymbol y = \pm -\boldsymbol{v}_{\rm f}$ is the velocity corresponding to generalized momentum $\boldsymbol p_{\rm f}$. 

<details><summary>Sign of input-output pair  </summary>

Inputs and output signs depend on the chosen convention for the vocal folds displacement. For the model to be coherent, positive forces must open the glottis. With the above sign convention, this means that:
if higher values of q_f result in bigger opening, then a plus must replace the above plus/minus signs, else, a minus must replace the above plus/minus signs. This propagates in the construction of the assembled model below

</details>


By construction of the flow $\boldsymbol f$ and power-conjugated effort $\boldsymbol e$ vectors, one gets the power balance, using the chain rule and \eqref{eq:dynamics_vf}:

$$
\begin{align}
      \frac{d}{dt} H + \boldsymbol u^\intercal \boldsymbol y \overset{(\boldsymbol e^\intercal \boldsymbol f)}{=}
      -(\boldsymbol M^{-1} \boldsymbol p_{\rm f})^\intercal 
      \boldsymbol{R_0}
        \boldsymbol M^{-1} \boldsymbol p_{\rm f}.
  \label{eq:power_balance}
\end{align}
$$

In the absence of an input ($\boldsymbol u = 0$), the energy is preserved (for $\boldsymbol R_0 = \boldsymbol R_{\rm uu} = 0$) or decreases owing to dissipation.

### Glottal flow

The flow description is composed of two equations. The first relates the pressure drop $\Delta P = P_{-} - P_{+}$ (where $P_-$ is the sub-glottal pressure and $P_+$ the supra-glottal pressure) to a resulting mean flow $Q_{\rm g}$ in the glottis

$$
  \begin{equation}
    Q_{\rm g} = A_{\rm min}(\boldsymbol q_{\rm f}) {\rm sign}(\Delta P) \sqrt{\frac{2 \vert \Delta P\vert}{\rho_0} },
    \label{eq:Bernoulli_vf}
  \end{equation}
$$

where $A_{\rm min} (\boldsymbol q_{\rm f})$ is the effective open area for the flow. The second equation describes the forces applied by the fluid pressure to the vocal folds. These forces can be written as linear combinations of the sub-glottal and supra-glottal pressures:

$$
\begin{equation}
\boldsymbol F_{\rm f} =  {\boldsymbol S}_{\rm f-} (\boldsymbol q_{\rm f}) P_- +  {\boldsymbol S}_{\rm f+} (\boldsymbol q_{\rm f}) P_+,
\label{eq:forces_vf}
\end{equation}
$$

where ${\boldsymbol S}_{\rm f-} (\boldsymbol q_{\rm f})$ and ${\boldsymbol S}_{\rm f+} (\boldsymbol q_{\rm f})$ are vectors of effective surfaces depending on the current geometry of the glottis. 


In most lumped models of the literature, the flow pumped by the vocal folds displacement is not included, as it is mostly negligible in front of the mean flow $Q_g$ in the glottis owing to the low range of oscillation frequencies. However, this choice results in non-power preserving and non-mass preserving models. Hence, it is not directly compatible with the port-Hamiltonian framework. In this work, we propose to minimally modify the model to include power conjugated flows induced by the vocal fold velocity vector $\boldsymbol v_{\rm f} = \boldsymbol M_{\rm f}^{-1} \boldsymbol p_{\rm f}$ as

$$
\begin{equation}
Q_{-} =  -Q_g -{\boldsymbol S}_{\rm f-}^\intercal (\boldsymbol q_{\rm f}) \boldsymbol v_{\rm f}, \quad Q_{+} =  Q_g - {\boldsymbol S}_{\rm f+}^\intercal (\boldsymbol q_{\rm f}) \boldsymbol v_{\rm f},
\label{eq:fold_flow}
\end{equation}
$$

with outgoing lower and upper flows $Q_{-}, Q_{+}$.
This writing effectively restores both the mass and power balances of the model.

The flow component may be written in a compact form as

$$
\begin{equation}
      \begin{bmatrix}
        y_0 = Q_{-} \\
        y_1 = Q_{+} \\
        \boldsymbol y_{2, 3} = \boldsymbol F_{\rm f}
      \end{bmatrix}
      =
      \begin{bmatrix}
        -R_{\rm g}         & R_{\rm g}          & -{\boldsymbol S}^{\intercal}_{\rm f-} \\
        R_{\rm g}          & -R_{\rm g}         & -{\boldsymbol S}^{\intercal}_{\rm f+} \\
        {\boldsymbol S}_{\rm f-} & {\boldsymbol S}_{\rm f+} & \boldsymbol 0
      \end{bmatrix}
      \begin{bmatrix}
        u_0 = P_{-} \\
        u_1 = P_{+} \\
        \boldsymbol u_{2, 3} = \boldsymbol v_{\rm f}
      \end{bmatrix},
  \label{eq:glottis}
\end{equation}
$$

with

$$
\begin{equation}
  R_{\rm g}(\Delta P, \boldsymbol q_{\rm f}) = \frac{Q_{\rm g}(\Delta P, A_{\rm min}(\boldsymbol q_{\rm f}))}{\Delta P}.
\end{equation}
$$

Note that the previous system has the vocal fold configuration $\boldsymbol q_{\rm f}$ as a parameter, making in non-autonomous. When connecting to the vocal fold model (see next section) $\boldsymbol q_{\rm f}$ is recovered as a state of the assembly, therefore closing the system. The first presentation as two separated sub-systems 

### Assembly

Assembly of the vocal folds and fluid flow models is done through the power ports by imposing forces and velocities equalities ($\boldsymbol{F}_{\rm f} = \boldsymbol{F}_{\rm f}$ and $\boldsymbol{v}_{\rm f} = \boldsymbol{v}_{\rm f}$). The assembled system then writes

$$
\begin{equation}
  \underbrace{
    \begin{bmatrix}
      \dot{\boldsymbol q_{\rm f}} \\
      \dot{\boldsymbol p_{\rm f}} \\
      \boldsymbol y
    \end{bmatrix}}_{\boldsymbol f}
  =
  \begin{bmatrix}
    \boldsymbol 0              & \boldsymbol I          & \boldsymbol 0           \\
    -\boldsymbol I & - \boldsymbol R_{0}      & \boldsymbol G           \\
    \boldsymbol 0              & -\boldsymbol G^\intercal & -\boldsymbol R_{\rm uu}
  \end{bmatrix}
  \underbrace{
    \begin{bmatrix}
      \nabla V           \\
      \boldsymbol M^{-1} \boldsymbol p_{\rm f} \\
      \boldsymbol u
    \end{bmatrix}}_{\boldsymbol e},
  \label{eq:dynamics_assembly}
\end{equation}
$$

with 

$$
\boldsymbol R_{\rm uu} = \begin{bmatrix}
        -R_{\rm g}         & R_{\rm g}  \\
        R_{\rm g}          & -R_{\rm g}         
      \end{bmatrix}, 
      \quad
\boldsymbol G = \begin{bmatrix}
    {\boldsymbol S}^{\intercal}_{\rm f-} \\\
    {\boldsymbol S}^{\intercal}_{\rm f+}
\end{bmatrix}.
$$

In the previous construction, the flow model can be viewed as a dissipative transformer between fluid flow/pressure and vocal folds velocity/force.


## Example models

### Body-cover (Story et Titze 1995[@story1995voice])

The body-cover model is a three mass model depicted in the following schematic.

<img src="../../../medias/larynx/body_cover_scheme.svg"; width=40%; style="display: block; margin: auto;"; alt="Vocal fold schematic: body cover model.">

The model can be fitted to the structure in \eqref{eq:dynamics_vf} quite easily, at least for the linear part of it. The state is composed of the displacements $\boldsymbol q_{\rm f} = [q_l, q_u, q_b]^\intercal$ around the rest positions and momenta $\boldsymbol p_{\rm f} = [p_l, p_u, p_b]^\intercal$ of the three masses. In the code, the convention is that higher displacement values corresponds to lower glottis opening (see remark above for implications on signs). The diagonal mass matrix writes 

$$
    \boldsymbol{M} = \begin{bmatrix}
        m_l & 0 & 0 \\\
        0 & m_u & 0 \\\
        0 & 0 & m_b 
    \end{bmatrix}.
$$

From the three d.o.f.s, the four spring elongations $\boldsymbol e = [e_l, e_u, e_b, e_{lu}]$ can be recovered as 

$$
\boldsymbol e = 
\underbrace{
\begin{bmatrix}
    1 & 0 & -1 \\\
    0 & 1 & -1 \\\
    0 & 0 & 1 \\\
    1 & -1 & 0 
\end{bmatrix}}_{\boldsymbol{W}}
\boldsymbol{q_{\rm f}}
$$

The stiffness matrix then writes 

$$
\boldsymbol{K} = 
\boldsymbol{W}^\intercal
\begin{bmatrix}
    k_l & 0 & 0 & 0 \\\
    0 & k_u & 0 & 0 \\\
    0 & 0 & k_b & 0\\\
    0 & 0 & 0 & k_{lu}
\end{bmatrix}
\boldsymbol{W},
$$

which is by construction semi-positive definite under the condition that all stiffness parameters are positive. Similarly, the dissipation matrix may be built as

$$
\boldsymbol{R_0} = 
\boldsymbol{W}^\intercal
\begin{bmatrix}
    r_l & 0 & 0 & 0 \\\
    0 & r_u & 0 & 0 \\\
    0 & 0 & r_b & 0\\\
    0 & 0 & 0 & r_{lu}
\end{bmatrix}
\boldsymbol{W}
$$
Note however than in the original model, and in this implementation, $r_{lu}$ is set to zero.

<details><summary>State dependent dissipation matrix  </summary>

The mechanical dissipation matrix is in fact state dependent, in a switching manner, as it should increase at contact. In the current state of the code, this is accounted for in a non power preserving manner at the discrete time level. In view of the relatively low values of the vocal fold model eigenvalue, this has no impact on stability. However, for bigger models, this may not be the case. The scheme can be slightly modified to allow for the dissipation matrix to be the sum of a fixed arbitrary part and of a state-dependent diagonal part, without reducing the computational efficiency. We used a similar trick for string simulations in Risse et. al 2025[@risse2025power].

</details>

With the above elements, the linear part of the vocal folds dynamics is fully described. In addition to that, $E_{\rm nl}(\boldsymbol q_{\rm f})$ encodes the mechanical nonlinearity. For a pair of vocal folds, this nonlinearity includes two elements:

- A cubic spring nonlinearity, which derives from a function of the elongations 
    $$
        E_{\rm material}(\boldsymbol e(\boldsymbol q_{\rm f})) = \frac{1}{4} [k_l, k_u, k_b, k_{lu}]^\intercal \left(\left(\frac{\boldsymbol e}{e_{\rm ref}}\right)^2  \odot \boldsymbol e^2 \right),
    $$
    where $e_{\rm ref}$ is a reference elongations such that the effective stiffness of the system is 4 times higher than around the rest position for $e = e_{\rm ref}$.
- A contact nonlinearity, consisting of an asymmetric power law, active only for interpenetrating vocal folds. Considering $\boldsymbol q_{c}(\boldsymbol{q}_{\rm f})$ to be a measure of the vocal fold interpenetration, the contact restoring force writes
    $$
        E_{\rm c}(\boldsymbol q_{c}(\boldsymbol{q}_{\rm f})) = \frac{1}{\alpha_c + 1} k_c  q_c^{\alpha_c} * (q_c > 0).
    $$
    Note that for a pair of non-symmetric vocal folds, $\boldsymbol q_{c}$ is a function of the concatenated state of both vocal folds.


### A symmetrical two mass model (Lous et. al 1998[@lous1998symmetrical])