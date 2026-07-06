---
title: Simple radiation model
---


At the mouth aperture, the boundary condition is classically represented by a radiation impedance relating pressure and volume flow for time-invariant geometries
$$
  Z_{\rm rad}(\omega) = \frac{P_{\rm lips}(\omega) - P_0}{Q_{\rm lips}(\omega)} = \frac{P_{\rm rad}(\omega)}{Q_{\rm rad}(\omega)}.
$$
Several models of the mouth radiation impedance exist in the literature (see e.g. Wakita 1987[@wakita1978toward], Chalker 1985 [@chalker1985models], Vojnovic 2005[@vojnovic2005improved]). From an analytical viewpoint, mathematical expressions for the radiation impedance of a piston in a spherical baffle are available (see e.g. Morse 1986[@morse1986theoretical]) and constitute a good approximation. However, these expressions involve infinite sums of Legendre and Bessel functions, making them hard to evaluate. Instead, the expression for a piston with infinite plane baffle can be used, as it yields acceptable results when the radius of the aperture (mouth) is small against the radius of the sphere (head), and is easier to compute. It writes
$$
\begin{equation}
  Z_{\rm rad} (\omega) = \frac{\rho_0 c_0}{A_{\rm lips}}\left(
    1
    - \frac{J_1 (2 k r_{\rm lips}) - j S_1 (2 k r_{\rm lips})}{kr_{\rm lips}}
  \right), \text{ with } k = \frac{\omega}{c_0},
  \label{eq:radiation_imp_infinite}
\end{equation}
$$
and where $A_{\rm lips} = \pi r_{\rm lips}^2$ are the mouth aperture area and radius, respectively. Function $J_1(x)$ is the first-order Bessel function (see Abramowitz 1965[@abramowitz1965handbook] (chapter 4)) and $S_1(x)$ the first-order Struve function (see Abramowitz 1965[@abramowitz1965handbook] (chapter 12)). This representation cannot easily be converted into a time-domain representation.
As a result, a low frequency approximation represented by a filter with a (small) finite number of coefficients is often used and characterized by an electrical circuit equivalent with coefficients depending on the mouth aperture. 

The simplest approximation consisting of a first-order high-pass filter is used here, and corresponds to a resistance $R_{\rm rad}$ in parallel with an acoustic mass $L_{\rm rad}$. The expression of the radiation impedance then reduces to
$$
\begin{equation}
  Z_{\rm rad}(\omega) = R_{\rm rad} \frac{j \frac{\omega}{\omega_0}}{1 + j\frac{\omega}{\omega_0}}, \text{ with } R_{\rm rad} \geq \text{ and } \omega_0 = \frac{R_{\rm rad}}{L_{\rm rad}}.
  \label{eq:imp_first_order}
\end{equation}
$$
<!-- \radiation -->
Following Flanagan 2013[@flanagan2013speech] (section 4.2, equation 4.4), coefficients are obtained from a Taylor series expansion of \eqref{eq:radiation_imp_infinite} as
$$
\begin{equation}
  R_{\rm rad} = Z_0 \frac{128}{9\pi^2}, \; L_{\rm rad} = Z_0 \frac{8r_{\rm lips}}{3\pi c_0}, \text{ with } Z_0= \frac{\rho_0 c_0}{\pi r_{\rm lips}^2},
  \label{eq:rad_coeff}
\end{equation}
$$
where $r_{\rm lips}$ is the radius of the mouth aperture. 

<img src="../../../medias/radiation/radiation.svg"; width=50%; style="display: block; margin: auto;"; alt="simulated radiation impedance curves">
<!-- \begin{figure}[hbtp]
  \centering
  \includegraphics{medias/voice/radiation.pdf}
  \caption{Normalized (w.r.t. $Z_0$) radiation impedances as computed from \eqref{eq:radiation_imp_infinite} (denoted ``Reference'') and from the first order approximation \eqref{eq:imp_first_order} (denoted ``First order''). Both the resistance (real part) and reactance (imaginary part) are shown, for realistic mouth openings corresponding to vowels `o' and `a' (with values from \cite{story1996vocal}). }
  \label{fig:radiation_imp}
\end{figure} -->

The above figure presents the impedance curves obtained from the analytical model and its first order approximation in a vocal tract like configuration. The validity range of the simplified model depends on the mouth aperture and is reduced for wide openings. For vowel 'o' with $A_{\rm lips} = 0.14$ cm$^2$, both models are nearly equivalent in the studied frequency range. However, significant deviation is visible for vowel 'a' with $A_{\rm lips} = 5.03$ cm$^2$ and frequencies above $4000$ Hz. Plane wave propagation is considered in the vocal tract to allow for the 3D to 1D reduction of the fluid model. As this approximation is expected to be valid up to $\sim 5000$ Hz (see e.g. Blandin 2015 [@blandin2015effects]), the simplified radiation model is suitable in the validity range of the approach for most of the mouth aperture values. 

Writing a port-Hamiltonian formulation of the simplified radiation impedance is straightforward from \eqref{eq:imp_first_order}. The state $x = [Q_{\rm reac}]$ corresponds to the fluid volume flow associated with the acoustic mass $L_{\rm rad}$. The corresponding Hamiltonian writes $H (Q_{\rm reac}) = \frac{1}{2} L_{\rm rad} Q_{\rm reac}^2$, with associated dynamics

$$
\begin{equation}
    \begin{bmatrix}
      \dot{Q}_{\rm reac}\\\
      y_r = - Q_{\rm rad}
    \end{bmatrix}
    =
    \begin{bmatrix}
       0 & \frac{1}{L_{\rm rad}} \\\
       -\frac{1}{L_{\rm rad}} & -\frac{1}{R_{\rm rad}}
    \end{bmatrix}
    \begin{bmatrix}
      L_{\rm rad} Q_{\rm reac} \\\
      u_r = P_{\rm rad}
    \end{bmatrix}.
    \label{eq:causality1}
\end{equation}
$$

Alternatively, a second formulation with reversed causality is obtained as

$$
\begin{equation}
    \begin{bmatrix}
      \dot{Q}_{\rm reac}\\\
      y_r = - P_{\rm rad}
    \end{bmatrix}
    =
    \begin{bmatrix}
       -\frac{R_{\rm rad}}{L_{\rm rad}^2} & \frac{R_{\rm rad}}{L_{\rm rad}} \\\
       \frac{R_{\rm rad}}{L_{\rm rad}} & - R_{\rm rad}
    \end{bmatrix}
    \begin{bmatrix}
      L_{\rm rad} Q_{\rm reac} \\\
      u_r = Q_{\rm rad}
    \end{bmatrix}.
    \label{eq:causality2}
\end{equation}
$$


<!-- \remark{Remark: Causality of the radiation component}{
  Some spatial discretization methods do not allow for the direct construction of mixed boundary conditions. This is the case for e.g. the partitioned finite element method \cite{cardoso2021partitioned}. While solutions are designed to tackle this issue (see e.g. \cite{brugnoli2020partitioned}), they either involve the introduction of Lagrange multipliers in the system, or an arbitrary subdivision of the domain which complexifies the assembly. 

  In prevision of the connection of the radiation component with the rest of the system, it is useful to have both representation \eqref{eq:causality1} and \eqref{eq:causality2} to suit match causality of the discrete fluid model.
} -->