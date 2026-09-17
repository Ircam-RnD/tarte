---
title: Z transform of a state space model in discrete time
---


Consider a LTI state space model of the form 

$$
\begin{align*}
    x^{n+1} &= A x^n + B u^n, \\
    y^{n+1} &= C x^n + D u^n.
\end{align*}
$$

Applying the Z-transform yields

$$
\begin{align*}
    z X(z) -z x_0&= A X(z) + B U(Z), \\
    Y(z) &= C X(z) + D U(z),
\end{align*}
$$

or equivalently

$$
\begin{align*}
    (z -A) X(z) &= z x_0 + B U(Z), \\
    Y(z) &= C (z -A)^{-1} z x_0 +  \left(C (z -A)^{-1}B + D\right) U(z).
\end{align*}
$$

In order to compute a transfer function, consider $x_0 = 0$. We then get 

$$
    \begin{align*}
    Y(z) =  \left(C (z -A)^{-1}B + D\right) U(z).
\end{align*}
$$

If $y$ is taken to be a power-balanced output to $u$ , the transfer function $H(z) = y/u$ should correspond to an impedance-like (or admittance) quantity. For a power-preserving discretization, this quantity is expected to be passive (such that the poles of  $\left(C (z -A)^{-1}B + D\right)$ must have positive real part).

