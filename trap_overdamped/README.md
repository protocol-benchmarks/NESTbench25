# Overdamped Particle in Harmonic Trap Benchmark for Protocol Design

This code implements a benchmark system for testing methods of protocol optimization in nonequilibrium statistical mechanics. The system consists of an overdamped Langevin particle in a moving harmonic trap, with the objective of minimizing work during trap translation.

## Problem Description

We consider an overdamped Langevin particle at position $x$ in a potential:

$$U(x,\lambda) = \frac{1}{2}(x-\lambda(t))^2$$

in units such that $k_B T = 1$. The particle undergoes the Langevin dynamics:

$$\gamma \dot{x} = -\frac{\partial U(x,\lambda(t))}{\partial x} + \left(2\gamma k_B T\right)^{1/2}\xi(t)$$

where $\gamma$ is the damping coefficient (set to 1), and $\xi(t)$ is Gaussian white noise with zero mean and correlation function $\langle\xi(t)\xi(t')\rangle = \delta(t-t')$.

The dynamics are integrated using a first-order Euler scheme:

$$x(t+\Delta t) = x(t) - \partial_x U(x,\lambda(t))\Delta t + \sqrt{2\Delta t}X$$

where $\Delta t = 10^{-3}$ and $X$ is a Gaussian random number with zero mean and unit variance.

## Optimization Objective

The control parameter of this problem is the trap center $\mathbf{c}(t) = \lambda(t)$. The aim is to move the trap center from an initial position $\lambda(0) = 0$ to a final position $\lambda(t_f) = 5$, in finite time $t_f = 1$, minimizing the work:

$$\langle W \rangle = \int_0^{t_f} dt\, \dot{c}_0(t) \left\langle \frac{\partial U(x, \lambda)}{\partial \lambda} \right\rangle$$

averaged over many realizations of the process.

## Optimal Protocol

The optimal (mean-work-minimizing) protocol is known analytically (Schmiedl and Seifert, 2007). It has a linear form:

$$\lambda^\star(t) = \lambda(t_f) \frac{t+1}{t_f+2}$$

for $0 < t < t_f$, with jump discontinuities at the start $(t=0)$ and end $(t=t_f)$. This protocol produces mean work:

$$\langle W \rangle^\star = \frac{\lambda(t_f)^2}{t_f+2} = \frac{25}{3} \approx 8.333$$

The jump discontinuities are a generic feature of many optimal protocols. Learning them requires a protocol ansatz that is not constrained to smoothly approach its boundary values.
