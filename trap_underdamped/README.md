# Underdamped Particle in a Harmonic Trap Benchmark

This code implements a benchmark system for testing methods of protocol optimization in nonequilibrium statistical mechanics. The system consists of an underdamped Langevin particle in a moving harmonic trap, with the objective of minimizing work during trap translation.

## Problem Description

We consider an underdamped particle with position $x$ experiencing the Langevin dynamics:

$$m \ddot{x} + \gamma \dot{x} = -\frac{\partial U(x,\lambda(t))}{\partial x} + \left(2 \gamma k_B T\right)^{1/2}\xi(t)$$

where $m$ is the particle mass, $\gamma$ is the damping coefficient, and $\xi(t)$ is Gaussian white noise with zero mean and correlation function $\langle\xi(t)\xi(t')\rangle=\delta(t-t')$.

The particle is confined in a harmonic potential:

$$U(x,\lambda) = \frac{k}{2}(x-\lambda(t))^2$$

where $k$ is the spring constant and $\lambda(t)$ is the time-dependent position of the trap center.

The natural scales of this problem are:
- Length: $\sigma_0 = \sqrt{k_B T/k} = 1$ nm
- Time: $\omega_0^{-1} = \sqrt{m/k}$ where $\omega_0 = 2\pi \times 1090$ Hz
- Quality factor: $Q = m\omega_0/\gamma = 7$

These values are motivated by micromechanical cantilever experiments, so the "particle" in this problem represents the position of a cantilever tip.

## Optimization Objective

The control parameter of this problem is the trap center position $\lambda(t)$. The aim is to move the trap from position $\lambda(0)=0$ to position $\lambda(t_f)=5$ in finite time $t_f = 0.15 \times (2\pi/\omega_0)$, which is short relative to the fundamental oscillation time.

The order parameter to be minimized is the mean work required to perform many independent realizations of this translation.

## Optimal Protocol

The optimal protocol for this problem is known analytically and is given by:

$$\lambda^*(t) = \lambda(t_f)\frac{Q\omega_0 t + 1}{Q\omega_0 t_f + 2} + \Lambda[\delta(t) - \delta(t-t_f)]$$

where $\Lambda = \lambda(t_f)/(\omega_0^2 t_f + 2\omega_0 Q^{-1})$.

This protocol contains delta-function kicks at the beginning and end, which effect an instantaneous change in the velocity of the particle. This poses a challenge for optimization algorithms since direct control of the particle velocity is not allowed, only control of the trap position $\lambda(t)$.

The mean work done under this optimal protocol is:

$$W^* = \frac{k\lambda(t_f)^2}{2 + Q\omega_0 t_f} \approx 2.90787$$

in units of $k_B T$.

