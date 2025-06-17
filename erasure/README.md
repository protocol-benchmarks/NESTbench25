# Logic Erasure Benchmark for Protocol Design

This code implements a benchmark system for testing methods of protocol optimization in nonequilibrium statistical mechanics. The system simulates logic erasure using an underdamped Langevin particle in a double-well potential, with the objective of maximizing the probability of successful erasure.

## Problem Description

We consider an underdamped Langevin particle that can start in either well of a double-well potential, representing a one-bit memory (left well = 0, right well = 1). The particle's motion is governed by:

$$m \ddot{x} + \gamma \dot{x} = -\frac{\partial U(x,\mathbf{c}(t))}{\partial x} + \left( 2 \gamma k_B T\right)^{1/2}\xi(t)$$

where $\xi(t)$ is Gaussian white noise with zero mean and correlation function $\langle\xi(t) \xi(t')\rangle = \delta(t-t')$. The potential is:

$$U(x, \mathbf{c}(t)) = \frac{k}{2} \left( x - S(x - c_0(t))\, c_1(t) \right)^2 + k c_0(t)\, c_1(t) \left( S(x - c_0(t)) + S(c_0(t)) \right)$$

where $S$ is the sign function. The system parameters are based on a micromechanical cantilever with:
- Fundamental frequency: $f_0 = 1090$ Hz
- Quality factor: $Q = 7$

## Optimization Objective

The control parameters of this problem are the potential shape parameters $\mathbf{c}(t) = (c_0(t), c_1(t))$. The aim is to manipulate these parameters to bring the particle to the left-hand well (logic state 0) at time $t_f$, regardless of its initial state, thus erasing the memory.

The protocol must operate within a very short time $t_f = 0.5 \times (2\pi/\omega_0)$, which is half the fundamental oscillation period of the system. This short duration makes perfect erasure impossible, and the objective is to minimize the erasure error rate.

The boundary conditions are $\mathbf{c}(0) = \mathbf{c}(t_f) = (0, 5)$, which correspond to a symmetric double-well potential with an energy barrier of $12.5$ (in units of $k_B T$).

## Default Protocol

The default protocol attempts to merge the double wells, translate the resulting single well to the left, and then reconstitute the double-well potential. While this approach works well for longer erasure times, it is ineffective for the short duration specified in this benchmark, achieving only about 53% successful erasure.

