# Ising Model Benchmark for Protocol Design

This code implements a benchmark system for testing methods of protocol optimization in nonequilibrium statistical mechanics. The system consists of a 2D Ising model undergoing state change, with the objective of minimizing entropy production.

## Problem Description

We consider the 2D Ising model on a square lattice of $N=32^2$ sites with periodic boundary conditions. On each site $i$ is a binary spin $S_i = \pm 1$. The lattice energy function is:

$$E = -J\sum_{\langle ij \rangle} S_i S_j -h\sum_{i=1}^N S_i$$

where $J=1$ is the Ising coupling, $h$ is the magnetic field, and the first sum runs over all nearest-neighbor bonds.

The model evolves by Glauber Monte Carlo dynamics. At each time step, a lattice site is chosen and a spin flip is proposed, which is accepted with probability:

$$P_{\text{Glauber}}(\Delta E)=\left( 1+ \exp(\Delta E/T) \right)^{-1}$$

where $\Delta E$ is the energy change under the proposed move.

## Optimization Objective

The control parameters of this problem are temperature and magnetic field, $\mathbf{c}(t) = (T(t), h(t))$. The aim is to change these parameters from the values $\mathbf{c}(0) = (0.7,-1)$ to the values $\mathbf{c}(t_f) = (0.7,1)$, in $t_f=100 N$ Monte Carlo steps (100 Monte Carlo sweeps). 

The order parameter to be minimized is the mean entropy production $\langle\sigma\rangle$, where the entropy produced over the course of a simulation is:

$$\sigma=E_f/T_f - E_0/T_0 - \sum_{k=1}^{t_f} \Delta E_k/T_k$$

Here $E_f$ and $E_0$ are the final and initial energies of the system, and $\Delta E_k$ and $T_k$ are the energy change and temperature at step $k$ of the simulation.

## Default Protocol

The default protocol is based on the work of Rotskoff and Crooks (2015), which is the optimal protocol within the near-equilibrium approximation. This protocol avoids the first-order phase transition line and the critical point of the Ising model.

