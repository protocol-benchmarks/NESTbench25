###This script demonstrates the use of automatic differentiation to compute the optimal (minimum work) protocol for dragging an overdamped Brownian particle with a harmonic potential. It uses a parameterization for the protocol that comprises N discrete points with linear interpolation between them. ###

import time
import functools
import numpy as onp
import tqdm
import time

import jax
import jax.numpy as jnp
from jax.scipy import ndimage
from jax import jit, grad, vmap, random, lax, ops, tree_util
import optax

from jax_md import space, minimize, simulate, energy, quantity, smap

f32 = jnp.float32
f64 = jnp.float64

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import ArtistAnimation
from matplotlib import rc
from matplotlib import colors
rc('animation', html='jshtml')
import pandas as pd
import pickle
import csv
import make_animation as animate
import utils as utils
  
###Parameters used for optimization###
N = 1 #particle number
dim = 1 #dimension
box_size = 20
end_time = 1#2.69
dt = 1e-3 #this needs to be sufficiently small.
simulation_steps = int((end_time)/dt)+1 #{0, dt, 2dt, ..., end time} add 1 for t=0
opt_steps = 100
batch_size = 5000
temperature = 1. #kBT
key = random.PRNGKey(0)
key, split = random.split(key)
#mass is default = 1
#gamma is defaulted to 1
#thus diffusion constant D is D = 1./(beta*gamma*mass) #s^(-1) = kT = 1.0

init_position = 5. #particle initial position
r0_init = 5.
r0_final = 10.
trap_stiffness = 1.
Neq=100 # number of steps to equilibrate for. 

energy_fn = utils.Vspring(box_size, trap_stiffness)
displacement_fn, shift_fn = space.periodic(box_size)
init_schedule = jnp.linspace(5., 10., 8) # start with a linear protocol


summaries = []
schedules = []
all_works = []
optimizer = optax.adam(.1)
opt_state = optimizer.init(init_schedule)
grad_fxn = utils.estimate_gradient(batch_size, energy_fn, init_position, r0_init, r0_final, Neq, shift_fn, simulation_steps, dt, temperature)
new_params = init_schedule
#save_every = 10 
###Optimization loop###
for j in tqdm.trange(opt_steps,position=0):
  params = new_params
  key, split = random.split(key)
  grad, (works, summary) = grad_fxn(params, split)
  updates, opt_state = optimizer.update(grad, opt_state)
  new_params = optax.apply_updates(params, updates)
  #if j % (opt_steps // save_every) == 0:
  all_works.append(works)
  summaries.append(summary)
  schedules.append((j,) + (new_params,))


all_works = tree_util.tree_map(lambda *args: jnp.stack(args), *all_works)
print("final average dissipated work: ", onp.mean(all_works[-1]))

savedir = 'your_directory_here_'
afile = open(savedir+'schedules.pkl', 'wb')
pickle.dump(schedules, afile)
afile.close()
afile = open(savedir+'works.pkl', 'wb')
pickle.dump(all_works, afile)
afile.close()

###Plot results###
labelfont=26
tickfont=26
_, (ax0, ax1) = plt.subplots(1, 2, figsize=[16, 8])
utils.plot_with_stddev(all_works.T, ax=ax0, alpha=0.1)
ax1.set_title('Loss')
ax0.set_xlabel("Optimization iteration",fontsize=labelfont)
ax0.set_ylabel("Dissipated work (kT)",fontsize=labelfont)
ax0.tick_params(axis='x', labelsize=tickfont)
ax0.tick_params(axis='y', labelsize=tickfont)
colors = ['']
ax0.set_ylim(0.,14.)

trap_fn = utils.make_trap_fxn(init_schedule,simulation_steps,r0_init,r0_final)
init_sched = trap_fn(jnp.arange(simulation_steps))
ax1.plot(jnp.arange(simulation_steps), init_sched, lw=5, label='Initial Guess')
plot_every = 10
for j, sched in schedules:
  if((j+1)%plot_every == 0):
    trap_fn = utils.make_trap_fxn(sched,simulation_steps,r0_init,r0_final)
    full_sched = trap_fn(jnp.arange(simulation_steps))
    ax1.plot(jnp.arange(simulation_steps), full_sched,lw=5, label=f'Step {j+1}')

#plot final:
j, sched = schedules[-1]
trap_fn = utils.make_trap_fxn(sched,simulation_steps,r0_init,r0_final)
full_sched = trap_fn(jnp.arange(simulation_steps))

#exact solution from Schmiedl & Seifert 2007:
lambda_t = utils.make_trap_fxn_theoretical(dt, simulation_steps, end_time, r0_init, r0_final)
trap=[]
for k in (jnp.arange(simulation_steps)):
  trap.append(lambda_t(k))
ax1.plot(jnp.arange(simulation_steps), trap, '--',lw=5, label='Theory')



ax1.legend(fontsize=labelfont)#
ax1.set_title('Schedule evolution')
ax1.tick_params(axis='x', labelsize=tickfont)
ax1.tick_params(axis='y', labelsize=tickfont)
ax1.set_xlabel("Time (s)",fontsize=labelfont)
ax1.set_ylabel("Trap displacement (nm)",fontsize=labelfont)
plt.show()

###optional: create animation of particle distribution evolving under final protocol###
anim = animate.create_histogram_animation(
    summaries=summaries,
    schedules=schedules,
    energy_fn=energy_fn,
    box_size=box_size,
    simulation_steps=simulation_steps,
    batch_size=batch_size,
    dt=dt,
    r0_init=r0_init,
    r0_final=r0_final,
    trap_stiffness=trap_stiffness,
    savedir=savedir,
    fps=30,
    skip_steps=10  # Adjust this to control animation speed/size
)
