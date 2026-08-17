## Code Structure

- `engine_trap_overdamped.c`: Main simulation engine that implements the overdamped Langevin dynamics
- `engine_trap_overdamped.h`: Header file for using the engine as a library
- `ga_process.c`: Implements a single member of a genetic population for protocol optimization
- `ga.c`: Runs the genetic algorithm for protocol learning
- `Makefile`: Compilation instructions
- `movie.py`: Python script for visualization
- `input_control_parameters_learned.dat`: A pre-learned protocol

## Usage

The user can specify the time-dependent protocol $\lambda(t)$ via the external file `input_control_parameters.dat`. The dimensions of this file are specified by the function `load_protocol()`.

To run the engine as a standalone code:
1. Uncomment the `main()` function in `engine_trap_overdamped.c`
2. Compile using `make standalone`
3. Run the executable

The code includes several useful functions:
- `load_default_protocol()`: Loads the optimal protocol
- `visualize_protocol()`: Calculates mean work using 10^5 trajectories and outputs visualization
- `final_answer()`: Calculates the order parameter (mean work) over 100 samples of 10^4 trajectories

For protocol learning, the code can be used as an external library as demonstrated in the neuroevolutionary approach included in this repository.
