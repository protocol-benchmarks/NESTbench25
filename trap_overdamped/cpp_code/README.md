## Code Structure

- `engine_trap_overdamped.cpp`: Main simulation engine that implements the overdamped Langevin dynamics
- `engine_trap_overdamped.h`: Header file for using the engine as a library
- `ga_process.cpp`: Implements a single member of a genetic population for protocol optimization
- `ga.cpp`: Runs the genetic algorithm for protocol learning
- `Makefile`: Compilation instructions
- `movie.py`: Python script for visualization
- `input_control_parameters_learned.dat`: A pre-learned protocol

## Usage

The user can specify the time-dependent protocol $\lambda(t)$ via the external file `input_control_parameters.dat`. The dimensions of this file are specified by the function `load_protocol()`.

To run the engine as a standalone code, compile with `make standalone` and run
`./sim` (no source modification is needed; `make library` automatically excludes
`main()`). Run `./sim -h` for command-line options: the protocol (`-p <file>`),
trajectory time (`-t`), boundary values (`-b`), and number of trajectories (`-n`)
can all be set without recompiling. Results are written to `report_answer.dat`.
To build the genetic-algorithm driver and worker, run `make ga`; SLURM settings
for the generated job scripts are read from the environment variables
`NESTBENCH_SLURM_PARTITION`, `NESTBENCH_SLURM_ACCOUNT`, `NESTBENCH_SLURM_QOS`,
and `NESTBENCH_SLURM_TIME`.

The code includes several useful functions:
- `load_default_protocol()`: Loads the optimal protocol
- `visualize_protocol()`: Calculates mean work using 10^5 trajectories and outputs visualization (run `./sim -v`; requires python with matplotlib/imageio and ffmpeg)
- `final_answer()`: Calculates the order parameter (mean work) over 100 samples of 10^4 trajectories

For protocol learning, the code can be used as an external library as demonstrated in the neuroevolutionary approach included in this repository.
