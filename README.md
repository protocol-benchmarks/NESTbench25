# NESTbench25

This is the official repository for NESTbench25. For details about the problems, please see [our paper](https://arxiv.org/abs/2506.15122)

The benchmark contains 5 problems. Each directory contains a C++ and a Python implementation
(and, for `trap_overdamped`, an additional JAX implementation).

| Problem | Objective |
|---|---|
| `trap_overdamped` | minimize mean work, overdamped particle in a moving trap |
| `trap_underdamped` | minimize mean work, underdamped particle in a moving trap |
| `erasure` | minimize the error rate of fast logic erasure |
| `ising` | minimize entropy production during magnetization reversal |
| `abp` | state-to-state transformation of an active Brownian particle |

## Quick start

C++ (each `<problem>/cpp_code` directory):

```bash
make standalone    # builds ./sim
./sim              # runs the default protocol; results are written to report_answer.dat
./sim -h           # full list of options
```

All engines accept the same core options — no recompilation or source
modification is needed to change them:

```bash
./sim -p <file>          # protocol from a file (e.g. input_control_parameters_learned.dat)
./sim -n <n_traj>        # number of trajectories
./sim -t <tf>            # trajectory time
./sim -b <c_i,c_f,...>   # protocol boundary values (e.g. -b 0.7,0.7,-1,1 for ising T_i,T_f,h_i,h_f)
./sim -v                 # generate movies/histograms instead of the final answer
make library             # builds a static library for use with your own optimizer (main() excluded)
```

Python (each `<problem>/python_code` directory):

```bash
pip install -r requirements.txt   # torch, numpy, matplotlib
python engine_<problem>.py        # default protocol; results printed and written to report_answer.dat
python engine_<problem>.py -h     # same options as the C++ code
```

The Python engines run on CPU or CUDA automatically.

## Dependencies

- C++: any C++11 compiler (g++/clang++); no external libraries.
- Python: see each `requirements.txt` (torch, numpy, matplotlib).
- Visualization only: `ffmpeg` (movies) and optionally LaTeX (plot labels; the
  code falls back to matplotlib's mathtext when LaTeX is not installed).
- The SLURM scripts in `trap_overdamped/cpp_code` read cluster-specific
  settings (partition/account/qos) from `NESTBENCH_SLURM_*` environment
  variables or `sbatch` command-line flags; nothing is hard-coded.

## Docker

Two containerized environments are provided for fully reproducible runs:

```bash
# Full image: C++ engines, Python (CPU PyTorch) engines, ffmpeg, and the test suite
docker build -t nestbench25 .
docker run --rm nestbench25                     # runs the quick test suite
docker run --rm -w /nestbench/erasure/cpp_code nestbench25 \
    sh -c "./sim -p input_control_parameters_learned.dat && cat report_answer.dat"
docker run --rm -w /nestbench/abp/python_code nestbench25 \
    python engine_abp.py -n 10000

# Minimal C++-only image (no Python; much smaller and faster to build)
docker build -f Dockerfile.cpp -t nestbench25-cpp .
docker run --rm -w /nestbench/ising/cpp_code nestbench25-cpp \
    sh -c "./sim -n 100 -s 10 && cat report_answer.dat"
```

All engines are pre-built inside the images. Results are written to
`report_answer.dat` in the engine's directory; mount a volume (`-v`) to keep
outputs on the host.

## Tests

An automated test suite (pytest) checks that every C++ engine builds and runs,
that every Python engine runs, and that the two implementations agree:

```bash
pip install pytest
pytest tests/ -v
```

The quick tests run in a few minutes on a laptop. Set `NESTBENCH_FULL_TESTS=1`
to also run the slower cross-language agreement checks at higher statistics.
