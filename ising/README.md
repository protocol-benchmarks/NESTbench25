# Ising Model Simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A high-performance implementation of the Ising model for studying phase transitions in ferromagnetic materials. This project provides both a standalone executable and a library that can be integrated into other applications.

## Table of Contents

- [Ising Model Simulation](#ising-model-simulation)
  - [Table of Contents](#table-of-contents)
  - [Overview](#overview)
  - [Features](#features)
  - [Directory Structure](#directory-structure)
  - [Requirements](#requirements)
    - [Core Requirements](#core-requirements)
    - [Visualization Requirements (Python)](#visualization-requirements-python)
  - [Installation](#installation)
    - [Quick Installation](#quick-installation)
    - [Manual Installation](#manual-installation)
  - [Usage](#usage)
    - [Standalone Executable](#standalone-executable)
    - [Library Integration](#library-integration)
    - [Visualization](#visualization)
    - [Examples](#examples)
    - [Testing](#testing)
  - [Protocol Files](#protocol-files)
  - [Output Files](#output-files)
  - [Contributing](#contributing)
  - [License](#license)

## Overview

The Ising model is a mathematical model of ferromagnetism in statistical mechanics. This implementation simulates a 2D lattice of spins that can be in one of two states (+1 or -1) and interact with their nearest neighbors. The simulation can be used to study phase transitions and critical phenomena.

## Features

- **High-Performance C/C++ Implementation**: Efficient simulation engine optimized for speed
- **Configurable Protocols**: Customize simulation parameters through protocol files
- **Visualization Tools**: Generate snapshots and animations of the simulation
- **Library Integration**: Use as a standalone executable or integrate as a library
- **Comprehensive Examples**: Learn how to use the library through examples
- **Test Suite**: Verify functionality through automated tests

## Directory Structure

```
ising/
├── src/                      # Source code directory
│   ├── engine/               # Core simulation engine
│   │   ├── engine_ising.c    # Main simulation code
│   │   ├── engine_ising.h    # Header file with API definitions
│   │   └── main.cpp          # Standalone executable entry point
│   └── visualization/        # Visualization tools
│       ├── movie.py          # Script to generate images and animations
│       └── movie_example.py  # Example visualization with synthetic data
├── data/                     # Data files directory
│   ├── protocols/            # Simulation protocols
│   │   ├── default_protocol.dat            # Default parameters
│   │   └── input_control_parameters_learned.dat  # Optimized parameters
│   └── output/               # Directory for simulation outputs
├── build/                    # Build artifacts (generated)
├── docs/                     # Documentation
│   └── README.md             # Detailed project documentation
├── examples/                 # Example usage of the library/executable
│   ├── simple_simulation.c   # Simple example of library usage
│   └── Makefile              # Build system for examples
├── tests/                    # Test cases
│   ├── test_basic.c          # Basic functionality tests
│   └── Makefile              # Build system for tests
├── venv/                     # Python virtual environment (generated)
└── Makefile                  # Main build system
```

## Requirements

### Core Requirements
- C++ compiler with C++11 support (g++ or clang++)
- GNU Make

### Visualization Requirements (Python)
- Python 3.6+
- matplotlib
- numpy
- imageio with ffmpeg support

## Installation

### Quick Installation

For a quick setup, use the provided installation script:

```bash
git clone https://github.com/yourusername/ising-model.git
cd ising-model
./setup.sh
```

This script will:
1. Create a Python virtual environment
2. Install required Python packages
3. Create necessary directories
4. Build the standalone executable

### Manual Installation

If you prefer to install manually:

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/ising-model.git
   cd ising-model
   ```

2. Set up Python environment for visualization:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```

3. Build the project:
   ```bash
   make standalone  # For standalone executable
   # OR
   make library     # For library integration
   ```

## Usage

### Standalone Executable

First, ensure the project is built. If you haven't already (e.g., via `./setup.sh`), build the standalone executable from the `ising` directory:

```bash
make standalone
```

To run the simulation:

1.  **Execute the simulation program:**
    ```bash
    ./build/sim
    ```
    This will:
    - Load the default protocol from `data/protocols/default_protocol.dat`.
    - Run the simulation.
    - Calculate the order parameter.
    - Attempt to generate visualization files (`output_picture.pdf`, `output_movie.mp4`) in `data/output/` by calling a Python script.

2.  **Manual Visualization (if step 1 did not produce visualization files):**
    If the visualization files were not created (e.g., due to the Python environment not being accessible to the compiled program), you can generate them manually. From the `ising` directory:
    ```bash
    source venv/bin/activate  # Activate the Python virtual environment
    python src/visualization/movie.py
    deactivate              # Optional: deactivate the virtual environment
    ```
    This will create `output_picture.pdf` and `output_movie.mp4` in the `data/output/` directory.

### Library Integration

Build the library:

```bash
make library
```

This creates `build/libengine_ising.a` which you can link with your own projects.

Include the header in your code:

```c
#include "path/to/engine_ising.h"
```

Link with the library:

```bash
g++ -o your_program your_program.c -L/path/to/build -lengine_ising -lm
```

### Visualization

The `./build/sim` executable attempts to automatically generate visualization outputs (`output_picture.pdf` and `output_movie.mp4`) in `data/output/` after a simulation completes.

If this automatic step fails (e.g., Python dependencies like `imageio` are not found by the C++ program's environment), you can generate the visualizations manually from the `ising` directory:

```bash
source venv/bin/activate  # Ensure the virtual environment is active
python src/visualization/movie.py
deactivate              # Optional: deactivate the virtual environment
```

This script uses the data generated by `./build/sim` to create:
- `output_picture.pdf`: A grid of snapshots from the simulation.
- `output_movie.mp4`: An animation of the simulation.

### Examples

The `examples/` directory contains sample code showing how to use the library:

```bash
cd examples
make
./simple_simulation
```

### Testing

Run the test suite to verify functionality:

```bash
cd tests
make run
```

## Protocol Files

The simulation behavior is controlled by protocol files:

- `data/protocols/default_protocol.dat`: Default simulation parameters
- `data/protocols/input_control_parameters_learned.dat`: Optimized parameters

Each protocol file contains rows with two columns:
1. Time parameter (0.0 to 1.0)
2. Control parameter value

## Output Files

The simulation generates several output files:

- `report_control_parameters_*.dat`: Control parameter values over time
- `report_lattice_time_*.dat`: Lattice state snapshots at different times
- `report_entprod_histogram.dat`: Entropy production histogram
- `output_picture.pdf`: Visualization of lattice snapshots
- `output_movie.mp4`: Animation of the simulation

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.