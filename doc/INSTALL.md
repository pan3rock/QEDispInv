# Installation Guide for QEDispInv

This document provides step-by-step instructions to install the QEDispInv program,
along with dependency installation for various systems

## Prerequisites

### 1. C++17 Compatible Compiler

Verify your compiler supports C++17:

```bash
# For g++
g++ --version  # Requires version 7.0+
echo '#include <string_view>' | g++ -x c++ -std=c++17 -c -o /dev/null -

# For clang++
clang++ --version  # Requires version 5.0+
echo '#include <string_view>' | clang++ -x c++ -std=c++17 -c -o /dev/null -
```

Install or update compilers if needed:

- **Ubuntu/Debian**: `sudo apt-get install gcc g++ clang`
- **Fedora/RHEL**: `sudo dnf install gcc gcc-c++ clang`
- **Arch Linux**: `sudo pacman -S gcc clang`
- **macOS (Homebrew)**: `brew install gcc llvm`

Typically, installing GCC will also install gfortran at the same time.

### 2. Required Libraries

#### OpenMPI

- **Ubuntu/Debian**: `sudo apt-get install libopenmpi-dev`
- **Fedora/RHEL**: `sudo dnf install openmpi`
- **Arch Linux**: `sudo pacman -S openmpi`
- **macOS (Homebrew)**: `brew install open-mpi`

#### Eigen3

- **Ubuntu/Debian**: `sudo apt-get install libeigen3-dev`
- **Fedora/RHEL**: `sudo dnf install eigen3-devel`
- **Arch Linux**: `sudo pacman -S eigen`
- **macOS (Homebrew)**: `brew install eigen`

#### CMake

Install CMake (required for building the project):

- **Ubuntu/Debian**: `sudo apt-get install cmake`
- **Fedora/RHEL**: `sudo dnf install cmake`
- **Arch Linux**: `sudo pacman -S cmake`
- **macOS (Homebrew)**: `brew install cmake`

## Build Steps with CMake

1. **Clone the Repository**

   ```bash
   git clone https://github.com/pan3rock/QEDispInv.git
   cd QEDispInv
   ```

2. **Initialize Submodules**

   ```bash
   git submodule update --init
   ```

3. **Create and Configure Build Directory**

   ```bash
   # Create build directory
   mkdir -p build && cd build

   # Configure with CMake (basic version)
   cmake ..
   ```

4. **Build the Project**
   ```bash
   # Use -j flag to enable parallel compilation (adjust number based on CPU cores)
   make -j4
   ```

## Verification

After successful build, verify the command:

```bash
../bin/forward -h
../bin/inversion -h
../bin/secfunc -h
```

If the program displays a help message, the installation was successful.

### MPI Verification

Verify your MPI installation is working correctly:

```bash
# Check MPI version
mpirun --version
mpiexec --version

# Test MPI execution with 2 processes
mpirun -n 2 ../bin/inversion -h
```

## MPI Runtime Usage

The `inversion` executable uses MPI for parallel execution across multiple processes. Only the `inversion` tool requires MPI; `forward` and `secfunc` run serially.

### Important: MPI Environment Compatibility

**⚠️ Common Issue: Multiple MPI Installations**

If you have **conda/mamba** installed, your `mpirun` may default to the conda environment's MPI instead of the system MPI used during compilation. This causes runtime errors or the program appearing to hang without output.

**Symptoms**:
- Multiple copies of output appearing (e.g., "lmin=..." printed twice)
- Program hangs without producing results
- Error: "MPI processes failed to start"
- Segmentation fault during execution

**Solution**: Always use the system MPI that matches your compiled libraries:

```bash
# 1. Check which mpirun is being used
which mpirun

# 2. If it points to conda/miniconda (e.g., ~/miniconda3/bin/mpirun),
#    use the system mpirun explicitly instead:
/usr/bin/mpirun -n 4 ./bin/inversion -c config.toml -d data.txt

# 3. Or temporarily remove conda from PATH:
export PATH=$(echo $PATH | tr ':' '\n' | grep -v conda | tr '\n' ':')
mpirun -n 4 ./bin/inversion -c config.toml -d data.txt
```

**Verification**:

```bash
# Ensure the program links to the correct MPI libraries
ldd ./bin/inversion | grep mpi

# Should show your system MPI path (e.g., /usr/lib/x86_64-linux-gnu/libmpi.so.40)
# NOT a conda path (e.g., ~/miniconda3/lib/libmpi.so)
```

### Running with MPI

To run the inversion with multiple MPI processes:

```bash
# Basic syntax
mpirun -n <num_processes> ./bin/inversion [options]

# Example: Run with 4 processes (1 master + 3 workers)
mpirun -n 4 ./bin/inversion -c config.toml -d data.txt -o inv.h5
```

### Performance Recommendations

- **Number of processes**: Use `min(available_cores, num_init/2)` for optimal performance
  - The `num_init` parameter in `config.toml` controls the total number of parallel tasks
  - Total tasks = `num_init × num_noise`
  - Each worker process handles one initial model at a time

- **Process roles**:
  - **Rank 0 (Master)**: Coordinates tasks, aggregates results
  - **Rank 1+ (Workers)**: Perform inversion computations

- **Example**: If you have 8 CPU cores and `num_init = 50` in your config:
  ```bash
  # Use 4 processes (leaves cores for system and master coordination)
  mpirun -n 4 ./bin/inversion -c config.toml -d data.txt
  ```

## Python Libraries Installation

To properly visualize and process the forward modeling and inversion results, you'll need to install these Python libraries: `numpy`, `matplotlib`, `scipy`, and `h5py`.

Choose one of the following methods based on your preference:

#### Using Conda

```bash
# Create a new environment with Python 3.11 (optional but recommended)
conda create -n qedispinv python=3.11
conda activate qedispinv

# Install dependencies
conda install numpy matplotlib scipy h5py
```

#### Using Mamba (faster alternative to conda)

```bash
# Create a new environment with Python 3.11 (optional but recommended)
mamba create -n qedispinv python=3.11
mamba activate qedispinv

# Install dependencies
mamba install numpy matplotlib scipy h5py
```

To verify successful installation:

```bash
python -c "import numpy, matplotlib, scipy, h5py, sys; print('All dependencies installed successfully with Python', sys.version.split()[0])"
```

This will confirm both the successful installation of required libraries.


[**Back to previous page**](../README.md)