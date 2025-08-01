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

### 2. Required Libraries

#### OpenMP

- **Ubuntu/Debian**: `sudo apt-get install libomp-dev`
- **Fedora/RHEL**: `sudo dnf install libomp-devel`
- **Arch Linux**: `sudo pacman -S openmp`
- **macOS (Homebrew)**: `brew install libomp`

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
python -c "import numpy, matplotlib, scipy, h5py; print('All dependencies installed successfully with Python', sys.version.split()[0])"
```

This will confirm both the successful installation of required libraries.
