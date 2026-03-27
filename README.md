# QEDispInv: Surface Wave Dispersion Curve Computation and Inversion Toolkit

## Overview

QEDispInv is an open-source C++ program developed for the forward modeling and inversion of surface wave dispersion curves,
accompanying this work on integrating quadratic extrema interpolation and randomized layering with multiple initial models.
This toolkit encompasses the complete workflow of surface wave analysis, from dispersion curve calculation to subsurface
structure inversion, designed to address the challenges of low-velocity layers (LVLs) and ``mode-kissing'' phenomena in
complex geological models.

## Key Features

- **Forward Modeling**: Implements the quadratic extrema interpolation method for accurate and efficient dispersion curve computation, specifically optimized to resolve mode-kissing zones in models with LVLs. This method significantly reduces root omission rates compared to conventional approaches (e.g., CPS) while achieving substantial speed improvements ($12.6\times$ for full-mode computation and $17.4\times$ for fundamental-mode calculation in benchmark tests).
- **Inversion Framework**: Incorporates wavelength-constrained layering schemes that adapt to the depth-dependent resolution of surface wave data, combined with depth-zoned empirical relationships for deriving P-wave velocity and density from shear-wave velocity. The inversion process utilizes the L-BFGS-B algorithm with **MPI-based parallel execution** across multiple initial models, with statistical aggregation of results to address non-uniqueness.
- **Parallel Computing**: MPI-parallelized inversion for distributed-memory systems, enabling efficient scaling across multiple CPU cores for large-scale parameter studies with multiple initial models and noise realizations.
- **Comprehensive Workflow**: Supports both forward computation (dispersion curve generation) and inversion (subsurface structure recovery), validated through numerical simulations (near-surface and crustal LVL models) and real data applications (Mirandola site and USArray).

## Installation

The program can be cloned from the GitHub repository:

```bash
git clone https://github.com/pan3rock/QEDispInv.git
cd QEDispInv
git submodule update --init
```

**Note**: This branch requires an OpenMPI installation for MPI-parallelized inversion. Please refer to the [`Install`](doc/INSTALL.md) document for detailed MPI setup instructions, compilation steps, and dependencies.

## Tutorial

A detailed tutorial covering step-by-step instructions for using QEDispInv, including forward
modeling setup, inversion parameter configuration, result visualization, and example walkthroughs,
is available in the [`Tutorial`](doc/TUTORIAL.md) document. This tutorial is designed to help new
users quickly get started with the toolkit and explore its full functionality.

## Usage

QEDispInv provides command-line interfaces for both forward and inversion modules:

- **Forward Modeling**: Compute dispersion curves for a given subsurface model using the quadratic extrema interpolation method.
- **Inversion**: Perform inversion of observed dispersion curves to derive shear-wave velocity profiles, incorporating the wavelength-constrained layering and empirical relationships described in this work.

Detailed usage examples and parameter specifications are available in the [`Tutorial`](doc/TUTORIAL.md).

## Citation

If you use QEDispInv in your research, please cite this work as follows:

- Pan, L., Wang, J., & Chen, X. (2026). Surface wave dispersion curve computation and inversion: A framework integrating quadratic extrema interpolation and randomized layering with multiple initial models. Bulletin of the Seismological Society of America, (Newly Published). https://doi.org/10.1785/0120250207 
- Pan, L., & Chen, X. (2025). Efficient Computation of Dispersion Curves in Low‐Velocity‐Layered Half‐Spaces. Bulletin of the Seismological Society of America, 115 (4), 1381–1391. https://doi.org/10.1785/0120250003
- Pan, L., Chen, X., Wang, J., Yang, Z., & Zhang, D. (2019). Sensitivity analysis of dispersion curves of Rayleigh waves with fundamental and higher modes. Geophysical Journal International, 216(2), 1276-1303. https://doi.org/10.1093/gji/ggy479

## Repository Link

https://github.com/pan3rock/QEDispINV

## License

This program is released under the [GNU LESSER GENERAL PUBLIC LICENSE (LGPL)](LICENSE) to promote open science and facilitate further development in the geophysical community.
