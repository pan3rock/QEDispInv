# Demo

To facilitate users in quickly understanding and utilizing the code, three examples are provided
below, each demonstrating the application of the framework in different scenarios based on the
content of the manuscript.

## Example 1: [Model 1](../demo/lvl-l4) (Four-Layer Model with a Low-Velocity Layer)

Model 1 is a four-layer numerical model designed to address the ``mode-kissing'' phenomenon, where adjacent dispersion curve branches converge closely in media containing low-velocity layers (LVLs). This model, with parameters detailed in Table 1, helps illustrate the limitations of conventional root-finding methods in detecting critical roots at mode-kissing zones. By applying the proposed quadratic extrema interpolation method, supplementary sampling points are strategically added based on the dispersion function’s morphology, resolving previously missed roots and ensuring complete dispersion curve computation. This example highlights the method’s efficiency in handling LVL-induced complexities, which is crucial for accurate subsequent inversion.

## Example 2: [Near-Surface Numerical Inversion](../demo/syn-nearsurface/)

This near-surface example uses a five-layer theoretical model (parameters in Table 2) to validate the inversion workflow. Simulated Rayleigh wave data are generated via the CPS package, processed using the Frequency-Bessel (F-J) transform for dispersion energy imaging, and multi-mode dispersion curves are extracted. The inversion employs the L-BFGS-B algorithm with wavelength-constrained layering, focusing on shear-wave velocity ($V_S$) optimization. Results show robust recovery of the target model, including accurate resolution of a shallow low-velocity layer (2–6 m depth) and consistent velocity trends. This example demonstrates the framework’s reliability in near-surface applications such as engineering site investigation.

## Example 3: [Crustal LVL and Uppermost Mantle Inversion](../demo/syn-crustmantle/)

This example validates the framework’s performance in deeper structures using a crustal model with a low-velocity layer (LVL), representative of tectonically active regions or rift zones. Theoretical dispersion curves are modified to include mode-dependent noise and frequency constraints, mimicking real-world data. Inversion with 100 randomized initial models (using Brocher’s empirical relationships for $V_P$ and density) successfully resolves key features, including the crustal LVL and uppermost mantle velocity trends. This example confirms the method’s ability to handle deep subsurface complexities, supporting crustal and lithospheric studies.

# Tutorial

## [Secular Function Computation](./secfunc.md)

## [Forward Computation](./forward.md)

## [Inversion](./inversion.md)
