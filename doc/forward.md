## Tutorial for Forward Computation

### Calculating Dispersion Curves

```bash
# Assume your terminal is currently in the demo/lvl-l4 folder.
../../bin/forward -h
```

The help information is shown below

<img src="figures/forward-h.png" width="600" alt="forward-h">

To compute dispersion curves using the provided tool, follow these commands based on your specific needs:

- **Compute fundamental-mode dispersion curves**:  
  Run the forward modeling executable with:

  ```bash
  ../../bin/forward
  ```

  Alternatively, explicitly specify the fundamental mode (mode 0) using:

  ```bash
  ../../bin/forward -m 0
  ```

- **Compute up to the 5th mode**:  
  To include modes from the fundamental up to the 5th order, use the `-m` flag to set the maximum mode:

  ```bash
  ../../bin/forward -m 5
  ```

- **Compute Love wave dispersion curves**:  
  By default, the tool calculates Rayleigh wave dispersion. To switch to Love waves, add the `--sh` flag:

  ```bash
  ../../bin/forward --sh
  ```

- **Custom Model Files**  
  To compute dispersion curves for a specific model, use the `--model` flag followed by your model file path:

  ```bash
  ../../bin/forward --model model_data.txt
  ```

  If no model is specified, the tool will default to the model file configured in `config.toml`.

- **Matching Frequency Sampling to Existing Dispersion Data**  
  To align the frequency sampling of your computed dispersion curves with an existing dispersion file (e.g., for inversion or comparison), use the `--disp` flag:

  ```bash
  ../../bin/forward --disp disp.txt
  ```

  This ensures the output uses the same frequency points as `disp.txt`. Without this flag, frequencies are determined by the default settings in `config.toml`.

- **Computing Sensitivity Kernels**

  The tool supports calculating sensitivity kernels for P-wave velocity ($V_P$), S-wave velocity ($V_S$), and density ($\rho$). The range of modes and frequencies for kernel computation is determined by the parameters used in the prior dispersion curve calculation.

  To compute sensitivity kernels, add the `--compute_kernel` flag to your command. For example, to calculate kernels corresponding to dispersion curves of modes 0 to 3:

  ```bash
  ../../bin/forward -m 3 --compute_kernel
  ```

  All sensitivity kernel results are saved to a unified output file named `kernel.h5`.

### Visualizing Results with Python Scripts

To visualize the computed dispersion curves, follow these steps using the provided Python script:

1. **Compute the complete dispersion curves** using the forward modeling tool with a sufficiently high mode order (e.g., up to mode 1000) to capture all relevant modes:

   ```bash
   ../../bin/forward -m 1000
   ```

2. **Generate the dispersion curve plot** by running the `plot_disp.py` script, specifying the output dispersion file (e.g., `disp.txt`):
   ```bash
   ../../python/plot_disp.py disp.txt
   ```
   This will produce a visualization of the dispersion curves, as shown in the figure below.

<img src="figures/disp_lvl4.jpg" width="600" alt="disp_lvl4">

The `plot_disp.py` script includes additional options for customizing the plot (e.g., scaling, annotations, or comparison). To explore these features, check the help documentation with:

```bash
../../python/plot_disp.py -h
```

Further details on each option are available in the help output, encouraging users to experiment with configurations tailored to their analysis needs.

To examine the computed sensitivity kernels using the Python visualization tool, follow these steps:

1. **Calculate the sensitivity kernels** for modes 0 to 2 (3 modes total) using the forward modeling tool with the `--compute_kernel` flag:

   ```bash
   ../../bin/forward -m 2 --compute_kernel
   ```

   This generates the `kernel.h5` file containing the kernels for $V_P$, $V_S$, and density ($\rho$).

2. **Visualize specific kernels** using the `plot_kernel.py` script. Use the `-m` flag to specify the mode, and `--comp` to select the parameter ($V_S$ is default if not specified):

   - View the $V_S$ kernel for the fundamental mode (mode 0):

     ```bash
     ../../python/plot_kernel.py kernel.h5 -m 0
     ```

   - View the $V_S$ kernel for the 1st mode:

     ```bash
     ../../python/plot_kernel.py kernel.h5 -m 1
     ```

   - View the density ($\rho$) kernel for the fundamental mode:

     ```bash
     ../../python/plot_kernel.py kernel.h5 -m 0 --comp rho
     ```

   - View the $V_P$ kernel for the fundamental mode:
     ```bash
     ../../python/plot_kernel.py kernel.h5 -m 0 --comp vp
     ```

3. **Additional visualization options**:

   - Add a colorbar to the kernel plot with `--show_cb`:

     ```bash
     ../../python/plot_kernel.py kernel.h5 -m 0 --show_cb
     ```

   - Overlay the corresponding dispersion curve on the kernel plot using `--plot_disp`:
     ```bash
     ../../python/plot_kernel.py kernel.h5 -m 0 --plot_disp
     ```

The `plot_kernel.py` script includes further customization features not covered here. To explore the full range of options (such as scaling adjustments or plot styling), refer to the built-in help documentation:

```bash
../../python/plot_kernel.py -h
```

This allows users to tailor visualizations to their specific analytical requirements through self-guided exploration.

By adjusting the calculation frequencies and modifying the parameters of the plotting script,
We can combine these elements to generate a figure as shown below.

<img src="figures/kernel.jpg" width="600" alt="kernel_lvl4">
