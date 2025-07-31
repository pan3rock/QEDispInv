## Tutorial for Inversion

### Understanding the Configuration File (`config.toml`)

The inversion workflow is guided by a configuration file named config.toml, which centralizes key parameters to streamline the inversion process. Below is an explanation of its core sections and settings, as illustrated in the example configuration (refer to the figure for visual reference):

<img src="figures/config.png" width="600" alt="config">

By adjusting `config.toml`, users can tailor the inversion to their specific dataset (e.g., near-surface vs. crustal scales) and balance computational efficiency with solution robustness.

The `[inversion]` block in `config.toml` (as shown in the figure) contains critical parameters that govern the inversion’s mathematical and statistical behavior. Here’s a breakdown of each setting:

1. **Model Conversion (`vs2model`)**

   - `vs2model`: Specifies how shear-wave velocity ($V_S$) is converted to other parameters (e.g., $V_P$ and density).
   - **Supported Modes**:
     - `fixvprho`: Fixes \(V_P\) and density to values from the reference model (`model_ref`), useful when these parameters are well-constrained.
     - `nearsurface`: Optimized for shallow models (meters to hundreds of meters), requiring the `vp2vs` parameter to define \(V_P/V_S\) ratios.
     - `gardner`: Applies Gardner’s relation for deeper models (kilometers), where density is derived from \(V_P\) (\(\rho = 0.31 V_P^{0.25}\)).
     - `brocher05`: Uses Brocher’s (2005) empirical relations for crustal to upper-mantle models (tens of kilometers), linking \(V_P\), \(V_S\), and density.

2. **Reference Model (`model_ref`)**

   - `model_ref = "mref.txt"`: Path to a text file containing a starting/reference model (e.g., a 1D layered model with depth, \(V_S\), \(V_P\), and density). This guides the inversion toward geologically reasonable solutions.

3. **Velocity Uncertainty (`vs_width`)**

   - `vs_width = 0.7`: Defines the allowable range for updating $V_S$ during inversion, calculated as $[V_{S,\text{ref}} - v_{\text{width}}/2, V_{S,\text{ref}} + v_{\text{width}}/2]$. Here, $V\_{S,\text{ref}}$ is the shear-wave velocity from the reference model (`model_ref`), linearly interpolated to match the depth nodes of the current initial model. This parameter restricts $V_S$ adjustments to a geophysically plausible range, preventing unrealistic deviations from the reference.

4. **Regularization Strength (`lambda`)**

   - `lambda = 1.0e-2`: Controls the tradeoff between data fit and model smoothness. Smaller values prioritize fitting data (risk of overfitting), while larger values enforce smoother models.

5. **Regularization Type (`reg_type`)**

   - `reg_type = 2`: Specifies the type of regularization (e.g., 1 = first-order Tikhonov regularization, 2 = An's (2020) adaptive first-order Tikhonov regularization), which penalizes rapid changes in model parameters.

6. **Initial Model Exploration (`num_init`, `num_noise`)**

   - `num_init = 100`: Number of randomized initial models to explore, reducing dependence on the starting model.
   - `num_noise = 1`: Number of noise realizations added to data for uncertainty analysis.

7. **Randomization (`rand_depth`, `rand_vs`)**

   - `rand_depth = true`: Enables random variation of layer depths in initial models.
   - `rand_vs = false`: Disables random variation of $V_S$ in initial models (useful if $V_S$ is well-constrained).

8. **Depth Constraints (`zmax`)**

   - `zmax = 0.025`: Serves as the lower bound for depth statistics in the inversion results. The actual inversion depth range is determined by the fundamental-mode dispersion curve (specifically, $\lambda_{\text{max}}/2$, where $\lambda_{\text{max}}$ is the maximum wavelength of the fundamental mode). Since the workflow uses a multi-initial-model inversion (generating multiple solution models), statistics like median, 10th percentile (P10), and 90th percentile (P90) are used to analyze results. `zmax` defines the deepest depth included in these statistical calculations, ensuring consistent depth bounds for comparing model uncertainties.

9. **layering Scaling Bounds (`r0`, `rmin`, `rmax`)**

   - `r0 = 0.5`, `rmin = 1.0`, `rmax = 1.5`: These parameters constrain the scaling factors governing the layering scheme of the inversion model, adapting to surface wave resolution limits:
     - `r0` corresponds to the reference value of `a`, a scaling factor for the shallowest interface ($d_1$). The first layer depth is defined as $d_1 = (\lambda_{\text{min}} / 2) \times a$, where $\lambda\_{\text{min}}$ is the minimum wavelength of the fundamental-mode dispersion curve, and `a` is constrained to $0 < a \leq 1$ (with `r0` as a baseline for fine-tuning near-surface resolution).
     - `rmin` and `rmax` bound `b`, the multiplicative factor for deeper layer thicknesses (below $d_2$). Each deeper layer's thickness is the product of the overlying layer’s thickness and a `b` value randomly selected within $(b_{\text{min}}, b_{\text{max}})$, with $b_{\text{min}}$ ensuring thicknesses increase gradually with depth—matching the reduced sensitivity of surface waves to deeper structures.

<img src="figures/layer.jpg" width="600" alt="layer">

10. **Parameter Ratio (`vs2vp`)**

    - `vs2vp = 3.0`: Fixed ratio of \(V_S\) to \(V_P\) (e.g., \(V_P = 3.0 V_S\)), simplifying inversion by reducing free parameters.

11. **Data Weights (`weight`)**
    - `weight = [4.0, 1.0, 1.0, 1.0]`: Weights applied to different mode (e.g., fundamental vs. higher modes), prioritizing high-confidence data.

### Further Explanation of `vs2model`

The parameterization strategy centers on S-wave velocity ($V_S$) as the core inversion parameter, leveraging its pivotal role in governing surface-wave dispersion. This design not only mitigates solution non-uniqueness but also preserves robust geophysical interpretability. In this framework, P-wave velocity ($V_P$) and density ($\rho$) are derived through empirical correlations with $V_S$---an approach that ensures petrophysical consistency across the model space while effectively reducing the dimensionality of the inverse problem. This section focuses on introducing empirical relationships under three different depth intervals, which are defined by inversion targets: the shallow subsurface, the middle layer, and the deep layer (encompassing the crust and uppermost mantle), to illustrate how $V_P$ and $\rho$ are specifically obtained via $V_S$-based empirical relationships.

1. **Empirical Relationships in the Shallow Subsurface (nearsurface)**

The shallow subsurface (typically $0$ to a few hundred meters) comprises unconsolidated sediments, weathered rocks, and shallow bedrock. Its empirical relationships between $V_S$-$V_P$ and $V_S$-$\rho$ are primarily governed by water saturation, lithology, and compaction, requiring site-specific calibration due to high sensitivity to near-surface conditions.

In the shallow subsurface, $V_P$ is derived from $V_S$ using $V_P/V_S$ ratio-based empirical relationships, and these ratios vary significantly with the saturation state of the materials and their lithological characteristics. For water-saturated soils, which are the most common in many shallow geological settings, the $V_P/V_S$ ratio ranges from approximately $3$ to over $10$. A common approximation in such cases is $V_P \approx 10\times V_S$. Here, $V_P$ is dominated by the pore water, which typically has a velocity of around $1.5$ km/s. For unsaturated or partially saturated soils, the $V_P/V_S$ ratio is lower, generally between $1.4$ and $3.3$. In this scenario, both $V_P$ and $V_S$ are controlled by the soil skeleton rather than pore fluids. For consolidated bedrock in the shallow subsurface, the $V_P/V_S$ ratio stabilizes at $1.6$-$2.0$, corresponding to a Poisson's ratio ($\nu$) of approximately $0.25$-$0.3$. Due to the consolidated nature of the rock, its velocity relationship is less sensitive to pore fluids compared to unconsolidated sediments.

In the shallow subsurface, distinct material types exhibit characteristic ranges of shear-wave velocity ($V_S$) and corresponding density, as summarized in the table. These ranges reflect the progressive increase in both stiffness and compactness from unconsolidated sediments to competent bedrock.

| **Material Type**  | **Vs Range (km/s)** | **Density (g/cm³)** |
| ------------------ | ------------------- | ------------------- |
| Soft clay and silt | 0.08-0.15           | 1.60-1.80           |
| Silty sand         | 0.15-0.25           | 1.70-1.90           |
| Sand               | 0.20-0.40           | 1.80-2.00           |
| Gravel             | 0.25-0.50           | 1.90-2.20           |
| Weathered rock     | 0.50-0.80           | 2.10-2.50           |
| Competent rock     | >0.80               | 2.50-2.80           |

**Table: Vs and Density Ranges for Different Shallow Materials**

To derive a quantitative relationship between $V_S$ and density ($\rho$) for shallow subsurface materials, we performed a quadratic function fitting using the observed velocity-density pairs across these lithological categories. This approach accounts for the non-linear nature of the relationship, particularly evident in the transition from unconsolidated sediments to weathered rock. The resulting empirical formula is expressed as:

$$\rho = -0.22374079 V_S^2 + 1.32248261 V_S + 1.54840433,$$

For competent rock ($V_S > 0.8$ km/s), the formula yields values greater than $2.38$ g/cm³, though lithology-specific adjustments are recommended due to potential deviations from the general trend. As with all empirical relationships, site-specific calibration using borehole density measurements is advised to improve accuracy in regional applications.

2. **Empirical Relationships in the Middle Layer (gardner)**

The middle layer, which is usually between a few hundred meters and several kilometers deep, consists of more consolidated rocks with lower porosity compared to the shallow subsurface. The rock types here are more stable, and the physical properties are mainly controlled by lithology and metamorphic grade. This stability leads to more consistent empirical relationships between $V_S$ and $V_P$, and between $V_S$ and $\rho$.

In terms of $V_P$, for crystalline rocks in the middle layer, such as granite and gneiss, empirical relationships have been developed based on extensive laboratory measurements and field data. A commonly used linear relationship is adopted here, specifically $V_P = 1.732\times V_S$. This relationship is applicable to most consolidated rocks in the middle layer, such as granite, gneiss, and andesite.

For the derivation of density ($\rho$) from $V_S$ in the middle layer, the Gardner formula is employed. The Gardner formula is a widely recognized empirical relationship in geophysics that describes the correlation between rock density and P-wave velocity. The corrected form used here is $\rho = a\times(1000V_P)^b$, where $\rho$ is the density in g/cm³, $V_P$ is the P-wave velocity in km/s, and the constants are defined as $a = 0.31$ and $b = 0.25$.

Since $V_P$ in the middle layer is derived from the linear relationship $V_P = 1.732\times V_S$, we substitute this expression into the Gardner formula to establish a direct link between $V_S$ and $\rho$. This substitution yields:

$$\rho = 0.31\times(1732\times V_S)^{0.25}$$

Notably, while this formula is applicable to most common lithologies in the middle layer (e.g., granite, metamorphic rocks), adjustments to the constants $a$ and $b$ may be required for lithologies with exceptional properties—such as high-porosity volcanic rocks or ultramafic intrusions—to account for their unique density-velocity relationships. Site-specific calibration with borehole density measurements remains advisable for critical applications.

3. **Empirical Relationships in the Deep Layer (brocher05)**

The deep layer, spanning the crust and the uppermost mantle, is characterized by extreme pressure and temperature conditions. The rocks here, including igneous, metamorphic, and residual sedimentary varieties, are strongly influenced by tectonic activity and magma intrusion. In this context, the empirical relationships between $V_S$ and $V_P$, as well as between $V_S$ and $\rho$, are adapted to these unique environmental factors. For this depth range, the empirical relationships proposed by [Brocher, 2005] are adopted, which have been widely validated in studies of crustal and uppermost mantle structures.

Brocher derived a fifth-order polynomial relationship for $V_S$ (km/s) as a function of $V_P$ (km/s), valid for the common crustal $V_P$ range of approximately $1.5$ km/s to $8.5$ km/s.

$$V_P = 0.9409 + 2.0947 V_S - 0.8206 V_S^2 + 0.2683 V_S^3 - 0.0251 V_S^4$$

The relationship was developed using extensive laboratory and in situ measurements for diverse crustal rock types.

He established a polynomial relationship for density $\rho$ (g/cm³) as a function of $V_P$ (km/s), also applicable to the common crustal $V_P$ range ($1.5$ - $8.5$ km/s).

$$\rho = 1.6612 V_P - 0.4721 V_P^2 + 0.0671 V_P^3 - 0.0043 V_P^4 + 0.000106 V_P^5$$

While density exhibits more scatter relative to $V_P$ than $V_S$ does, the equation provides a robust estimate.

Exercise significant caution when applying Brocher's empirical relations for estimating S-wave velocity ($V_S$) or density ($\rho$) from P-wave velocity ($V_P$) to rocks within the uppermost mantle. These relations were explicitly developed and calibrated using data primarily from crustal rocks within a specific $V_P$ range of $1.5$ to $8.5$ km/s. The mineralogical composition and physical state of ultramafic uppermost mantle rocks (dominantly peridotite) differ fundamentally from crustal lithologies. Consequently, applying the crustal relations to mantle $V_P$ values (typically > $8.0$ km/s) represents an extrapolation beyond the valid range defined in the original study. While the impact on density estimates might be somewhat mitigated in surface wave dispersion curve modeling due to the relatively lower sensitivity of dispersion to density compared to shear-wave velocity, significant deviations remain likely, especially for density itself. These deviations could still compromise interpretations where accurate density is crucial, such as in gravity modeling, lithospheric studies, or combined seismic-gravity inversions. Therefore, for the shallow upper mantle, seek and utilize empirical relations specifically derived for mantle conditions or clearly acknowledge the extrapolation and its potential impact.
