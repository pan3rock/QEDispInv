### Further Explanation of `vs2model`

The parameterization strategy centers on S-wave velocity ($V_S$) as the core inversion parameter, leveraging its pivotal role in governing surface-wave dispersion. This design not only mitigates solution non-uniqueness but also preserves robust geophysical interpretability. In this framework, P-wave velocity ($V_P$) and density ($\rho$) are derived through empirical correlations with $V_S$-an approach that ensures petrophysical consistency across the model space while effectively reducing the dimensionality of the inverse problem. This section focuses on introducing empirical relationships under three different depth intervals, which are defined by inversion targets: the shallow subsurface, the middle layer, and the deep layer (encompassing the crust and uppermost mantle), to illustrate how $V_P$ and $\rho$ are specifically obtained via $V_S$-based empirical relationships.

1. **Empirical Relationships in the Shallow Subsurface (nearsurface)**

The shallow subsurface (typically $0$ to a few hundred meters) comprises unconsolidated sediments, weathered rocks, and shallow bedrock. Its empirical relationships between $V_S-V_P$ and $V_S-\rho$ are primarily governed by water saturation, lithology, and compaction, requiring site-specific calibration due to high sensitivity to near-surface conditions.

In the shallow subsurface, $V_P$ is derived from $V_S$ using $V_P/V_S$ ratio-based empirical relationships, and these ratios vary significantly with the saturation state of the materials and their lithological characteristics. For water-saturated soils, which are the most common in many shallow geological settings, the $V_P/V_S$ ratio ranges from approximately $3$ to over $10$. A common approximation in such cases is $V_P \approx 10\times V_S$. Here, $V_P$ is dominated by the pore water, which typically has a velocity of around $1.5$ km/s. For unsaturated or partially saturated soils, the $V_P/V_S$ ratio is lower, generally between $1.4$ and $3.3$. In this scenario, both $V_P$ and $V_S$ are controlled by the soil skeleton rather than pore fluids. For consolidated bedrock in the shallow subsurface, the $V_P/V_S$ ratio stabilizes at $1.6-2.0$, corresponding to a Poisson's ratio ($\nu$) of approximately $0.25-0.3$. Due to the consolidated nature of the rock, its velocity relationship is less sensitive to pore fluids compared to unconsolidated sediments.

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

[**Back to previous page**](./inversion.md)