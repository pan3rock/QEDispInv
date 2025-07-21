#ifndef MODEL_H_
#define MODEL_H_

#include <Eigen/Dense>
#include <map>

class Vs2Model {
public:
  Vs2Model(const Eigen::Ref<const Eigen::ArrayXXd> model);
  Vs2Model() {}
  virtual Eigen::ArrayXXd generate(const Eigen::ArrayXd &z,
                                   const Eigen::ArrayXd &vs) = 0;
  virtual std::map<std::string, Eigen::ArrayXd>
  derivative(const Eigen::ArrayXd &vs) = 0;
  Eigen::ArrayXd interp_vs(const Eigen::ArrayXd &z);
  virtual void get_vs_limits(const Eigen::ArrayXd &z, double vs_width,
                             Eigen::ArrayXd &vs_ref, Eigen::ArrayXd &lb,
                             Eigen::ArrayXd &ub);
  virtual void set_param(const std::vector<double> &) {}

protected:
  Eigen::ArrayXd z2interpdepth(const Eigen::ArrayXd &z);
  Eigen::ArrayXd z_, rho_, vs_, vp_;
};

class FixVpRho : public Vs2Model {
public:
  using Vs2Model::Vs2Model;
  Eigen::ArrayXXd generate(const Eigen::ArrayXd &z,
                           const Eigen::ArrayXd &vs) override;
  std::map<std::string, Eigen::ArrayXd>
  derivative(const Eigen::ArrayXd &vs) override;
  void get_vs_limits(const Eigen::ArrayXd &z, double vs_width,
                     Eigen::ArrayXd &vs_ref, Eigen::ArrayXd &lb,
                     Eigen::ArrayXd &ub) override;

private:
};

class Brocher05 : public Vs2Model {
  // Brocher, T. M. (2005). Empirical relations between elastic wavespeeds and
  // density in the Earth's crust. Bulletin of the seismological Society of
  // America, 95(6), 2081-2092.
public:
  using Vs2Model::Vs2Model;
  Eigen::ArrayXXd generate(const Eigen::ArrayXd &z,
                           const Eigen::ArrayXd &vs) override;
  std::map<std::string, Eigen::ArrayXd>
  derivative(const Eigen::ArrayXd &vs) override;

private:
};

class Gardner : public Vs2Model {
  // Gardner, G. H. F., Gardner, L. W., & Gregory, A. (1974). Formation velocity
  // and density—The diagnostic basics for stratigraphic traps. Geophysics,
  // 39(6), 770-780.
public:
  using Vs2Model::Vs2Model;
  Eigen::ArrayXXd generate(const Eigen::ArrayXd &z,
                           const Eigen::ArrayXd &vs) override;
  std::map<std::string, Eigen::ArrayXd>
  derivative(const Eigen::ArrayXd &vs) override;

private:
  const double vp2vs_ = 1.7321;
  const double alpha_ = 0.31;
  const double beta_ = 0.25;
};

class NearSurface : public Vs2Model {
  /**
   * Shallow Surface Wave Exploration (Vs-Vp Empirical Relationships)
   *
   * Core Relationships:
   * 1. Water-Saturated Soils (most common):
   *    - Vp/Vs ≈ 3-10+  (Typically Vs ≈ 0.1*Vp)
   *    - Vp dominated by pore water (∼1500 m/s)
   *    - Vs reflects soil skeleton strength
   *
   * 2. Unsaturated/Partially-Saturated Soils:
   *    - Vp/Vs ≈ 1.4-3  (Typically Vs ≈ 0.3-0.7*Vp)
   *    - Both controlled by soil skeleton
   *
   * 3. Consolidated Bedrock:
   *    - Vp/Vs ≈ 1.6-2.0  (Poisson's ratio ν ≈ 0.25-0.3)
   *
   * Based on field and laboratory measurements:
   *   Material Type        | Dry         | Partial Sat | Saturated
   *   ------------------------------------------------------------
   *   Soft clay/silt       | 2.0-2.5     | 2.8-3.5     | 3.5-5.0
   *   Sand                 | 1.8-2.2     | 2.2-2.8     | 3.0-4.0
   *   Gravel               | 1.7-2.0     | 2.0-2.5     | 2.5-3.0
   *   Weathered rock       | 1.7-1.9     | 1.8-2.0     | 1.8-2.2
   *   Solid rock           | 1.6-1.8     | -           | 1.7-1.9
   *
   * Critical Notes:
   * 1. Water saturation is the dominant factor
   * 2. Vs is the primary indicator of geomechanical properties
   * 3. Vp insensitive to skeleton changes in saturated zones
   * 4. Avoid universal formulas (site variations >30%)
   */

  /*
   * Shallow Surface Wave Exploration (Vs-Density Empirical Relationships)
   *   rho = a * vs^2 + b * vs + c
   *   Where:
   *     rho = Density in g/cm^3
   *     vs = Shear-wave velocity in km/s
   *     a = -0.22374079
   *     b = 1.32248261
   *     c = 1.54840433
   *
   * Validation Range (vs in km/s, rho in g/cm^3):
   *   Material Type         Vs Range      rho Range       Formula Fit
   *   -------------------------------------------------------------
   *   Soft clay/silt        0.08-0.15     1.60-1.80     1.61-1.78
   *   Silty sand            0.15-0.25     1.70-1.90     1.77-1.92
   *   Sand                  0.20-0.40     1.80-2.00     1.90-1.98
   *   Gravel                0.25-0.50     1.90-2.20     1.97-2.13
   *   Weathered rock        0.50-0.80     2.10-2.50     2.12-2.38
   *   Competent rock        >0.80         2.50-2.80     >2.38 (caution)
   *
   * Usage Recommendations:
   *   1. Preferred for unconsolidated to semi-consolidated materials (Vs < 0.8
   * km/s)
   *   2. For competent rock (Vs > 0.8 km/s), combine with lithology-specific
   * adjustments
   *   3. Always validate with site-specific borehole data when available
   */

public:
  using Vs2Model::Vs2Model;
  Eigen::ArrayXXd generate(const Eigen::ArrayXd &z,
                           const Eigen::ArrayXd &vs) override;
  std::map<std::string, Eigen::ArrayXd>
  derivative(const Eigen::ArrayXd &vs) override;

  void set_param(const std::vector<double> &param) override;

private:
  double vp2vs_;
  const double a_ = -0.22374079;
  const double b_ = 1.32248261;
  const double c_ = 1.54840433;
};

Eigen::ArrayXd generate_random_depth(int N, double zmax, double min_gap);

Eigen::ArrayXd generate_depth_by_layer_ratio(double lmin, double lmax,
                                             double r0, double rmin,
                                             double rmax, double zmax);

#endif