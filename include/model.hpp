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

protected:
  Eigen::ArrayXd z2depth(const Eigen::ArrayXd &z, int nl);
  Eigen::ArrayXd z_, rho_, vs_, vp_;
};

class FixVpRho : public Vs2Model {
public:
  using Vs2Model::Vs2Model;
  Eigen::ArrayXXd generate(const Eigen::ArrayXd &z,
                           const Eigen::ArrayXd &vs) override;
  std::map<std::string, Eigen::ArrayXd>
  derivative(const Eigen::ArrayXd &vs) override;

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

class NearSurface : public Vs2Model {
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

Eigen::ArrayXd generate_random_depth(int N, double zmax, double min_gap);

#endif