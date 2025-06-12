#ifndef MODEL_H_
#define MODEL_H_

#include <Eigen/Dense>
#include <vector>

class Model {
public:
  Model(const Eigen::Ref<const Eigen::ArrayXXd> model);
  Model() {}
  virtual Eigen::ArrayXXd generate(const std::vector<double> &z,
                                   const std::vector<double> &vs) = 0;

protected:
  Eigen::ArrayXd z2depth(const std::vector<double> &z, int nl);
  Eigen::ArrayXd z_, rho_, vs_, vp_;
};

class FixVpRho : public Model {
public:
  using Model::Model;
  Eigen::ArrayXXd generate(const std::vector<double> &z,
                           const std::vector<double> &vs) override;

private:
};

class Brocher05 : public Model {
  // Brocher, T. M. (2005). Empirical relations between elastic wavespeeds and
  // density in the Earth's crust. Bulletin of the seismological Society of
  // America, 95(6), 2081-2092.
public:
  using Model::Model;
  Eigen::ArrayXXd generate(const std::vector<double> &z,
                           const std::vector<double> &vs) override;

private:
};

class NearSurface : public Model {
  // Gardner, G. H. F., Gardner, L. W., & Gregory, A. (1974). Formation velocity
  // and density—The diagnostic basics for stratigraphic traps. Geophysics,
  // 39(6), 770-780.
public:
  using Model::Model;
  Eigen::ArrayXXd generate(const std::vector<double> &z,
                           const std::vector<double> &vs) override;

private:
  const double vp2vs_ = 1.7321;
  const double alpha_ = 0.31;
  const double beta_ = 0.25;
};

#endif