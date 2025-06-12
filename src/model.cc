#include "model.hpp"

#include <Eigen/Dense>
#include <fmt/format.h>

using namespace Eigen;

Model::Model(const Eigen::Ref<const Eigen::ArrayXXd> model)
    : z_(model.col(1)), rho_(model.col(2)), vs_(model.col(3)),
      vp_(model.col(4)) {}

Eigen::ArrayXd Model::z2depth(const std::vector<double> &z, int nl) {
  ArrayXd dep(nl);
  dep(0) = 0.0;
  for (int i = 1; i < nl; ++i) {
    dep(i) = z[i];
  }
  return dep;
}

Eigen::ArrayXXd FixVpRho::generate(const std::vector<double> &z,
                                   const std::vector<double> &vs) {
  const int nl = vs.size();
  ArrayXd dep = z2depth(z, nl);

  ArrayXd z_ref(nl);
  z_ref.head(nl - 1) = (dep.head(nl - 1) + dep.tail(nl - 1)) / 2.0;
  z_ref(nl - 1) = (dep(nl - 1) + z_(z_.rows() - 1)) / 2.0;

  ArrayXXd model(nl, 5);
  model(0, 0) = 1;
  model(0, 1) = 0.0;
  model(0, 2) = rho_(0);
  model(0, 3) = vs[0];
  model(0, 4) = vp_(0);

  int iloc = 0;
  for (int i = 1; i < nl; ++i) {
    model(i, 0) = i + 1;
    model(i, 1) = dep(i);
    model(i, 3) = vs[i];

    // linear interpolation
    while (z_(iloc) < z_ref(i)) {
      ++iloc;
    }
    double z0 = z_(iloc - 1);
    double z1 = z_(iloc);
    model(i, 2) = rho_(iloc - 1) +
                  (rho_(iloc) - rho_(iloc - 1)) * (z_ref(i) - z0) / (z1 - z0);
    model(i, 4) = vp_(iloc - 1) +
                  (vp_(iloc) - vp_(iloc - 1)) * (z_ref(i) - z0) / (z1 - z0);
  }
  return model;
}

Eigen::ArrayXXd Brocher05::generate(const std::vector<double> &z,
                                    const std::vector<double> &vs) {
  const int nl = vs.size();
  ArrayXd dep = z2depth(z, nl);

  ArrayXXd model(nl, 5);
  for (int i = 0; i < nl; ++i) {
    model(i, 0) = i + 1;
    model(i, 1) = dep(i);
    model(i, 3) = vs[i];
    double vp = 0.9409 + 2.0947 * vs[i] - 0.8206 * pow(vs[i], 2) +
                0.2683 * pow(vs[i], 3) - 0.0251 * pow(vs[i], 4);
    double rho = 1.6612 * vp - 0.4721 * pow(vp, 2) + 0.0671 * pow(vp, 3) -
                 0.0043 * pow(vp, 4) + 0.000106 * pow(vp, 5);
    model(i, 2) = rho;
    model(i, 4) = vp;
  }
  return model;
}

Eigen::ArrayXXd NearSurface::generate(const std::vector<double> &z,
                                      const std::vector<double> &vs) {
  const int nl = vs.size();
  ArrayXd dep = z2depth(z, nl);

  ArrayXXd model(nl, 5);
  for (int i = 0; i < nl; ++i) {
    model(i, 0) = i + 1;
    model(i, 1) = dep(i);
    model(i, 3) = vs[i];
    model(i, 4) = vp2vs_ * vs[i];
    model(i, 2) = alpha_ * pow(model(i, 4) * 1000.0, beta_);
  }
  return model;
}