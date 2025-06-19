#include "model.hpp"

#include <Eigen/Dense>
#include <fmt/format.h>
#include <map>
#include <random>

using namespace Eigen;
using Dict = std::map<std::string, Eigen::ArrayXd>;

Vs2Model::Vs2Model(const Eigen::Ref<const Eigen::ArrayXXd> model)
    : z_(model.col(1)), rho_(model.col(2)), vs_(model.col(3)),
      vp_(model.col(4)) {}

Eigen::ArrayXd Vs2Model::z2depth(const Eigen::ArrayXd &z, int nl) {
  ArrayXd dep(nl);
  dep(0) = 0.0;
  for (int i = 1; i < nl; ++i) {
    dep(i) = z(i);
  }
  return dep;
}

Eigen::ArrayXXd FixVpRho::generate(const Eigen::ArrayXd &z,
                                   const Eigen::ArrayXd &vs) {
  const int nl = vs.size();
  ArrayXd dep = z2depth(z, nl);

  ArrayXd z_ref(nl);
  z_ref.head(nl - 1) = (dep.head(nl - 1) + dep.tail(nl - 1)) / 2.0;
  z_ref(nl - 1) = (dep(nl - 1) + z_(z_.rows() - 1)) / 2.0;

  ArrayXXd model(nl, 5);
  model(0, 0) = 1;
  model(0, 1) = 0.0;
  model(0, 2) = rho_(0);
  model(0, 3) = vs(0);
  model(0, 4) = vp_(0);

  int iloc = 0;
  for (int i = 1; i < nl; ++i) {
    model(i, 0) = i + 1;
    model(i, 1) = dep(i);
    model(i, 3) = vs(i);

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

Dict FixVpRho::derivative([[maybe_unused]] const Eigen::ArrayXd &vs) {
  int nl = vs.rows();
  ArrayXd drho = ArrayXd::Zero(nl);
  ArrayXd dvp = ArrayXd::Zero(nl);
  ArrayXd dvs = ArrayXd::Ones(nl);

  Dict res;
  res["rho"] = drho;
  res["vp"] = dvp;
  res["vs"] = dvs;
  return res;
}

Eigen::ArrayXXd Brocher05::generate(const Eigen::ArrayXd &z,
                                    const Eigen::ArrayXd &vs) {
  const int nl = vs.size();
  ArrayXd dep = z2depth(z, nl);

  ArrayXXd model(nl, 5);
  for (int i = 0; i < nl; ++i) {
    model(i, 0) = i + 1;
    model(i, 1) = dep(i);
    model(i, 3) = vs(i);
    double vp = 0.9409 + 2.0947 * vs[i] - 0.8206 * pow(vs[i], 2) +
                0.2683 * pow(vs[i], 3) - 0.0251 * pow(vs[i], 4);
    double rho = 1.6612 * vp - 0.4721 * pow(vp, 2) + 0.0671 * pow(vp, 3) -
                 0.0043 * pow(vp, 4) + 0.000106 * pow(vp, 5);
    model(i, 2) = rho;
    model(i, 4) = vp;
  }
  return model;
}

Dict Brocher05::derivative(const Eigen::ArrayXd &vs) {
  ArrayXd dvp = 2.0947 - 0.8206 * vs * 2 + 0.2683 * vs.pow(2) * 3 -
                0.0251 * vs.pow(3) * 4;
  ArrayXd vp = 0.9409 + 2.0947 * vs - 0.8206 * vs.pow(2) + 0.2683 * vs.pow(3) -
               0.0251 * vs.pow(4);
  ArrayXd drho = 1.6612 - 0.4721 * vp * 2 + 0.0671 * vp.pow(2) * 3 -
                 0.0043 * vp.pow(3) * 4 + 0.000106 * vp.pow(4) * 5;
  drho *= dvp;
  ArrayXd dvs = ArrayXd::Ones(vs.rows());

  Dict res;
  res["rho"] = drho;
  res["vp"] = dvp;
  res["vs"] = dvs;
  return res;
}

Eigen::ArrayXXd NearSurface::generate(const Eigen::ArrayXd &z,
                                      const Eigen::ArrayXd &vs) {
  const int nl = vs.size();
  ArrayXd dep = z2depth(z, nl);

  ArrayXXd model(nl, 5);
  for (int i = 0; i < nl; ++i) {
    model(i, 0) = i + 1;
    model(i, 1) = dep(i);
    model(i, 3) = vs(i);
    model(i, 4) = vp2vs_ * vs[i];
    model(i, 2) = alpha_ * pow(model(i, 4) * 1000.0, beta_);
  }
  return model;
}

Dict NearSurface::derivative(const Eigen::ArrayXd &vs) {
  ArrayXd dvs = ArrayXd::Ones(vs.rows());
  ArrayXd dvp = vp2vs_ * dvs;
  ArrayXd vp = vp2vs_ * vs;
  ArrayXd drho = alpha_ * beta_ * pow(1000.0 * vp, beta_ - 1.0) * 1000.0 * dvp;

  Dict res;
  res["rho"] = drho;
  res["vp"] = dvp;
  res["vs"] = dvs;
  return res;
}

Eigen::ArrayXd generate_random_depth(int nl, double zmax, double min_gap) {
  int N = nl - 1;
  if (N < 2) {
    throw std::invalid_argument("N must be at least 2");
  }
  if (min_gap <= 0) {
    throw std::invalid_argument("min_gap must be positive");
  }
  if (min_gap * (N - 1) >= 1.0) {
    throw std::invalid_argument("min_gap too large for given N");
  }

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  double total_min_gap = min_gap * (N - 1);
  double leftover = 1.0 - total_min_gap - 1e-9;

  std::vector<double> increments(N - 1);
  double sum = 0.0;
  for (int i = 0; i < N - 1; ++i) {
    increments[i] = dist(rng);
    sum += increments[i];
  }

  if (sum > 0) {
    for (int i = 0; i < N - 1; ++i) {
      increments[i] = (increments[i] / sum) * leftover;
    }
  } else {
    double avg_increment = leftover / (N - 1);
    for (int i = 0; i < N - 1; ++i) {
      increments[i] = avg_increment;
    }
  }

  std::vector<double> sequence;
  sequence.reserve(N);
  sequence.push_back(0.0);

  double current = 0.0;
  for (int i = 0; i < N - 1; ++i) {
    current += min_gap + increments[i];
    sequence.push_back(current);
  }

  std::vector<double> depth;
  for (int i = 0; i < N; ++i) {
    depth.push_back(sequence[i] * zmax);
  }
  depth.push_back(zmax);

  ArrayXd ret = Map<ArrayXd, Unaligned>(depth.data(), depth.size());

  return ret;
}