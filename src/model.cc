#include "model.hpp"

#include <Eigen/Dense>
#include <fmt/format.h>
#include <iostream>
#include <map>
#include <random>

using namespace Eigen;
using Dict = std::map<std::string, Eigen::ArrayXd>;

namespace {
Eigen::ArrayXd interp1d(const Eigen::ArrayXd &x_old,
                        const Eigen::ArrayXd &y_old,
                        const Eigen::ArrayXd &x_new) {
  const Eigen::Index n_old = x_old.size();
  assert(n_old >= 2 && "x_old must contain at least 2 points");
  assert(y_old.size() == n_old && "x_old and y_old must be the same size");

  assert(((x_old.tail(n_old - 1) - x_old.head(n_old - 1)).minCoeff() > 0) &&
         "x_old must be strictly increasing");

  Eigen::ArrayXd y_new(x_new.size());

  const double x_min = x_old.minCoeff();
  const double x_max = x_old.maxCoeff();

  Eigen::ArrayXd slopes(n_old - 1);
  for (Eigen::Index i = 0; i < n_old - 1; ++i) {
    slopes[i] = (y_old[i + 1] - y_old[i]) / (x_old[i + 1] - x_old[i]);
  }

  for (Eigen::Index i = 0; i < x_new.size(); ++i) {
    const double x = x_new[i];

    // extrapolate
    if (x <= x_min) {
      y_new[i] = y_old[0];
    } else if (x >= x_max) {
      y_new[i] = y_old[n_old - 1];
    }
    // interpolate
    else {
      Eigen::Index pos =
          std::upper_bound(x_old.data(), x_old.data() + n_old, x) -
          x_old.data() - 1;
      y_new[i] = y_old[pos] + (x - x_old[pos]) * slopes[pos];
    }
  }

  return y_new;
}
} // namespace

Vs2Model::Vs2Model(const Eigen::Ref<const Eigen::ArrayXXd> model)
    : z_(model.col(1)), rho_(model.col(2)), vs_(model.col(3)),
      vp_(model.col(4)) {}

Eigen::ArrayXd Vs2Model::z2interpdepth(const Eigen::ArrayXd &z) {
  const int nl = z.rows();
  ArrayXd dep(nl);
  dep(0) = 0.0;
  for (int i = 1; i < nl - 1; ++i) {
    dep(i) = (z(i) + z(i + 1)) / 2.0;
  }
  dep(nl - 1) = z(nl - 1);
  return dep;
}

Eigen::ArrayXd Vs2Model::interp_vs(const Eigen::ArrayXd &z) {
  ArrayXd dep = z2interpdepth(z);
  return interp1d(z_, vs_, dep);
}

void Vs2Model::get_vs_limits(const Eigen::ArrayXd &z, double vs_width,
                             Eigen::ArrayXd &vs_ref, Eigen::ArrayXd &lb,
                             Eigen::ArrayXd &ub) {
  vs_ref = interp_vs(z);
  ArrayXd vs_min = vs_ref - vs_width / 2.0;
  ArrayXd vs_max = vs_ref + vs_width / 2.0;
  for (int i = 0; i < vs_min.rows(); ++i) {
    vs_min(i) = std::max(0.01, vs_min(i));
  }
  lb = vs_min;
  ub = vs_max;
}

Eigen::ArrayXXd FixVpRho::generate(const Eigen::ArrayXd &z,
                                   const Eigen::ArrayXd &vs) {
  ArrayXd dep = z2interpdepth(z);
  const int nl = vs.size();
  ArrayXXd model(nl, 5);
  model.col(0) = ArrayXd::LinSpaced(nl, 1, nl);
  model.col(1) = z;
  model.col(2) = interp1d(z_, rho_, dep);
  model.col(3) = vs;
  model.col(4) = interp1d(z_, vp_, dep);
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

void FixVpRho::get_vs_limits(const Eigen::ArrayXd &z, double vs_width,
                             Eigen::ArrayXd &vs_ref, Eigen::ArrayXd &lb,
                             Eigen::ArrayXd &ub) {
  vs_ref = interp_vs(z);
  ArrayXd dep = z2interpdepth(z);
  ArrayXd vp_ref = interp1d(z_, vp_, dep);
  ArrayXd vs_min = vs_ref - vs_width / 2.0;
  ArrayXd vs_max = vs_ref + vs_width / 2.0;
  double tiny = 1.0e-2;
  for (int i = 0; i < vs_min.rows(); ++i) {
    vs_min(i) = std::max(0.01, vs_min(i));
    if (vs_max(i) > vp_ref(i) - tiny) {
      vs_max(i) = vp_ref(i) - tiny;
    }
  }
  lb = vs_min;
  ub = vs_max;
}

Eigen::ArrayXXd Brocher05::generate(const Eigen::ArrayXd &z,
                                    const Eigen::ArrayXd &vs) {
  const int nl = vs.size();

  ArrayXXd model(nl, 5);
  for (int i = 0; i < nl; ++i) {
    model(i, 0) = i + 1;
    model(i, 1) = z(i);
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
void NearSurface::set_param(const std::vector<double> &param) {
  vp2vs_ = param[0];
}

Eigen::ArrayXXd Gardner::generate(const Eigen::ArrayXd &z,
                                  const Eigen::ArrayXd &vs) {
  const int nl = vs.size();

  ArrayXXd model(nl, 5);
  for (int i = 0; i < nl; ++i) {
    model(i, 0) = i + 1;
    model(i, 1) = z(i);
    model(i, 3) = vs(i);
    model(i, 4) = vp2vs_ * vs[i];
    model(i, 2) = alpha_ * pow(model(i, 4) * 1000.0, beta_);
  }
  return model;
}

Dict Gardner::derivative(const Eigen::ArrayXd &vs) {
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

Eigen::ArrayXXd NearSurface::generate(const Eigen::ArrayXd &z,
                                      const Eigen::ArrayXd &vs) {
  const int nl = vs.size();

  ArrayXXd model(nl, 5);
  for (int i = 0; i < nl; ++i) {
    model(i, 0) = i + 1;
    model(i, 1) = z(i);
    model(i, 3) = vs(i);
    model(i, 4) = vp2vs_ * vs[i];
    model(i, 2) = a_ * pow(vs[i], 2) + b_ * vs[i] + c_;
  }
  return model;
}

Dict NearSurface::derivative(const Eigen::ArrayXd &vs) {
  ArrayXd dvs = ArrayXd::Ones(vs.rows());
  ArrayXd dvp = vp2vs_ * dvs;
  ArrayXd drho = (2 * a_ * vs + b_) * dvs;

  Dict res;
  res["rho"] = drho;
  res["vp"] = dvp;
  res["vs"] = dvs;
  return res;
}

Eigen::ArrayXd generate_random_depth(int N, double zmax, double min_gap_raw) {
  double min_gap = min_gap_raw / zmax;
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
    depth.push_back(sequence[i]);
  }

  ArrayXd ret = Map<ArrayXd, Unaligned>(depth.data(), depth.size());
  ret *= zmax;

  return ret;
}

Eigen::ArrayXd generate_depth_by_layer_ratio(double lmin, double lmax,
                                             double ratio, double zmax) {
  double depmax = lmax / 2.0;
  if (zmax > depmax)
    depmax = zmax;
  std::vector<double> depth{0.0};
  depth.push_back(ratio * lmin / 3.0);
  double d = depth.back() + ratio * depth.back();
  depth.push_back(d);

  while (depth.back() < depmax) {
    d = depth.back() + ratio * (depth.back() - depth[depth.size() - 2]);
    depth.push_back(d);
  }

  ArrayXd depth_a = Map<ArrayXd, Unaligned>(depth.data(), depth.size());
  return depth_a;
}