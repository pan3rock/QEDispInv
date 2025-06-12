/*
 Copyright (c) 2024 Lei Pan

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 of the Software, and to permit persons to whom the Software is furnished to do
 so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#include "disp.hpp"
#include "secfunc.hpp"
#include "toms748.h"

#include <Eigen/Dense>
#include <fmt/format.h>
#include <numeric>

using namespace Eigen;

namespace {
double calculate_newton_step(const double root, const double vp,
                             const double vs) {
  double p = pow(root, -1);
  double p2 = pow(root, -2);
  double vs_2 = pow(vs, -2);
  double vp_2 = pow(vp, -2);
  double xi = sqrt(p2 - vs_2);
  double eta = sqrt(p2 - vp_2);
  double func = pow(vs_2 - 2.0 * p2, 2) - 4.0 * xi * eta * p2;
  double deri = p2 * (8.0 * p * (vs_2 - 2 * p2) + 8.0 * p * xi * eta +
                      4.0 * p2 * p * (xi / eta + eta / xi));
  return func / deri;
}

double find_extremum(double x1, double x2, double x3, double f1, double f2,
                     double f3) {
  double B1 = (pow(x2, 2) - pow(x3, 2)) * f1;
  double B2 = (pow(x3, 2) - pow(x1, 2)) * f2;
  double B3 = (pow(x1, 2) - pow(x2, 2)) * f3;
  double C1 = (x2 - x3) * f1;
  double C2 = (x3 - x1) * f2;
  double C3 = (x1 - x2) * f3;
  double D = (x1 - x2) * (x2 - x3) * (x3 - x1);
  double b = (B1 + B2 + B3) / D;
  double a = -(C1 + C2 + C3) / D;
  return -b / (2.0 * a);
}

std::vector<double> sort_by_index(const std::vector<double> &x,
                                  const std::vector<size_t> &index) {
  std::vector<double> ret;
  for (auto i : index) {
    ret.push_back(x[i]);
  }
  return ret;
}
} // namespace

std::vector<int> find_required_nl(const Eigen::ArrayXd &vs) {
  const int nl = vs.rows();
  std::vector<int> ireq;
  for (int i = 1; i < nl - 1; ++i) {
    if (vs(i) < vs(i - 1) && vs(i) < vs(i + 1)) {
      ireq.push_back(i);
    }
  }
  ireq.push_back(nl);
  return ireq;
}

Dispersion::Dispersion(const Eigen::Ref<const Eigen::ArrayXXd> model, bool sh)
    : nl_(model.rows()),
      thk_(model.col(1).tail(nl_ - 1) - model.col(1).head(nl_ - 1)),
      vs_(model.col(3)), vp_(model.col(4)), sh_(sh),
      itop_(model(0, 3) == 0 ? 1 : 0),
      sf_(std::make_unique<SecularFunction>(model, sh)) {
  vs0_ = vs_(0);
  vp0_ = vp_(0);
  vs_min_ = vs_.minCoeff();
  vs_max_ = vs_.maxCoeff();
  vs_hf_ = vs_(nl_ - 1);
  rayv_ = evaluate_rayleigh_velocity();
}

Dispersion::~Dispersion() = default;

double Dispersion::evaluate_rayleigh_velocity() {
  double root = 0.8 * vs0_;
  int count = 0;
  while (true) {
    double update_step = calculate_newton_step(root, vp0_, vs0_);
    root -= update_step;
    if (count > 10 || std::abs(update_step) < ctol_)
      break;
    ++count;
  }
  return root;
}

double Dispersion::approx(double f, double c) {
  double sum = 0.0;
  double c_2 = pow(c, -2);
  for (int i = itop_; i < nl_ - 1; ++i) {
    if (c > vs_(i)) {
      sum += sqrt(pow(vs_(i), -2) - c_2) * thk_(i);
    }
  }

  if (!sh_) {
    for (int i = 0; i < nl_ - 1; ++i) {
      if (c > vp_(i)) {
        sum += sqrt(pow(vp_(i), -2) - c_2) * thk_(i);
      }
    }
  }
  sum *= 2.0 * f;
  return sum;
}

std::vector<double> Dispersion::get_samples(double f) {
  std::vector<double> pred;
  int nmax = static_cast<int>(std::floor(approx(f, vs_hf_))) + 1;
  double dc = (vs_hf_ - vs_min_) / nmax;
  double c1 = vs_min_;
  double e1 = approx(f, c1);
  while (c1 < vs_hf_) {
    double c2 = c1 + dc;
    double e2 = approx(f, c2);
    while (abs(e2 - e1) > ednn_) {
      c2 = c1 + (c2 - c1) * 0.618;
      e2 = approx(f, c2);
    }
    c1 = c2;
    e1 = e2;
    if (c2 < vs_hf_) {
      pred.push_back(c2);
    }
  }

  pred.push_back(0.8 * vs_min_);
  pred.push_back(vs_hf_ - ctol_);
  std::sort(pred.begin(), pred.end());
  for (int i = 0; i < nfine_; ++i) {
    int num = pred.size();
    for (int j = 0; j < num - 1; ++j) {
      double mid = (pred[j] + pred[j + 1]) / 2.0;
      pred.push_back(mid);
    }
    std::sort(pred.begin(), pred.end());
  }

  std::vector<double> samples;
  samples.insert(samples.end(), pred.begin(), pred.end());
  samples.push_back(vs_min_ - ctol_ * 10);
  samples.push_back(vs_min_ + ctol_ * 10);
  std::sort(samples.begin(), samples.end());

  if (!sh_) {
    samples.push_back(rayv_);
    std::sort(samples.begin(), samples.end());
  }

  return samples;
}

std::vector<std::pair<double, double>>
Dispersion::find_coarse_intv(double f, int num_mode) {
  std::vector<double> samples = get_samples(f);

  std::vector<double> x_sample{samples[0]};
  std::vector<double> y_sample{sf_->evaluate(f, samples[0])};

  int count_root = 0;
  for (size_t i = 1; i < samples.size(); ++i) {
    double x = samples[i];
    double y = sf_->evaluate(f, x);
    x_sample.push_back(x);
    y_sample.push_back(y);
    if (y_sample[i - 1] * y_sample[i] < 0)
      ++count_root;
    if (count_root >= num_mode + nbias_)
      break;
  }

  std::vector<double> x_ext, y_ext;
  if (x_sample.size() >= 3) {
    locate_extremum(f, x_sample, y_sample, x_ext, y_ext);
  }

  std::vector<double> x_coarse(x_sample.begin(), x_sample.end());
  std::vector<double> y_coarse(y_sample.begin(), y_sample.end());
  x_coarse.insert(x_coarse.end(), x_ext.begin(), x_ext.end());
  y_coarse.insert(y_coarse.end(), y_ext.begin(), y_ext.end());
  std::vector<size_t> iasc(x_coarse.size());
  std::iota(iasc.begin(), iasc.end(), 0);
  std::sort(iasc.begin(), iasc.end(),
            [&x_coarse](int i, int j) { return x_coarse[i] < x_coarse[j]; });
  x_coarse = sort_by_index(x_coarse, iasc);
  y_coarse = sort_by_index(y_coarse, iasc);

  std::vector<std::pair<double, double>> find_intv;
  for (size_t i = 0; i < x_coarse.size() - 1; ++i) {
    if (y_coarse[i] * y_coarse[i + 1] < 0) {
      find_intv.push_back({x_coarse[i], x_coarse[i + 1]});
    }
  }

  if (int(find_intv.size()) > num_mode) {
    std::vector<std::pair<double, double>> s(find_intv.data(),
                                             find_intv.data() + num_mode);
    return s;
  } else {
    return find_intv;
  }
}

std::vector<double> Dispersion::search(double f, int num_mode) {
  std::function<double(double)> func = [&](double c) {
    double val = sf_->evaluate(f, c);
    return val;
  };
  auto find_intv = find_coarse_intv(f, num_mode);

  std::vector<double> find;
  for (auto it = find_intv.cbegin(); it != find_intv.end(); ++it) {
    double c1 = it->first;
    double c2 = it->second;
    double root = toms748(func, c1, c2);
    if (!std::isnan(root))
      find.push_back(root);
  }

  if (int(find.size()) > num_mode) {
    std::vector<double> s(find.data(), find.data() + num_mode);
    return s;
  } else {
    return find;
  }
}

double Dispersion::search_mode(double f, int mode) {
  auto cs = search(f, mode + 1);
  if (static_cast<int>(cs.size()) < mode + 1) {
    return std::numeric_limits<double>::quiet_NaN();
  } else {
    return cs[mode];
  }
}

std::vector<double>
Dispersion::predict_extremum(double f, const std::vector<double> &samples) {
  std::vector<double> fx;
  for (size_t i = 0; i < samples.size(); ++i) {
    fx.push_back(sf_->evaluate(f, samples[i]));
  }

  std::vector<double> xe;
  for (size_t i = 0; i < samples.size() - 2; ++i) {
    double x1 = samples[i];
    double x2 = samples[i + 1];
    double x3 = samples[i + 2];
    double f1 = fx[i];
    double f2 = fx[i + 1];
    double f3 = fx[i + 2];
    double x = find_extremum(x1, x2, x3, f1, f2, f3);
    xe.push_back(x);
  }
  return xe;
}

void Dispersion::locate_extremum(double f, const std::vector<double> &x,
                                 const std::vector<double> &y,
                                 std::vector<double> &x_ext,
                                 std::vector<double> &y_ext) {
  auto x_tmp = x;
  auto y_tmp = y;

  for (int it = 0; it < niter_ext_; ++it) {
    x_ext.clear();
    y_ext.clear();
    for (size_t i = 0; i < x_tmp.size() - 2; ++i) {
      double x1 = x_tmp[i];
      double x2 = x_tmp[i + 1];
      double x3 = x_tmp[i + 2];
      double f1 = y_tmp[i];
      double f2 = y_tmp[i + 1];
      double f3 = y_tmp[i + 2];
      double xe_i = find_extremum(x1, x2, x3, f1, f2, f3);
      if (xe_i < vs_hf_ && xe_i > 0.8 * vs_min_) {
        x_ext.push_back(xe_i);
        y_ext.push_back(sf_->evaluate(f, xe_i));
      }
    }
    if (it == niter_ext_ - 1 || x_ext.size() < 3)
      break;

    x_tmp = x_ext;
    y_tmp = y_ext;

    std::vector<size_t> iasc(x_tmp.size());
    std::iota(iasc.begin(), iasc.end(), 0);
    std::sort(iasc.begin(), iasc.end(),
              [&x_tmp](int i, int j) { return x_tmp[i] < x_tmp[j]; });
    x_tmp = sort_by_index(x_tmp, iasc);
    y_tmp = sort_by_index(y_tmp, iasc);
  }
}