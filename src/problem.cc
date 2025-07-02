#include "problem.hpp"
#include "disp.hpp"
#include "model.hpp"
#include "swegn96.hpp"
#include "utils.hpp"

#include <Eigen/Dense>
#include <fmt/format.h>
#include <random>
#include <utility>

using namespace Eigen;

Data::Data(const Eigen::Ref<const Eigen::ArrayXXd> disp) : count(disp.rows()) {
  for (int i = 0; i < disp.rows(); ++i) {
    int m = static_cast<int>(disp(i, 2));
    freq[m].push_back(disp(i, 0));
    c[m].push_back(disp(i, 1));
    mode.insert(m);
  }

  if (disp.cols() == 4) {
    for (int i = 0; i < disp.rows(); ++i) {
      int m = static_cast<int>(disp(i, 2));
      sigma[m].push_back(disp(i, 3));
    }
  }
}

void Data::add_sigma(const std::vector<double> &sigma_in) {
  for (size_t m = 0; m < sigma_in.size(); ++m) {
    int num = freq[m].size();
    for (int j = 0; j < num; ++j) {
      sigma[m].push_back(sigma_in[m]);
    }
  }
}

Data resample(const Data &data) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  Data data_resampled = data;
  auto modes = data.mode;
  for (auto it = modes.begin(); it != modes.end(); ++it) {
    int mode = *it;
    if (data.sigma.count(mode)) {
      size_t nf = data_resampled.freq[mode].size();
      auto &c = data_resampled.c[mode];
      auto &s = data_resampled.sigma[mode];
      for (size_t i = 0; i < nf; ++i) {
        c[i] = c[i] + dist(gen) * s[i];
      }
    }
  }

  return data_resampled;
}

namespace lbfgspp {
DispersionCurves::DispersionCurves(const VectorXd &z_model,
                                   const VectorXd &vs_ref,
                                   std::shared_ptr<Vs2Model> vs2model,
                                   const std::vector<double> &weight,
                                   double lamb_vs, bool sh, int rtype)
    : z_model_(z_model), vs_ref_(vs_ref), nx_(z_model.rows()),
      vs2model_(vs2model), weight_(ArrayXd::Zero(20)),
      matM_(MatrixXd::Zero(nx_, nx_)), lamb_vs_(lamb_vs), sh_(sh),
      nl_(z_model.rows()), water_(vs_ref(0) == 0.0) {
  for (size_t i = 0; i < weight.size(); ++i) {
    weight_(i) = weight[i];
  }
  if (rtype == 1) {
    update_matM_tr1();
  } else if (rtype == 2) {
    update_matM_An2020();
  }
}

DispersionCurves::~DispersionCurves() = default;

void DispersionCurves::update_matM_An2020() {
  MatrixXd matL = MatrixXd::Zero(nx_ - 1, nx_);
  ArrayXd diff = (vs_ref_.head(nl_ - 1) - vs_ref_.tail(nl_ - 1)).array().abs();
  double a = diff.maxCoeff() * 0.1;
  ArrayXd w = a / (a + diff);

  for (int i = 0; i < nl_ - 1; ++i) {
    matL(i, i) = w(i);
    matL(i, i + 1) = -w(i);
  }
  matL.topRows(nl_ - 1) *= sqrt(lamb_vs_);
  matM_ = matL.transpose() * matL;
}

void DispersionCurves::update_matM_tr1() {
  MatrixXd matL = MatrixXd::Zero(nx_ - 1, nx_);

  for (int i = 0; i < nl_ - 1; ++i) {
    matL(i, i) = 1.0;
    matL(i, i + 1) = -1.0;
  }
  matL.topRows(nl_ - 1) *= sqrt(lamb_vs_);
  matM_ = matL.transpose() * matL;
}

void DispersionCurves::load_data(Data &data) {
  data_ = &data;
  f_obs_.clear();
  c_obs_.clear();
  m_obs_.clear();

  for (auto it = data.mode.begin(); it != data.mode.end(); ++it) {
    int mode = *it;
    if (weight_(mode) == 0.) {
      continue;
    }
    std::vector<double> freq = data.freq[mode];
    std::vector<double> c = data.c[mode];
    std::vector<int> modes(freq.size(), mode);
    f_obs_.insert(f_obs_.end(), freq.begin(), freq.end());
    c_obs_.insert(c_obs_.end(), c.begin(), c.end());
    m_obs_.insert(m_obs_.end(), modes.begin(), modes.end());
  }

  num_forward_ = 0;
}

int DispersionCurves::update_fitness(const VectorXd &x) {
  ++num_forward_;
  auto equal = [&](const VectorXd &x2) { return (x - x2).norm() < 1.0e-8; };
  auto x_exist = std::find_if(x_computed_.begin(), x_computed_.end(), equal);
  if (x_exist != x_computed_.end()) {
    int index = x_exist - x_computed_.begin();
    return index;
  }

  ArrayXXd model = x2model(x);
  auto der = vs2model_->derivative(x);

  Dispersion disp(model, sh_);
  SwEgn96 se(model, sh_);

  ArrayXd c_syn(f_obs_.size());
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (size_t i = 0; i < f_obs_.size(); ++i) {
    c_syn(i) = disp.search_mode(f_obs_[i], m_obs_[i]);
  }

  ArrayXi counts = ArrayXi::Zero(weight_.size());
  for (size_t i = 0; i < f_obs_.size(); ++i) {
    if (!std::isnan(c_syn(i)))
      counts(m_obs_[i]) += 1;
  }

  double f = 0.0;
  VectorXd grad = VectorXd::Zero(x.rows());

  for (size_t i = 0; i < f_obs_.size(); ++i) {
    if (!std::isnan(c_syn(i))) {
      int mode = m_obs_[i];
      f += weight_(mode) / counts(mode) * std::pow(c_syn(i) - c_obs_[i], 2);
      auto ker = se.kernel(f_obs_[i], c_syn(i));
      grad += (2.0 * weight_(mode) / counts(mode) * (c_syn[i] - c_obs_[i]) *
               (ker["rho"] * der["rho"] + ker["vs"] * der["vs"] +
                ker["vp"] * der["vp"]))
                  .matrix();
    }
  }

  f += f_reg(x);
  grad += g_reg(x);

  if (x_computed_.size() == max_deque_) {
    x_computed_.pop_front();
    f_computed_.pop_front();
    g_computed_.pop_front();
  }
  x_computed_.push_back(x);
  f_computed_.push_back(f);
  g_computed_.push_back(grad);
  return x_computed_.size() - 1;
}

std::pair<double, double>
DispersionCurves::get_last_fval(const Eigen::VectorXd &x) {
  update_fitness(x);
  double freg = f_reg(x);
  double fres = f_computed_.back() - freg;
  return std::make_pair(fres, freg / lamb_vs_);
}

Eigen::ArrayXXd
DispersionCurves::forward(const Eigen::Ref<const Eigen::ArrayXXd> model) {
  Dispersion disp(model, sh_);
  std::vector<ArrayXd> list_disp;
  for (auto it = data_->mode.begin(); it != data_->mode.end(); ++it) {
    int mode = *it;
    std::vector<double> freq = data_->freq[mode];
    std::vector<double> c_data = data_->c[mode];
    for (size_t i = 0; i < freq.size(); ++i) {
      double c_find = disp.search_mode(freq[i], mode);
      if (!std::isnan(c_find)) {
        ArrayXd d(3);
        d << freq[i], c_find, static_cast<int>(mode);
        list_disp.push_back(d);
      }
    }
  }

  ArrayXXd d3(list_disp.size(), 3);
  for (size_t i = 0; i < list_disp.size(); ++i) {
    d3.row(i) = list_disp[i];
  }
  return d3;
}

Eigen::ArrayXXd DispersionCurves::x2model(const Eigen::VectorXd &vs) {
  return vs2model_->generate(z_model_, vs);
}

double DispersionCurves::operator()(const Eigen::VectorXd &x,
                                    Eigen::VectorXd &grad) {
  int index = update_fitness(x);
  grad = g_computed_[index];
  return f_computed_[index];
}

double DispersionCurves::f_reg(const Eigen::VectorXd &vs) {
  auto mat = (vs - vs_ref_).transpose() * matM_ * (vs - vs_ref_);
  double f = 1.0 / nx_ * mat(0, 0);
  return f;
}

Eigen::VectorXd DispersionCurves::g_reg(const Eigen::VectorXd &vs) {
  VectorXd grad = 2.0 / nx_ * matM_ * (vs - vs_ref_);
  return grad;
}

} // namespace lbfgspp

Eigen::ArrayXXd compute_hist2d(const std::vector<Eigen::ArrayXd> &z_inv,
                               const std::vector<Eigen::ArrayXd> &vs_inv,
                               const std::vector<double> &fitness, double vsmin,
                               double vsmax, double zmax, int num_hist,
                               Eigen::ArrayXd &z_samples,
                               Eigen::ArrayXd &vs_samples) {
  double dz = zmax / (num_hist - 1);
  double dvs = (vsmax - vsmin) / (num_hist - 1);

  ArrayXd f_val = Map<const ArrayXd, Unaligned>(fitness.data(), fitness.size());
  f_val = 1.0 / f_val;
  f_val /= f_val.maxCoeff();

  z_samples = ArrayXd::LinSpaced(num_hist, 0, zmax);
  vs_samples = ArrayXd::LinSpaced(num_hist, vsmin, vsmax);

  ArrayXXd hist2d = ArrayXXd::Zero(num_hist, num_hist);
  int nl = z_inv[0].rows();
  for (size_t n = 0; n < z_inv.size(); ++n) {
    ArrayXd z = z_inv[n];
    ArrayXd vs = vs_inv[n];
    int i_z = 0;
    for (int i = 0; i < nl; ++i) {
      int i1_v = static_cast<int>(floor((vs(i) - vsmin) / dvs));

      double zub;
      if (i == nl - 1) {
        zub = zmax;
      } else {
        zub = z(i + 1);
      }

      while (i_z * dz <= zub) {
        hist2d(i_z, i1_v) += f_val[n];
        ++i_z;
      }
      // if (i != nl - 1) {
      //   int i2_v = i1_v;
      //   if (vs(i) < vs(i + 1)) {
      //     while (vsmin + i2_v * dvs < vs(i + 1) && i2_v < num_hist) {
      //       hist2d(i_z, i2_v) += 1;
      //       ++i2_v;
      //     }
      //   } else {
      //     while (vsmin + i2_v * dvs > vs(i + 1) && i2_v >= 0) {
      //       hist2d(i_z, i2_v) += 1;
      //       --i2_v;
      //     }
      //   }
      // }

      // i_z += 1;
    }
  }
  return hist2d;
}

void compute_statistics(const Eigen::ArrayXd &z, const Eigen::ArrayXd &vs,
                        const Eigen::Ref<const Eigen::ArrayXXd> hist,
                        Eigen::ArrayXd &vs_mean, Eigen::ArrayXd &vs_median,
                        Eigen::ArrayXd &vs_mode, Eigen::ArrayXd &vs_cred10,
                        Eigen::ArrayXd &vs_cred90) {
  const int n = vs.size();
  for (int i_z = 0; i_z < z.rows(); ++i_z) {
    const double total_weight = hist.row(i_z).sum();
    if (total_weight <= 0) {
      throw std::invalid_argument("Total weight must be positive");
    }

    ArrayXd hist1d = hist.row(i_z);
    vs_mean(i_z) = (vs * hist1d).sum() / total_weight;

    int max_idx = 0;
    hist1d.maxCoeff(&max_idx);
    vs_mode(i_z) = vs(max_idx);

    Eigen::VectorXd cum_weights(n);
    cum_weights(0) = hist1d(0);
    for (int i = 1; i < n; ++i) {
      cum_weights(i) = cum_weights(i - 1) + hist1d(i);
    }

    auto computeQuantile = [&](double quantile) -> double {
      const double target_weight = total_weight * quantile;

      int idx = 0;
      while (idx < n && cum_weights(idx) < target_weight) {
        idx++;
      }

      if (idx == 0) {
        return vs(0);
      } else if (idx == n) {
        return vs(n - 1);
      }

      const double weight_before = cum_weights(idx - 1);
      const double weight_here = cum_weights(idx);
      const double x_before = vs(idx - 1);
      const double x_here = vs(idx);

      if (weight_here <= weight_before) {
        return x_before;
      }

      const double t =
          (target_weight - weight_before) / (weight_here - weight_before);
      return x_before + t * (x_here - x_before);
    };

    vs_median(i_z) = computeQuantile(0.5);
    vs_cred10(i_z) = computeQuantile(0.1);
    vs_cred90(i_z) = computeQuantile(0.9);
  }
}

std::vector<size_t> detect_outliers(const std::vector<double> &fitness,
                                    double multiplier) {

  std::vector<double> qs = percentiles(fitness, {25.0, 75.0});
  double iqr = qs[1] - qs[0];
  double lower_bound = qs[0] - multiplier * iqr;
  double upper_bound = qs[1] + multiplier * iqr;
  std::vector<size_t> find;
  for (int i = fitness.size() - 1; i >= 0; --i) {
    double f = fitness[i];
    if (f < lower_bound or f > upper_bound) {
      find.push_back(i);
    }
  }
  return find;
}