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

#include "model.hpp"
#include "problem.hpp"
#include "tqdm.hpp"
#include "utils.hpp"

#include <CLI11.hpp>
#include <Eigen/Dense>
#include <LBFGSB.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <highfive/H5Easy.hpp>
#include <random>
#include <stdexcept>
#include <string>
#include <toml.hpp>
#include <vector>

using namespace Eigen;
using namespace LBFGSpp;

int main(int argc, char const *argv[]) {
  CLI::App app{"Calculating dispersion curves."};

  std::string file_config = "config.toml";
  app.add_option("-c,--config", file_config, "toml-type configure file");
  std::string file_data;
  app.add_option("-d,--data", file_data, "filename of dispersion curves");
  std::string file_mref;
  app.add_option("--model_ref", file_mref, "filename of reference model");
  bool sh = false;
  app.add_flag("--sh", sh, "whether are Love waves");
  std::string file_out = "inv.h5";
  app.add_option("-o,--out", file_out, "filename of output");

  CLI11_PARSE(app, argc, argv);

  const auto config = toml::parse(file_config);
  auto conf_inv = toml::find(config, "inversion");

  if (file_mref == "")
    file_mref = toml::find<std::string>(conf_inv, "model_ref");
  auto vs2model = toml::find<std::string>(conf_inv, "vs2model");
  const auto vs_width = toml::find<double>(conf_inv, "vs_width");
  const auto lamb_vs = toml::find<double>(conf_inv, "lambda");
  const auto rtype = toml::find<int>(conf_inv, "reg_type");
  const auto weight = toml::find<std::vector<double>>(conf_inv, "weight");

  ArrayXXd model_ref = loadtxt(file_mref);

  ArrayXXd data_input = loadtxt(file_data);
  Data data(data_input);
  if (data_input.cols() < 4) {
    const auto sigma = toml::find<std::vector<double>>(conf_inv, "sigma");
    if (sigma.size() < weight.size()) {
      std::string msg = fmt::format(
          "the size of sigma ({:d}) is less than that of weight ({:d})",
          sigma.size(), weight.size());
      throw std::runtime_error(msg);
    }
    data.add_sigma(sigma);
  }

  std::transform(vs2model.begin(), vs2model.end(), vs2model.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  std::shared_ptr<Vs2Model> pmodel;
  if (vs2model == "nearsurface") {
    pmodel = std::make_shared<NearSurface>(model_ref);
  } else if (vs2model == "fixvprho") {
    pmodel = std::make_shared<FixVpRho>(model_ref);
  } else if (vs2model == "brocher05") {
    pmodel = std::make_shared<Brocher05>(model_ref);
  } else {
    std::string msg = fmt::format(
        "invalid vs2model: {:s}, not in (NearSurface, FixVpRho, Brocher05)",
        vs2model);
    throw std::runtime_error(msg);
  }

  LBFGSBParam<double> param;
  param.max_iterations = 100;
  LBFGSBSolver<double> solver(param);

  const auto num_init = toml::find<int>(conf_inv, "num_init");
  const auto num_noise = toml::find<int>(conf_inv, "num_noise");
  const auto rand_depth = toml::find<bool>(conf_inv, "rand_depth");
  const auto dintv_min = toml::find<double>(conf_inv, "dintv_min");
  const auto nl = toml::find<int>(conf_inv, "nlayer");
  const auto zmax = toml::find<double>(conf_inv, "zmax");

  std::vector<Data> data_noise;
  if (num_noise == 1) {
    data_noise.push_back(data);
  } else {
    for (int i_d = 0; i_d < num_noise; ++i_d) {
      Data data_resampled = resample(data);
      data_noise.push_back(data_resampled);
    }
  }

  std::vector<ArrayXd> z_init;
  if (num_init == 1) {
    z_init.push_back(model_ref.col(1));
  } else {
    for (int i_m = 0; i_m < num_init; ++i_m) {
      if (rand_depth) {
        z_init.push_back(generate_random_depth(nl, zmax, dintv_min));
      } else {
        z_init.push_back(model_ref.col(1));
      }
    }
  }

  std::vector<ArrayXd> vs_ref, vs_lb, vs_ub;
  for (int i_m = 0; i_m < num_init; ++i_m) {
    ArrayXd z = z_init[i_m];
    ArrayXd vsr(nl), lb(nl), ub(nl);
    pmodel->get_vs_limits(z, vs_width, vsr, lb, ub);
    vs_ref.push_back(vsr);
    vs_lb.push_back(lb);
    vs_ub.push_back(ub);
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0.0, 1.0);
  auto gen_rand_vsinit = [&](const ArrayXd &lb, const ArrayXd &ub) -> ArrayXd {
    ArrayXd rand = ArrayXd::NullaryExpr(nl, [&]() { return dist(gen); });
    ArrayXd x = lb.array() + (ub - lb).array() * rand;
    return x;
  };

  std::vector<ArrayXd> vs_init;
  if (num_init == 1) {
    vs_init.push_back(model_ref.col(3));
  } else {
    for (int i_m = 0; i_m < num_init; ++i_m) {
      vs_init.push_back(gen_rand_vsinit(vs_lb[i_m], vs_ub[i_m]));
    }
  }

  // for (size_t i = 0; i < z_init.size(); ++i) {
  //   auto z1 = z_init[i];
  //   auto vs1 = vs_init[i];
  //   fmt::print("{:5d}{:5d}{:5d}\n", i, z1.rows(), vs1.rows());
  //   for (int j = 0; j < z1.rows(); ++j) {
  //     fmt::print("{:5d}{:12.5f}{:12.5f}\n", j, z1(j), vs1(j));
  //   }
  //   fmt::print("\n");
  // }

  Tqdm bar;
  int num_total = num_init * num_noise;
  std::vector<ArrayXd> z_inv, vs_inv;
  std::vector<double> fitness;
  std::vector<int> niter;
  std::vector<ArrayXXd> disp_syn;
  for (int i = 0; i < num_total; ++i) {
    bar.progress(i, num_total);
    int i_d = i / num_init;
    int i_m = i % num_init;
    auto data_resampled = data_noise[i_d];
    ArrayXd z_model = z_init[i_m];
    lbfgspp::DispersionCurves prob(z_model, vs_ref[i_m], pmodel, weight,
                                   lamb_vs, sh, rtype);
    prob.load_data(data_resampled);

    VectorXd x = vs_init[i_m];
    double feval = 0.0;
    try {
      int it = solver.minimize(prob, x, feval, vs_lb[i_m], vs_ub[i_m]);
      fitness.push_back(feval);
      niter.push_back(it);
      z_inv.push_back(z_model);
      vs_inv.push_back(x);
      auto model = pmodel->generate(z_model, x);
      auto disp = prob.forward(model);
      disp_syn.push_back(disp);
    } catch (const std::exception &exc) {
      // std::cout << std::endl;
      // auto model = pmodel->generate(z_model, x);
      // std::cout << model << std::endl;
      // std::cerr << exc.what() << std::endl;
    }
  }
  bar.finish();

  std::vector<size_t> idx_outlier = detect_outliers(fitness);
  remove_by_indices(fitness, idx_outlier);
  remove_by_indices(niter, idx_outlier);
  remove_by_indices(z_inv, idx_outlier);
  remove_by_indices(vs_inv, idx_outlier);
  remove_by_indices(disp_syn, idx_outlier);

  // hist of inversion model
  const int num_hist = 100;
  double vsmin = 1.0e10;
  double vsmax = 0.0;
  for (size_t i = 0; i < vs_lb.size(); ++i) {
    double lbmin = vs_lb[i].minCoeff();
    vsmin = std::min(lbmin, vsmin);
    double ubmax = vs_ub[i].maxCoeff();
    vsmax = std::max(ubmax, vsmax);
  }
  ArrayXd z_samples(num_hist), vs_samples(num_hist);
  ArrayXXd hist = compute_hist2d(z_inv, vs_inv, fitness, vsmin, vsmax, zmax,
                                 num_hist, z_samples, vs_samples);

  ArrayXd vs_mean(num_hist), vs_mode(num_hist), vs_median(num_hist),
      vs_cred10(num_hist), vs_cred90(num_hist);
  compute_statistics(z_samples, vs_samples, hist, vs_mean, vs_median, vs_mode,
                     vs_cred10, vs_cred90);

  H5Easy::File out_h5(file_out, H5Easy::File::Overwrite);
  H5Easy::dump(out_h5, "fitness", fitness);
  H5Easy::dump(out_h5, "niter", niter);
  H5Easy::dump(out_h5, "z_sample", z_samples);
  H5Easy::dump(out_h5, "vs_sample", vs_samples);
  H5Easy::dump(out_h5, "vs_hist2d", hist);
  H5Easy::dump(out_h5, "data", data_input);
  H5Easy::dump(out_h5, "vs_mean", vs_mean);
  H5Easy::dump(out_h5, "vs_median", vs_median);
  H5Easy::dump(out_h5, "vs_mode", vs_mode);
  H5Easy::dump(out_h5, "vs_cred10", vs_cred10);
  H5Easy::dump(out_h5, "vs_cred90", vs_cred90);

  ArrayXXd model_mean = pmodel->generate(z_samples, vs_mean);
  H5Easy::dump(out_h5, "model_mean", model_mean);
  ArrayXd vs_ref_save = pmodel->interp_vs(z_samples);
  H5Easy::dump(out_h5, "vs_ref", vs_ref_save);

  std::vector<int> mode_used;
  for (size_t i = 0; i < weight.size(); ++i) {
    if (weight[i] > 0) {
      mode_used.push_back(i);
    }
  }
  H5Easy::dump(out_h5, "mode_used", mode_used);

  int num_valid = fitness.size();
  H5Easy::dump(out_h5, "num_valid", num_valid);
  for (int i = 0; i < num_valid; ++i) {
    std::string key = fmt::format("disp/{:d}", i);
    H5Easy::dump(out_h5, key, disp_syn[i]);
  }

  return 0;
}
