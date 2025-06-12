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
#include "utils.hpp"

#include <CLI11.hpp>
#include <Eigen/Dense>
#include <LBFGSB.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
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
  const auto lamb_vs = toml::find<double>(conf_inv, "lamb_vs");
  const auto rtype = toml::find<int>(conf_inv, "reg_type");
  const auto weight = toml::find<std::vector<double>>(conf_inv, "weight");

  ArrayXXd model_ref = loadtxt(file_mref);

  std::transform(vs2model.begin(), vs2model.end(), vs2model.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  std::shared_ptr<Vs2Model> pmodel;
  if (vs2model == "nearsurface") {
    pmodel = std::make_shared<NearSurface>();
  } else if (vs2model == "fixvprho") {
    pmodel = std::make_shared<FixVpRho>(model_ref);
  } else if (vs2model == "brocher05") {
    pmodel = std::make_shared<Brocher05>();
  } else {
    std::string msg = fmt::format(
        "invalid vs2model: {:s}, not in (NearSurface, FixVpRho, Brocher05)",
        vs2model);
    throw std::runtime_error(msg);
  }

  // vs boundary
  ArrayXd vs_ref = model_ref.col(3);
  ArrayXd vs_min = vs_ref - vs_width / 2.0;
  ArrayXd vs_max = vs_ref + vs_width / 2.0;
  double tiny = 1.0e-2;
  ArrayXd vp_ref = model_ref.col(4);
  for (int i = 0; i < vs_min.rows(); ++i) {
    if (vs_min(i) < tiny) {
      vs_min(i) = tiny;
    }
    if (vs2model == "fixvprho" && vs_max(i) > vp_ref(i) - tiny) {
      vs_max(i) = vp_ref(i) - tiny;
    }
  }
  ArrayXd lb = vs_min;
  ArrayXd ub = vs_max;

  ArrayXXd data_input = loadtxt(file_data);
  Data data(data_input);

  const int nl = model_ref.rows();
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0.0, 1.0);
  auto gen_rand_minit = [&]() -> ArrayXd {
    ArrayXd rand = ArrayXd::NullaryExpr(nl, [&]() { return dist(gen); });
    ArrayXd x = lb.array() + (ub - lb).array() * rand;
    return x;
  };

  LBFGSBParam<double> param;
  param.max_iterations = 100;
  LBFGSBSolver<double> solver(param);

  const auto num_init = toml::find<int>(conf_inv, "num_init");
  const auto num_noise = toml::find<int>(conf_inv, "num_noise");

  for (int i = 0; i < num_noise * num_init; ++i) {
  }
  for (int i_d = 0; i_d < num_noise; ++i_d) {
    Data data_resampled = resample(data);
    for (int i_m = 0; i_m < num_init; ++i_m) {
      ArrayXd z_model = model_ref.col(1);
      lbfgspp::DispersionCurves prob(z_model, vs_ref, pmodel, weight, lamb_vs,
                                     sh, rtype);
      VectorXd x = gen_rand_minit();
      double feval = 0.0;
      int num_iter;
      try {
        num_iter = solver.minimize(prob, x, feval, lb, ub);
      } catch (const std::exception &exc) {
        ;
        // std::cerr << exc.what() << std::endl;
      }
    }
  }

  Timer clock;
  clock.tick();
  clock.tock();
  fmt::print("elasped time: {:.3f} s\n", clock.duration().count() / 1.0e3);

  return 0;
}
