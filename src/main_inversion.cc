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
#include "rjmcmc.hpp"
#include "utils.hpp"

#include <CLI11.hpp>
#include <Eigen/Dense>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <stdexcept>
#include <string>
#include <toml.hpp>
#include <vector>

using namespace Eigen;

int main(int argc, char const *argv[]) {
  CLI::App app{"Calculating dispersion curves."};

  std::string file_config = "config.toml";
  app.add_option("-c,--config", file_config, "toml-type configure file");
  std::string file_data;
  app.add_option("-d,--data", file_data, "filename of dispersion curves");
  bool sh = false;
  app.add_flag("--sh", sh, "whether are Love waves");
  std::string file_out = "inv.h5";
  app.add_option("-o,--out", file_out, "filename of output");

  CLI11_PARSE(app, argc, argv);

  const auto config = toml::parse(file_config);
  auto conf_inv = toml::find(config, "inversion");
  const auto zmax = toml::find<double>(conf_inv, "zmax");
  const auto vsmin = toml::find<double>(conf_inv, "vsmin");
  const auto vsmax = toml::find<double>(conf_inv, "vsmax");
  const auto thkmin = toml::find<double>(conf_inv, "thkmin");
  const auto weight = toml::find<std::vector<double>>(conf_inv, "weight");
  const auto burnin_samples = toml::find<int>(conf_inv, "burnin_samples");
  const auto total_samples = toml::find<int>(conf_inv, "total_samples");
  const auto linear_model = toml::find<bool>(conf_inv, "linear_model");

  auto vs2model = toml::find<std::string>(conf_inv, "vs2model");
  std::transform(vs2model.begin(), vs2model.end(), vs2model.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  std::shared_ptr<Model> pmodel;
  if (vs2model == "nearsurface") {
    pmodel = std::make_shared<NearSurface>();
  } else if (vs2model == "fixvprho") {
    auto file_model = toml::find<std::string>(conf_inv, "model_reference");
    ArrayXXd model_ref = loadtxt(file_model);
    double zmax_ref = model_ref(model_ref.rows() - 1, 1);
    if (zmax_ref < zmax) {
      std::string msg = fmt::format("maximum depth ({:.4f}) of reference model "
                                    "should be greater than zmax ({:.4f})",
                                    zmax_ref, zmax);
      throw std::runtime_error(msg);
    }
    pmodel = std::make_shared<FixVpRho>(model_ref);
  } else if (vs2model == "brocher05") {
    pmodel = std::make_shared<Brocher05>();
  } else {
    std::string msg = fmt::format(
        "invalid vs2model: {:s}, not in (NearSurface, FixVpRho, Brocher05)",
        vs2model);
    throw std::runtime_error(msg);
  }

  Parameter par;
  par.z_max = zmax;
  par.thk_min = thkmin;
  par.vs_min = vsmin;
  par.vs_max = vsmax;
  par.burnin_samples = burnin_samples;
  par.total_samples = total_samples;
  par.linear_model = linear_model;

  Timer clock;
  clock.tick();
  Data data(file_data);
  RjMcMC rm(par, pmodel, data, weight, sh);
  rm.run(file_out);
  clock.tock();
  fmt::print("elasped time: {:.3f} s\n", clock.duration().count() / 1.0e3);

  return 0;
}
