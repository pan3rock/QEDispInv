/*
 QEDispInv: Surface Wave Dispersion Curve Computation and Inversion Toolkit

 GNU General Public License, Version 3, 29 June 2007

 Copyright (c) 2025 Lei Pan

 Xiaofei Chen Research Group,
 Department of Earth and Space Sciences,
 Southern University of Science and Technology, China.
 */

#include "disp.hpp"
#include "secfunc.hpp"
#include "utils.hpp"

#include <CLI11.hpp>
#include <Eigen/Dense>
#include <fmt/format.h>
#include <highfive/H5Easy.hpp>
#include <string>
#include <toml.hpp>

using namespace Eigen;

int main(int argc, char const *argv[]) {
  CLI::App app{"Scanning secular function given a model and frequency."};

  double freq;
  app.add_option("freq", freq, "frequency")->required();
  std::string file_config = "config.toml";
  app.add_option("-c,--config", file_config, "toml-type configure file");
  bool sh = false;
  app.add_flag("--sh", sh, "whether to compute Love waves");
  std::string file_model = "";
  app.add_option("--model", file_model, "filename of model");
  std::string file_out = "secfunc.h5";
  app.add_option("-o,--out", file_out, "filename of output");

  CLI11_PARSE(app, argc, argv);

  const auto config = toml::parse(file_config);
  auto conf_secfunc = toml::find(config, "secfunc");
  if (file_model == "") {
    file_model = toml::find<std::string>(conf_secfunc, "file_model");
  }
  const int nc = toml::find<int>(conf_secfunc, "nc");

  ArrayXXd model = loadtxt(file_model);
  SecularFunction sf(model, sh);

  double cmin = model.col(3).minCoeff() * 0.8;
  double cmax = model.col(3).maxCoeff();
  ArrayXd c = ArrayXd::LinSpaced(nc, cmin, cmax);

  ArrayXd sfunc(nc);
  for (int i_c = 0; i_c < nc; ++i_c) {
    sfunc(i_c) = sf.evaluate(freq, c(i_c));
  }

  Dispersion disp(model, sh);
  auto samples = disp.get_samples(freq);
  ArrayXd N(samples.size());
  for (size_t i = 0; i < samples.size(); ++i) {
    N(i) = disp.approx(freq, samples[i]);
  }

  std::vector<double> val(samples.size());
  for (size_t i_c = 0; i_c < samples.size(); ++i_c) {
    val[i_c] = sf.evaluate(freq, samples[i_c]);
  }

  std::vector<double> x_ext, y_ext;
  disp.locate_extremum(freq, samples, val, x_ext, y_ext);

  const int maxmode = 1000;
  auto roots = disp.search(freq, maxmode);

  H5Easy::File fout(file_out, H5Easy::File::Overwrite);
  H5Easy::dump(fout, "f", freq);
  H5Easy::dump(fout, "c", c);
  H5Easy::dump(fout, "sfunc", sfunc);
  H5Easy::dump(fout, "samples", samples);
  H5Easy::dump(fout, "samples_ext", x_ext);
  H5Easy::dump(fout, "roots", roots);
  H5Easy::dump(fout, "N", N);

  return 0;
}
