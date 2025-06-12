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

  H5Easy::File fout(file_out, H5Easy::File::Overwrite);
  H5Easy::dump(fout, "f", freq);
  H5Easy::dump(fout, "c", c);
  H5Easy::dump(fout, "sfunc", sfunc);
  H5Easy::dump(fout, "samples", samples);

  return 0;
}
