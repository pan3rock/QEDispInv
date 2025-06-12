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
#include "utils.hpp"

#include <CLI11.hpp>
#include <Eigen/Dense>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <string>
#include <toml.hpp>
#include <vector>

using namespace Eigen;

int main(int argc, char const *argv[]) {
  CLI::App app{"Calculating dispersion curves."};

  std::string file_config = "config.toml";
  app.add_option("-c,--config", file_config, "toml-type configure file");
  int mode_max = 0;
  app.add_option("-m,--mode", mode_max, "maximum mode up to");
  std::string file_disp = "";
  app.add_option("--disp", file_disp,
                 "input file containing dispersion curves at target "
                 "frequencies for computation");
  bool sh = false;
  app.add_flag("--sh", sh, "whether to compute Love waves");
  std::string file_model = "";
  app.add_option("--model", file_model, "filename of model");
  std::string file_out = "disp.txt";
  app.add_option("-o,--out", file_out, "filename of output");

  CLI11_PARSE(app, argc, argv);

  const auto config = toml::parse(file_config);
  auto dispersion = toml::find(config, "forward");
  if (file_model == "") {
    file_model = toml::find<std::string>(dispersion, "file_model");
  }

  auto model = loadtxt(file_model);

  Dispersion disp(model, sh);
  std::ofstream out(file_out);

  if (file_disp == "") {
    const auto fmin = toml::find<double>(dispersion, "fmin");
    const auto fmax = toml::find<double>(dispersion, "fmax");
    const auto nf = toml::find<int>(dispersion, "nf");
    ArrayXd freqs = ArrayXd::LinSpaced(nf, fmin, fmax);
    for (int i = 0; i < freqs.size(); ++i) {
      auto c = disp.search(freqs(i), mode_max + 1);
      for (size_t m = 0; m < c.size(); ++m) {
        fmt::print(out, "{:15.5f}{:15.7f}{:15d}\n", freqs(i), c[m], m);
      }
    }
  } else {
    ArrayXXd disp_data = loadtxt(file_disp);
    for (int i = 0; i < disp_data.rows(); ++i) {
      double f = disp_data(i, 0);
      int mode = static_cast<int>(disp_data(i, 2));
      double c = disp.search_mode(f, mode);
      if (!std::isnan(c)) {
        fmt::print(out, "{:15.5f}{:15.7f}{:15d}\n", f, c, mode);
      }
    }
  }

  out.close();

  return 0;
}