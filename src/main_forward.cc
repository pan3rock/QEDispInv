/*
 QEDispInv: Surface Wave Dispersion Curve Computation and Inversion Toolkit

 GNU General Public License, Version 3, 29 June 2007

 Copyright (c) 2025 Lei Pan

 Xiaofei Chen Research Group,
 Department of Earth and Space Sciences,
 Southern University of Science and Technology, China.
 */

#include "disp.hpp"
#include "swegn96.hpp"
#include "utils.hpp"

#include <CLI11.hpp>
#include <Eigen/Dense>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <highfive/H5Easy.hpp>
#include <string>
#include <toml.hpp>
#include <vector>

using namespace Eigen;

int main(int argc, char const *argv[]) {
  CLI::App app{"Calculate dispersion curves given a model."};

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
  bool compute_kernel = false;
  app.add_flag("--compute_kernel", compute_kernel, "whether to compute kernel");
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

  if (compute_kernel) {
    ArrayXXd disp = loadtxt(file_out);
    const int nd = disp.rows();
    const int nl = model.rows();
    ArrayXXd kvp(nl, nd), kvs(nl, nd), krho(nl, nd);
    SwEgn96 se(model, sh);
    for (int i = 0; i < disp.rows(); ++i) {
      double f = disp(i, 0);
      double c = disp(i, 1);
      auto ker = se.kernel(f, c);
      kvp.col(i) = ker["vp"];
      kvs.col(i) = ker["vs"];
      krho.col(i) = ker["rho"];
    }
    ArrayXd z = model.col(1);

    std::string file_ker = "kernel.h5";
    H5Easy::File fout(file_ker, H5Easy::File::Overwrite);
    H5Easy::dump(fout, "disp", disp);
    H5Easy::dump(fout, "kvp", kvp);
    H5Easy::dump(fout, "kvs", kvs);
    H5Easy::dump(fout, "krho", krho);
    H5Easy::dump(fout, "z", z);
  }

  return 0;
}