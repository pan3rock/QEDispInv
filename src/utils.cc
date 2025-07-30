/*
 QEDispInv: Surface Wave Dispersion Curve Computation and Inversion Toolkit

 GNU General Public License, Version 3, 29 June 2007

 Copyright (c) 2025 Lei Pan

 Xiaofei Chen Research Group,
 Department of Earth and Space Sciences,
 Southern University of Science and Technology, China.
 */

#include "utils.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <vector>

const double INF = std::numeric_limits<double>::infinity();

using namespace Eigen;

ArrayXXd loadtxt(const std::string &filename, const std::string &delim) {
  std::filesystem::path path(filename);
  std::string msg = "file " + filename + " doesn't exist.";
  if (!exists(path)) {
    throw std::runtime_error(msg);
  }

  std::ifstream infile(filename);
  std::string line;
  std::string strnum;

  auto data = std::vector<std::vector<double>>();

  data.clear();

  while (getline(infile, line)) {
    if (line.empty()) {
      continue;
    }

    data.push_back(std::vector<double>());

    for (std::string::const_iterator i = line.begin(); i != line.end(); ++i) {
      if (delim.find(*i) == std::string::npos) {
        strnum += *i;
        if (i + 1 != line.end())
          continue;
      }

      if (strnum.empty())
        continue;

      double number;

      std::istringstream(strnum) >> number;
      data.back().push_back(number);

      strnum.clear();
    }
  }

  ArrayXXd eArray(data.size(), data[0].size());
  for (size_t i = 0; i < data.size(); ++i)
    eArray.row(i) = Eigen::ArrayXd::Map(&data[i][0], data[0].size());
  return eArray;
}

std::vector<size_t> argsort(const Eigen::ArrayXd &v) {
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

  return idx;
}

std::vector<double> percentiles(const std::vector<double> &data,
                                const std::vector<double> &ps) {
  if (data.size() == 0) {
    throw std::invalid_argument("Data cannot be empty");
  }

  std::vector<double> sorted_data = data;
  std::sort(sorted_data.begin(), sorted_data.end());

  std::vector<double> results;
  for (double p : ps) {
    if (p < 0 || p > 100) {
      throw std::invalid_argument("Percentiles must be between 0 and 100");
    }

    if (p == 100.0) {
      results.push_back(sorted_data.back());
      continue;
    }

    const double n = sorted_data.size();
    const double pos = p * (n - 1) / 100.0;
    const size_t k = static_cast<size_t>(pos);
    const double f = pos - k;

    if (k + 1 >= n) {
      results.push_back(sorted_data[k]);
    } else {
      results.push_back(sorted_data[k] +
                        f * (sorted_data[k + 1] - sorted_data[k]));
    }
  }

  return results;
}

double min_varray(const std::vector<Eigen::ArrayXd> &va) {
  double vmin = INF;
  for (size_t i = 0; i < va.size(); ++i) {
    double min = va[i].minCoeff();
    vmin = std::min(min, vmin);
  }
  return vmin;
}

double max_varray(const std::vector<Eigen::ArrayXd> &va) {
  double vmax = -INF;
  for (size_t i = 0; i < va.size(); ++i) {
    double max = va[i].maxCoeff();
    vmax = std::max(max, vmax);
  }
  return vmax;
}