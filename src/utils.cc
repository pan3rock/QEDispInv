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

#include "utils.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <vector>

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