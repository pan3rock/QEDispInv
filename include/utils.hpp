#ifndef UTILS_H_
#define UTILS_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

Eigen::ArrayXXd loadtxt(const std::string &file_model,
                        const std::string &delim = " \t");

std::vector<size_t> argsort(const Eigen::ArrayXd &x);

std::vector<double> percentiles(const std::vector<double> &data,
                                const std::vector<double> &ps);

template <typename T>
void remove_by_indices(std::vector<T> &data, std::vector<size_t> indices) {
  for (size_t index : indices) {
    data.erase(data.begin() + index);
  }
}

double min_varray(const std::vector<Eigen::ArrayXd> &va);
double max_varray(const std::vector<Eigen::ArrayXd> &va);
#endif