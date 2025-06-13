#ifndef UTILS_H_
#define UTILS_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

Eigen::ArrayXXd loadtxt(const std::string &file_model,
                        const std::string &delim = " \t");

std::vector<size_t> argsort(const Eigen::ArrayXd &x);

#endif