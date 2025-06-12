#ifndef UTILS_H_
#define UTILS_H_

#include <Eigen/Dense>
#include <chrono>
#include <string>
#include <vector>

Eigen::ArrayXXd loadtxt(const std::string &file_model,
                        const std::string &delim = " \t");

std::vector<size_t> argsort(const Eigen::ArrayXd &x);

template <class DT = std::chrono::milliseconds,
          class ClockT = std::chrono::steady_clock>
class Timer {
  using timep_t = typename ClockT::time_point;
  timep_t _start = ClockT::now(), _end = {};

public:
  void tick() {
    _end = timep_t{};
    _start = ClockT::now();
  }

  void tock() { _end = ClockT::now(); }

  template <class T = DT> auto duration() const {
    return std::chrono::duration_cast<T>(_end - _start);
  }
};
#endif