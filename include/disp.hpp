#ifndef DISP_H_
#define DISP_H_

#include <Eigen/Dense>
#include <memory>
#include <vector>

class SecularFunction;

class Dispersion {
public:
  Dispersion(const Eigen::Ref<const Eigen::ArrayXXd> model, bool sh);
  ~Dispersion();
  std::vector<double> search(double f, int num_mode);
  double search_mode(double f, int mode);
  std::vector<double> get_samples(double f);
  std::vector<double> predict_extremum(double f,
                                       const std::vector<double> &samples);
  std::vector<std::pair<double, double>> find_coarse_intv(double f,
                                                          int num_mode);
  double approx(double f, double c);

private:
  double evaluate_rayleigh_velocity();
  void locate_extremum(double f, const std::vector<double> &x,
                       const std::vector<double> &y, std::vector<double> &x_ext,
                       std::vector<double> &y_ext);

  // parameters for initial samplings
  const double ednn_ = 0.50;
  const int nfine_ = 2;
  const double ctol_ = 1.0e-5;

  // parameters for intervals bracketing roots
  const int nbias_ = 2;
  const int niter_ext_ = 2;

  const int nl_;
  Eigen::ArrayXd thk_, vs_, vp_;
  const bool sh_;
  const int itop_;
  std::unique_ptr<SecularFunction> sf_;
  double vs0_, vp0_, vs_min_, vs_max_, vs_hf_, rayv_;
};

#endif