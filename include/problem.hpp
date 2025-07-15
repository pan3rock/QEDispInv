#ifndef PROBLEM_H_
#define PROBLEM_H_

#include <Eigen/Dense>
#include <deque>
#include <memory>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

class Data {
public:
  Data(const Eigen::Ref<const Eigen::ArrayXXd> disp);
  void add_sigma(const std::vector<double> &sigma);
  std::set<int> mode;
  std::unordered_map<int, std::vector<double>> freq, c, sigma;
  int count;

  double lmin, lmax;
  double cmin, cmax;
};

Data resample(const Data &data);

class Vs2Model;

namespace lbfgspp {

class DispersionCurves {
public:
  DispersionCurves(const Eigen::VectorXd &z_model,
                   const Eigen::VectorXd &vs_ref,
                   std::shared_ptr<Vs2Model> vs2model,
                   const std::vector<double> &weight, double lamb_vs, bool sh,
                   int rtype);
  // regularization:
  // rtype == 1: d/dz
  // rtype == 2: An2020
  ~DispersionCurves();
  void load_data(Data &data);
  int update_fitness(const Eigen::VectorXd &vs);
  Eigen::ArrayXXd x2model(const Eigen::VectorXd &vs);
  double operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad);
  int num_forward() { return num_forward_; }

  Eigen::ArrayXXd forward(const Eigen::Ref<const Eigen::ArrayXXd> model);
  std::pair<double, double> get_last_fval(const Eigen::VectorXd &vs);

private:
  double f_reg(const Eigen::VectorXd &vs);
  Eigen::VectorXd g_reg(const Eigen::VectorXd &vs);

  void update_matM_An2020();
  void update_matM_tr1();

  const Eigen::VectorXd z_model_, vs_ref_;
  const int nx_;
  std::shared_ptr<Vs2Model> vs2model_;
  Eigen::ArrayXd weight_;
  Eigen::MatrixXd matM_;
  const double lamb_vs_;
  const bool sh_;
  const int nl_;
  const bool water_;
  Data *data_;
  std::vector<double> f_obs_, c_obs_;
  std::vector<int> m_obs_;

  std::deque<double> f_computed_;
  std::deque<Eigen::VectorXd> x_computed_, g_computed_;
  int num_forward_ = 0;

  const size_t max_deque_ = 10;
};
} // namespace lbfgspp

Eigen::ArrayXXd compute_hist2d(const std::vector<Eigen::ArrayXd> &z_inv,
                               const std::vector<Eigen::ArrayXd> &vs_inv,
                               const std::vector<double> &fitness, double vsmin,
                               double vsmax, double zmax, int num_hist,
                               Eigen::ArrayXd &z_samples,
                               Eigen::ArrayXd &vs_samples);

void compute_statistics(const Eigen::ArrayXd &z, const Eigen::ArrayXd &vs,
                        const Eigen::Ref<const Eigen::ArrayXXd> hist,
                        Eigen::ArrayXd &vs_mean, Eigen::ArrayXd &vs_median,
                        Eigen::ArrayXd &vs_mode, Eigen::ArrayXd &vs_cred10,
                        Eigen::ArrayXd &vs_cred90);

std::vector<size_t> detect_outliers(const std::vector<double> &fitness,
                                    double multiplier = 1.5);
#endif