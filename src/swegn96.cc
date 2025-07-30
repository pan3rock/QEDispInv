/*
 QEDispInv: Surface Wave Dispersion Curve Computation and Inversion Toolkit

 GNU General Public License, Version 3, 29 June 2007

 Copyright (c) 2025 Lei Pan

 Xiaofei Chen Research Group,
 Department of Earth and Space Sciences,
 Southern University of Science and Technology, China.
 */
#include "swegn96.hpp"

#include <Eigen/Dense>
#include <map>
#include <string>

using namespace Eigen;
using Dict = std::map<std::string, Eigen::ArrayXd>;

SwEgn96::SwEgn96(const Eigen::Ref<const Eigen::ArrayXXd> &model, bool sh,
                 bool sphere)
    : nl_(model.rows()), tn_(ArrayXf::Zero(nl_)),
      rho_(model.col(2).cast<float>()), vs_(model.col(3).cast<float>()),
      vp_(model.col(4).cast<float>()), sh_(sh), sphere_(sphere) {
  for (int i = 0; i < nl_ - 1; ++i) {
    tn_(i) = model(i + 1, 1) - model(i, 1);
  }
}

Dict SwEgn96::kernel(double freq, double c) {
  double period = 1.0 / freq;
  ArrayXd dc_dvp = ArrayXd::Zero(nl_);
  ArrayXd dc_dvs = ArrayXd::Zero(nl_);
  ArrayXd dc_drho = ArrayXd::Zero(nl_);
  ArrayXd dc_dh = ArrayXd::Zero(nl_);

  int iflsph = sphere_ ? 1 : 0;
  double cg;
  if (sh_) {
    ArrayXd ut(nl_), tt(nl_);
    slegn96_(tn_.data(), vs_.data(), rho_.data(), nl_, &period, &c, &cg,
             ut.data(), tt.data(), dc_dvs.data(), dc_dh.data(), dc_drho.data(),
             iflsph);
  } else {
    ArrayXd ur(nl_), uz(nl_), tr(nl_), tz(nl_);
    sregn96_(tn_.data(), vp_.data(), vs_.data(), rho_.data(), nl_, &period, &c,
             &cg, ur.data(), uz.data(), tr.data(), tz.data(), dc_dvp.data(),
             dc_dvs.data(), dc_dh.data(), dc_drho.data(), iflsph);
  }

  Dict res;
  res["vp"] = dc_dvp;
  res["vs"] = dc_dvs;
  res["rho"] = dc_drho;
  res["h"] = dc_dh;
  return res;
}
