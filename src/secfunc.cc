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

#include "secfunc.hpp"

#include <Eigen/Dense>
#include <complex>

using namespace Eigen;
using namespace std::complex_literals;
using Cmplx = std::complex<double>;

SecularFunction::SecularFunction(const Eigen::Ref<const Eigen::ArrayXXd> model,
                                 bool sh)
    : nl_(model.rows()), thk_(nl_ - 1), dns_(model.col(2)), vs_(model.col(3)),
      vp_(model.col(4)), sh_(sh), is_water_(vs_(0) == 0.0) {
  ArrayXd z = model.col(1);
  thk_ = z.tail(nl_ - 1) - z.head(nl_ - 1);

  if (is_water_) {
    vs_(0) = 1.0e-8;
    iwater_ = 1;
  } else {
    iwater_ = 0;
  }
}

double SecularFunction::evaluate(double f, double c) {
  if (sh_) {
    return evaluate_sh(f, c);
  } else {
    return evaluate_psv(f, c);
  }
}

double SecularFunction::evaluate_psv(double f, double c) {
  double omega = 2.0 * M_PI * f;
  double wvno = omega / c;
  ArrayXd e = ArrayXd::Zero(5);
  ArrayXd ee = ArrayXd::Zero(5);

  double wvno2 = wvno * wvno;
  double xka = omega / vp_[nl_ - 1];
  double xkb = omega / vs_[nl_ - 1];
  double wvnop = wvno + xka;
  double wvnom = abs(wvno - xka);
  double ra = sqrt(wvnop * wvnom);
  wvnop = wvno + xkb;
  wvnom = abs(wvno - xkb);
  double rb = sqrt(wvnop * wvnom);
  double t = vs_[nl_ - 1] / omega;

  // E matrix for the bottom half-space
  double gammk = 2.0 * t * t;
  double gam = gammk * wvno2;
  double gamm1 = gam - 1.0;
  double rho1 = dns_[nl_ - 1];
  e[0] = rho1 * rho1 * (gamm1 * gamm1 - gam * gammk * ra * rb);
  e[1] = -rho1 * ra;
  e[2] = rho1 * (gamm1 - gammk * ra * rb);
  e[3] = rho1 * rb;
  e[4] = wvno2 - ra * rb;

  // Matrix multiplication from bottom layer upward
  for (int m = nl_ - 2; m >= iwater_; --m) {
    xka = omega / vp_[m];
    xkb = omega / vs_[m];
    t = vs_[m] / omega;
    gammk = 2.0 * t * t;
    gam = gammk * wvno2;
    wvnop = wvno + xka;
    wvnom = abs(wvno - xka);
    ra = sqrt(wvnop * wvnom);
    wvnop = wvno + xkb;
    wvnom = abs(wvno - xkb);
    rb = sqrt(wvnop * wvnom);

    double dpth = thk_[m];
    rho1 = dns_[m];
    double p = ra * dpth;
    double q = rb * dpth;

    // Evaluate cosP, cosQ...
    double w, cosp, a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz;
    var(p, q, ra, rb, wvno, xka, xkb, dpth, w, cosp, a0, cpcq, cpy, cpz, cqw,
        cqx, xy, xz, wy, wz);

    // Evaluate Dunkin's matrix
    ArrayXXd ca = dnka(wvno2, gam, gammk, rho1, a0, cpcq, cpy, cpz, cqw, cqx,
                       xy, xz, wy, wz);

    ee = (ca.matrix().transpose() * e.matrix()).array();
    e = ee;
  }

  double dlt;
  if (iwater_ == 1) {
    xka = omega / vp_[0];
    wvnop = wvno + xka;
    wvnom = abs(wvno - xka);
    ra = sqrt(wvnop * wvnom);
    double dpth = thk_[0];
    rho1 = dns_[0];
    double p = ra * dpth;
    double zul = 1.0e-5;
    double w, cosp, a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz;
    var(p, zul, ra, zul, wvno, xka, xkb, dpth, w, cosp, a0, cpcq, cpy, cpz, cqw,
        cqx, xy, xz, wy, wz);
    dlt = cosp * e[0] - rho1 * w * e[1];
  } else {
    dlt = e[0];
  }
  return dlt;
}

void SecularFunction::var(double p, double q, double ra, double rb, double wvno,
                          double xka, double xkb, double dpth, double &w,
                          double &cosp, double &a0, double &cpcq, double &cpy,
                          double &cpz, double &cqw, double &cqx, double &xy,
                          double &xz, double &wy, double &wz) {
  // Examine P - wave eigenfunctions
  // Checking whether c > vp, c = vp or c < vp
  double pex = 0.0;
  double x = 0.0;
  if (wvno < xka) {
    double sinp = sin(p);
    w = sinp / ra;
    x = -ra * sinp;
    cosp = cos(p);
  } else if (wvno == xka) {
    cosp = 1.0;
    w = dpth;
    x = 0.0;
  } else {
    pex = p;
    double fac = 0.0;
    if (p < 16.0) {
      fac = exp(-2.0 * p);
    }
    cosp = (1.0 + fac) * 0.5;
    double sinp = (1.0 - fac) * 0.5;
    w = sinp / ra;
    x = ra * sinp;
  }

  // Examine S-wave eigenfunctions
  // Checking whether c > vs, c = vs or c < vs
  double sex = 0.0;
  double sinq = 0.0, y = 0.0, z = 0.0, cosq = 0.0;
  if (wvno < xkb) {
    sinq = sin(q);
    y = sinq / rb;
    z = -rb * sinq;
    cosq = cos(q);
  } else if (wvno == xkb) {
    cosq = 1.0;
    y = dpth;
    z = 0.0;
  } else {
    sex = q;
    double fac = 0.0;
    if (q < 16.0) {
      fac = exp(-2.0 * q);
    }
    cosq = (1.0 + fac) * 0.5;
    sinq = (1.0 - fac) * 0.5;
    y = sinq / rb;
    z = rb * sinq;
  }

  // Form eigenfunction products for use with compound matrices
  double exa = pex + sex;
  if (exa < 60.0) {
    a0 = exp(-exa);
  } else {
    a0 = 0.0;
  }
  cpcq = cosp * cosq;
  cpy = cosp * y;
  cpz = cosp * z;
  cqw = cosq * w;
  cqx = cosq * x;
  xy = x * y;
  xz = x * z;
  wy = w * y;
  wz = w * z;
}

ArrayXXd SecularFunction::dnka(double wvno2, double gam, double gammk,
                               double rho, double a0, double cpcq, double cpy,
                               double cpz, double cqw, double cqx, double xy,
                               double xz, double wy, double wz) {

  ArrayXXd ca = ArrayXXd::Zero(5, 5);

  double gamm1 = gam - 1.0;
  double twgm1 = gam + gamm1;
  double gmgmk = gam * gammk;
  double gmgm1 = gam * gamm1;
  double gm1sq = gamm1 * gamm1;

  double rho2 = rho * rho;
  double a0pq = a0 - cpcq;
  double t = -2.0 * wvno2;

  ca(0, 0) = cpcq - 2.0 * gmgm1 * a0pq - gmgmk * xz - wvno2 * gm1sq * wy;
  ca(0, 1) = (wvno2 * cpy - cqx) / rho;
  ca(0, 2) = -(twgm1 * a0pq + gammk * xz + wvno2 * gamm1 * wy) / rho;
  ca(0, 3) = (cpz - wvno2 * cqw) / rho;
  ca(0, 4) = -(2.0 * wvno2 * a0pq + xz + wvno2 * wvno2 * wy) / rho2;

  ca(1, 0) = (gmgmk * cpz - gm1sq * cqw) * rho;
  ca(1, 1) = cpcq;
  ca(1, 2) = gammk * cpz - gamm1 * cqw;
  ca(1, 3) = -wz;
  ca(1, 4) = ca(0, 3);

  ca(3, 0) = (gm1sq * cpy - gmgmk * cqx) * rho;
  ca(3, 1) = -xy;
  ca(3, 2) = gamm1 * cpy - gammk * cqx;
  ca(3, 3) = ca(1, 1);
  ca(3, 4) = ca(0, 1);

  ca(4, 0) =
      -(2.0 * gmgmk * gm1sq * a0pq + gmgmk * gmgmk * xz + gm1sq * gm1sq * wy) *
      rho2;
  ca(4, 1) = ca(3, 0);
  ca(4, 2) = -(gammk * gamm1 * twgm1 * a0pq + gam * gammk * gammk * xz +
               gamm1 * gm1sq * wy) *
             rho;
  ca(4, 3) = ca(1, 0);
  ca(4, 4) = ca(0, 0);

  ca(2, 0) = t * ca(4, 2);
  ca(2, 1) = t * ca(3, 2);
  ca(2, 2) = a0 + 2.0 * (cpcq - ca(0, 0));
  ca(2, 3) = t * ca(1, 2);
  ca(2, 4) = t * ca(0, 2);

  return ca;
}

double SecularFunction::evaluate_sh(double f, double c) {
  double omega = 2.0 * M_PI * f;
  double wvno = omega / c;
  double beta1 = vs_(nl_ - 1);
  double rho1 = dns_(nl_ - 1);
  double xkb = omega / beta1;
  double wvnop = wvno + xkb;
  double wvnom = abs(wvno - xkb);
  double rb = sqrt(wvnop * wvnom);
  double e1 = rho1 * rb;
  double e2 = 1.0 / (beta1 * beta1);

  for (int m = nl_ - 2; m >= iwater_; --m) {
    beta1 = vs_(m);
    rho1 = dns_(m);
    double xmu = rho1 * beta1 * beta1;
    xkb = omega / beta1;
    wvnop = wvno + xkb;
    wvnom = abs(wvno - xkb);
    rb = sqrt(wvnop * wvnom);
    double q = thk_(m) * rb;

    double sinq, y, z, cosq;
    if (wvno < xkb) {
      sinq = sin(q);
      y = sinq / rb;
      z = -rb * sinq;
      cosq = cos(q);
    } else if (wvno == xkb) {
      cosq = 1.0;
      y = thk_(m);
      z = 0.0;
    } else {
      double fac = 0.0;
      if (q < 16.0) {
        fac = exp(-2.0 * q);
      }
      cosq = (1.0 + fac) * 0.5;
      sinq = (1.0 - fac) * 0.5;
      y = sinq / rb;
      z = rb * sinq;
    }

    double e10 = e1 * cosq + e2 * xmu * z;
    double e20 = e1 * y / xmu + e2 * cosq;
    double norm = std::max(abs(e10), abs(e20));
    if (norm > 1.0e-40) {
      e1 = e10 / norm;
      e2 = e20 / norm;
    }
  }
  return e1;
}