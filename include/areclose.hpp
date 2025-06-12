// AreClose
// ========
//
// C++ equivalent of Python's math.isclose function
// [well that should be called areclose too! :-)]
// Has one small improvement from Boost.
//
// Copyright © 2021, Shriramana Sharma, samjnaa-at-gmail-dot-com
//
// Use, modification and distribution are permitted subject to the
// "BSD-2-Clause"-type license stated in the accompanying file LICENSE.txt

// NOTES:
// Code adapted from
// https://github.com/python/cpython/blob/master/Modules/mathmodule.c#L2932 The
// logic of the formulation is well argued out at
// https://www.python.org/dev/peps/pep-0485/ The default values are also from
// there. See also: https://docs.python.org/dev/library/math.html#math.isclose

// USAGE:
//     AreClose areClose1, areClose2{2e-12, 4 * DBL_EPSILON};
//     if (areClose1(a, b) || areClose2(a, b))
//         ....

#include <cfloat>
#include <cmath>

struct AreClose {
public:
  AreClose() : abs_tol{0}, rel_tol{1e-09} {}

  AreClose(double absTol, double relTol) : abs_tol{absTol}, rel_tol{relTol} {
    if (invalidTolerance(abs_tol) || invalidTolerance(rel_tol))
      abs_tol = rel_tol = NAN;
    // this makes operator() always return false
  }

  bool isValid() const {
    return abs_tol == abs_tol && rel_tol == rel_tol;
  } // quick NAN test

  bool operator()(double a, double b) const {
    // The following is needed to catch infinities of the same sign since their
    // difference is NAN. And rarely, finite inputs may also be equal.
    if (a == b)
      return true;

    // Since an infinity would have an infinite relative tolerance,
    // any finite number would be considered relatively close to an infinity.
    // Further we need to catch infinities of opposite signs whose
    // difference is infinite and would compare equal to their relative
    // tolerance. The following is needed for both those cases:
    if (std::isinf(a) || std::isinf(b))
      return false;

    // The below is effectively the same as:
    // abs(a - b) <= max(abs_tol, rel_tol * max(abs(a), abs(b)))
    double diff = std::fabs(a - b);
    return diff <= abs_tol || diff <= std::fabs(rel_tol * a) ||
           diff <= std::fabs(rel_tol * b);
  }

  double absoluteTolerance() const { return abs_tol; }

  bool setAbsoluteTolerance(double absTol) // returns success status
  {
    if (invalidTolerance(absTol))
      return false;
    abs_tol = absTol;
    return true;
  }

  double relativeTolerance() const { return rel_tol; }

  bool setRelativeTolerance(double relTol) // returns success status
  {
    if (invalidTolerance(relTol))
      return false;
    rel_tol = relTol;
    return true;
  }

  static bool invalidTolerance(double tol) {
    return !(tol == 0.0 || (tol >= 4 * DBL_EPSILON && std::isfinite(tol)));
    // A tolerance may be zero to switch off that particular type of tolerance
    // checking but if it is non-zero then it should be finite and at least 4 ×
    // machine epsilon See:
    // https://www.boost.org/doc/libs/1_75_0/libs/math/doc/html/math_toolkit/roots_noderiv/root_termination.html
  }

private:
  double abs_tol, rel_tol;
};