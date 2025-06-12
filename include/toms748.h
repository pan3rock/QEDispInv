// toms748
// =======
//
// TOMS 748 algorithm of Alefeld, Potra and Shi to find a root of a function in
// a given interval using inverse cubic and “Newton-quadratic” interpolations.
//
// This is a C++ implementation with a C-compatible interface and an optional
// C++ helper. Copyright © 2021, Shriramana Sharma, samjnaa-at-gmail-dot-com
//
// Use, modification and distribution are permitted subject to the
// "BSD-2-Clause"-type license stated in the accompanying file LICENSE.txt

#ifdef __cplusplus
extern "C" {
#endif

// type declarations

enum Toms748ErrorCode {
  // errors in function input, may be OR-ed together
  TOMS748_NO_ERROR = 0,
  TOMS748_INVALID_INTERVAL_START = 1,
  TOMS748_INVALID_INTERVAL_END = 2,
  TOMS748_INVALID_INTERVAL = 4,
  TOMS748_INVALID_ABSOLUTE_TOLERANCE = 8,
  TOMS748_INVALID_RELATIVE_TOLERANCE = 16,
  TOMS748_INVALID_MAXIMUM_ITERATIONS = 32,
  TOMS748_INVALID_INTERPOLATIONS_PER_ITERATION = 64,
  // errors in function execution, mutually exclusive
  TOMS748_INVALID_FUNCTION_VALUE = 128,
  TOMS748_INTERVAL_DOES_NOT_BRACKET_A_ROOT = 129,
  TOMS748_MAXIMUM_ITERATIONS_REACHED = 130
};

typedef struct {
  int iterations, functionCalls;
  double bracketStart, bracketEnd;
  int errorCode;
} Toms748ResultStatus;

typedef double (*Toms748InputFunction)(double, void *);

// function declarations

double toms748Custom(Toms748InputFunction function, void *otherInput,
                     double intervalStart, double intervalEnd,
                     Toms748ResultStatus *resultStatus,
                     double absoluteTolerance, double relativeTolerance,
                     int maximumIterations, int interpolationsPerIteration);

double toms748(Toms748InputFunction function, void *otherInput,
               double intervalStart, double intervalEnd,
               Toms748ResultStatus *resultStatus);

#ifdef __cplusplus
} // extern "C"

template <typename Functor>
double _toms748customFunction(double var, void *fixed) {
  return static_cast<Functor *>(fixed)->operator()(var);
}

template <typename Functor>
double toms748(Functor functor, double intervalStart, double intervalEnd,
               Toms748ResultStatus *resultStatus = nullptr) {
  return toms748(&_toms748customFunction<Functor>,
                 static_cast<void *>(&functor), intervalStart, intervalEnd,
                 resultStatus);
}

template <typename Functor>
double toms748(Functor functor, double intervalStart, double intervalEnd,
               Toms748ResultStatus *resultStatus, double absoluteTolerance,
               double relativeTolerance, int maximumIterations,
               int interpolationsPerIteration) {
  return toms748Custom(
      &_toms748customFunction<Functor>, static_cast<void *>(&functor),
      intervalStart, intervalEnd, resultStatus, absoluteTolerance,
      relativeTolerance, maximumIterations, interpolationsPerIteration);
}
#endif