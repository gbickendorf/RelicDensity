#pragma once

#include "dofHelper.h"
#include "legendre_rule.h"
#include "matrixelements.h"
#include "physicalConstants.h"
#include "settings.h"
#include "thermalAverage.h"
#include <cmath>
#include <stdio.h>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

class relic {
public:
  relic(Settings *sett);
  ~relic();
  void Init();
  double CalculateXf();
  double CalculateRelicDensity(int writeFiles, int verbose);
  double CalculateRelicDensity2(int writeFiles, int verbose);
  double CalculateRelicDensityAutoX0(int writeFiles, int verbose);
  double Yeq(double x);
  double fprime_gsl(double xrun, double f);
  void relicFunc();
  Settings *settings;
  thermalAverage *thermAvg;

private:
  vector<double> x;
  vector<double> w;
  vector<double> y;

  double fprime(int i, double f);
  void numericSolver(double f0, double xi, double xf, int N);
};
