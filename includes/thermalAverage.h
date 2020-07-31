#pragma once

#include "dofHelper.h"
#include "legendre_rule.h"
#include "matrixelements.h"
#include "settings.h"
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <stdio.h>
#include <vector>

class thermalAverage {
public:
  thermalAverage(Settings *sett);
  ~thermalAverage();
  double thermalAverageCrossection(double T);
  double thermalAverageCrossection2(double T, int iProcess);

  double Wij(double s);
  double totalCrossection(double s, int iProcess);
  Settings *settings;
  double Neq(double T);

private:
  double *wb1;
  double *xb1;
  double *wt1;
  double *xt1;
  double BesselExpand(double x, double xbar);
  double lambda(double a, double b, double c);
};
