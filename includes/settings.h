#pragma once

#include "matrixelements.h"
#include "physicalConstants.h"
#include <cmath>
#include <functional>
#include <stdio.h>
#include <vector>

class Settings {
public:
  Settings();

  double x0;
  double x1;
  double xinf;
  double xmin;

  double dxfine;
  double dxrough;

  double mchi;
  double gchi;

  int IsFermion_chi;

  double mmed;
  double wmed;
  double gmed;
  int IsFermion_med;
  double cchi;
  double csm;
  double IsAntiSame;
  double FastModePrecision;

  int NAngularIntegral;
  int NBetaIntegral;

  double eps_ode;
  double eps_xf;
  void PrintSettings();
};
