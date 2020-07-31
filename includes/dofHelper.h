#pragma once

#include "legendre_rule.h"
#include "settings.h"
#include <math.h>
#include <stdio.h>
#include <vector>

class dofHelper {
public:
  static void Init();
  static double DOF(double T, Settings *sett);
  static double DOFs(double T, Settings *sett);
  static double dofSMContribution(double T);
  static double dofsSMContribution(double T);

private:
  // static Settings* settings;
  static int GetIndexClosest(double x, double arr[], int N);
  static double Interpolate(double x, double arrx[], double arry[], int N);
  static double xconts[];
  static double contrFermion[];
  static double contrBoson[];
  static double scontrFermion[];
  static double scontrBoson[];
  static double gEeff[];
  static double TBenchmark[];
  static double gSeff[];
  static double TBenchmarks[];
};
