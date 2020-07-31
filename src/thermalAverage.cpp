#include "thermalAverage.h"

thermalAverage::thermalAverage(Settings *sett) {
  settings = sett;
  //Set structure for the numeric integration
  wb1 = new double[settings->NBetaIntegral];
  xb1 = new double[settings->NBetaIntegral];
  cgqf(settings->NBetaIntegral, 0, 1, xb1, wb1);

  wt1 = new double[settings->NAngularIntegral];
  xt1 = new double[settings->NAngularIntegral];
  cgqf(settings->NAngularIntegral, 0, M_PI, xt1, wt1);
}

thermalAverage::~thermalAverage() {
  delete[] wb1;
  delete[] xb1;
  delete[] wt1;
  delete[] xt1;
}

double thermalAverage::thermalAverageCrossection2(double T, int iProcess) {
  double res = 0.0;
  double s = 0.0;
  for (size_t i = 0; i < (size_t)settings->NBetaIntegral; i++) {
    s = 4.0 * pow(settings->mchi, 2) / (1 - pow(xb1[i], 2));
    res += wb1[i] * xb1[i] / pow(1 - pow(xb1[i], 2), 2) * sqrt(s) *
           (s - pow(2.0 * settings->mchi, 2)) *
           std::cyl_bessel_k(1, sqrt(s) / T) * totalCrossection(s, iProcess);
  }

  return res * 2.0 * M_PI * M_PI * T * 8.0 * pow(settings->mchi, 2) /
         pow(4.0 * M_PI * pow(settings->mchi, 2) * T *
                 std::cyl_bessel_k(2, settings->mchi / T),
             2);
}

//Calculate the weighting funktion Wij
double thermalAverage::Wij(double s) {
  double matEleVal;
  double res = 0.0;
  for (size_t i = 0; i < (size_t)settings->NAngularIntegral; i++) {
    res += wt1[i] * sin(xt1[i]) / 2.0 *
           MatrixElements::TotalElement(s, xt1[i], settings);
  }

  return res / (8.0 * M_PI * s);
}
//Approximate form of the bessel funktion appearing in the diff. eq.
double thermalAverage::BesselExpand(double x, double xbar) {
  double res = 6.726783879201919 + 0.5479626202923281 * (-15. + x) +
               0.0009591879029919937 * pow(-15. + x, 2) -
               0.00005655647189783126 * pow(-15. + x, 3);
  return res / sqrt(x) * exp(-xbar);
}

//Calculate the thermalaveraged Crossection 
double thermalAverage::thermalAverageCrossection(double T) {

  double res = 0.0;
  double s = 0.0;
  double buff = 0.0;
  double x = settings->mchi / T;
  for (size_t i = 0; i < (size_t)settings->NBetaIntegral; i++) {
    s = 4.0 * pow(settings->mchi, 2) / (1 - pow(xb1[i], 2));
    buff = pow(settings->mchi, 2) * wb1[i] * xb1[i] /
           pow(1 - xb1[i] * xb1[i], 2) *
           sqrt(lambda(s, pow(settings->mchi, 2), pow(settings->mchi, 2)) / s) *
           Wij(s);
    buff *= sqrt(2.0 / (M_PI * sqrt(s) * T)) /
            (2.0 * T * pow(settings->gchi * pow(settings->mchi, 2), 2)) *
            settings->mchi *
            exp(2.0 * x * (1.0 - 1.0 / sqrt(1 - xb1[i] * xb1[i])));
    buff *= (1.0 + 3.0 * T / 8.0 / sqrt(s) - 15.0 * T * T / 128.0 / s +
             105.0 * T * T * T / 1024.0 / pow(s, 1.5));
    buff /= pow(1.0 + 15.0 / (8.0 * x) + 105.0 / (128.0 * x * x) -
                    315.0 / (1024.0 * x * x * x),
                2);
    res += buff;
  }
  return res;
}

double thermalAverage::lambda(double a, double b, double c) {
  return a * a + b * b + c * c - 2.0 * (a * b + b * c + a * c);
}

//Calculate the equilibrium number density at temperature T
double thermalAverage::Neq(double T) {
  return settings->gchi * pow(settings->mchi * T / (2 * M_PI), 1.5) *
         exp(-settings->mchi / (T));
}

//Calculate the total crossection at energy s
double thermalAverage::totalCrossection(double s, int iProcess) {
  double res = 0.0;
  for (size_t i = 0; i < (size_t)settings->NAngularIntegral; i++) {
    res += wt1[i] * MatrixElements::TotalElement(s, xt1[i], settings) *
           sin(xt1[i]);
  }
  return res / (s * s * 32.0 * M_PI);
}
