#pragma once
class Settings;
#include "physicalConstants.h"
#include "settings.h"
#include <cmath>
#include <complex>

class MatrixElements {
public:
  static double TotalElement(double s, double t, Settings *sett);
  static double Width(Settings *sett);

  static double Kinfactor1(Settings *sett, double s);
  static double Kinfactor2(Settings *sett, double s);
  static double Kinfactor3(Settings *sett, double s);
  static double MatrixElementToElectron(double s, double t, Settings *sett);
  static double MatrixElementToMuon(double s, double t, Settings *sett);
  static double MatrixElementToMed(double s, double t, Settings *sett);
  static double MatrixElementToPhoton(double s, double t, Settings *sett);
};
