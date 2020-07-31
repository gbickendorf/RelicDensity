#include "matrixelements.h"

/*
++++++++++++++++++++++++++++++++++++++++++++++++
This File contains Matrix elements automatically
constructed in Mathematica.
++++++++++++++++++++++++++++++++++++++++++++++++
*/



//Define Kaehlerfunction
double lambda(double a, double b, double c) {
  return a * a + b * b + c * c - 2.0 * (a * b + b * c + a * c);
}

double MatrixElements::MatrixElementToElectron(double s, double t,
                                               Settings *sett) {
  if (s < pow(2.0 * me, 2))
    return 0.0;
  double mmed = sett->mmed;
  double wmed = sett->wmed;
  double mchi = sett->mchi;
  double cchi = sett->cchi;
  double csm = sett->csm;
  return (-4 * pow(cchi, 2) * pow(csm, 2) * (4 * pow(mchi, 2) - s) *
          (-4 * pow(me, 2) + s)) /
         pow(pow(mmed, 2) - s, 2);
}

double MatrixElements::MatrixElementToMuon(double s, double t, Settings *sett) {
  if (s < pow(2.0 * mmu, 2))
    return 0.0;
  double mmed = sett->mmed;
  double wmed = sett->wmed;
  wmed = 0.0;
  double mchi = sett->mchi;
  double cchi = sett->cchi;
  double csm = sett->csm;
  return (-4 * pow(cchi, 2) * pow(csm, 2) * (4 * pow(mchi, 2) - s) *
          (-4 * pow(mmu, 2) + s)) /
         pow(pow(mmed, 2) - s, 2);
}

double MatrixElements::MatrixElementToMed(double s, double t, Settings *sett) {
  if (s <= pow(2.0 * sett->mmed, 2))
    return 0.0;
  double mmed = sett->mmed;
  double wmed = sett->wmed;
  double mchi = sett->mchi;
  double cchi = sett->cchi;
  double csm = sett->csm;
  return -(
      (pow(cchi, 4) * (4 * pow(mchi, 2) - s) *
       (192 * pow(mchi, 2) * pow(mmed, 4) -
        224 * pow(mchi, 2) * pow(mmed, 2) * s + 16 * pow(mmed, 4) * s +
        76 * pow(mchi, 2) * pow(s, 2) - 8 * pow(mmed, 2) * pow(s, 2) +
        pow(s, 3) -
        16 * pow(mchi, 2) *
            (16 * pow(mmed, 4) - 16 * pow(mmed, 2) * s + 3 * pow(s, 2)) *
            cos(2 * t) +
        (4 * pow(mchi, 2) - s) * pow(-4 * pow(mmed, 2) + s, 2) * cos(4 * t))) /
      pow(pow(-2 * pow(mmed, 2) + s, 2) +
              (4 * pow(mchi, 2) - s) * (-4 * pow(mmed, 2) + s) * pow(cos(t), 2),
          2));
}

std::complex<double> F12(double t) {
  using namespace std::complex_literals;
  if (t >= 1)
    return -2.0 * t * (1.0 + (1.0 - t) * pow(asin(1 / sqrt(t)), 2.0));
  else {
    std::complex<double> base, cpow;
    base = -1.0i * M_PI + log((1.0 + sqrt(1.0 - t)) / (1.0 - sqrt(1.0 - t)));
    cpow = 2.0;
    return -2.0 * t * (1.0 - (1.0 - t) / 4.0 * pow(base, cpow));
  }
}

double InducedCoupling(double s, Settings *sett) {
  double csm = sett->csm;
  using namespace std::complex_literals;
  return abs(csm / mmu * F12(4.0 * mmu * mmu / s) +
             csm / mmu * F12(4.0 * me * me / s) +
             csm / mmu * F12(4.0 * mt * mt / s)) /
         (137.0 * 2.0 * M_PI);
}

double MatrixElements::MatrixElementToPhoton(double s, double t,
                                             Settings *sett) {
  double g = InducedCoupling(s, sett);
  double mchi = sett->mchi;
  double mmed = sett->mmed;
  return (-16 * pow(g, 2) * (4 * pow(mchi, 2) - s) * pow(s, 2)) /
         pow(pow(mmed, 2) - s, 2);
}

double MatrixElements::Width(Settings *sett) {
  double csm = sett->csm;
  double mmed = sett->mmed;
  return (pow(csm, 2) * pow(-4 * pow(me, 2) + pow(mmed, 2), 1.5)) /
         (12. * pow(mmed, 2) * M_PI);
}

double MatrixElements::TotalElement(double s, double t, Settings *sett) {
  double res = 0.0;
  double buff;

  buff = MatrixElementToMed(s, t, sett);
  if (buff != 0.0)
    res += sqrt(lambda(s, pow(sett->mmed, 2), pow(sett->mmed, 2))) * buff;

  buff = MatrixElementToElectron(s, t, sett);
  if (buff != 0.0)
    res += me * me / mmu / mmu * sqrt(lambda(s, pow(me, 2), pow(me, 2))) * buff;

  buff = MatrixElementToMuon(s, t, sett);
  if (buff != 0.0)
    res += sqrt(lambda(s, pow(mmu, 2), pow(mmu, 2))) * buff;

  buff = MatrixElementToPhoton(s, t, sett);
  if (buff > 0.0)
    res += sqrt(lambda(s, 0, 0)) * buff;
  return res;
}

double MatrixElements::Kinfactor1(Settings *sett, double s) {
  return sqrt(lambda(s, pow(sett->mmed, 2), pow(sett->mmed, 2)));
}

double MatrixElements::Kinfactor2(Settings *sett, double s) {
  return sqrt(lambda(s, pow(me, 2), pow(me, 2)));
}
