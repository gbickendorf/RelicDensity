#include "settings.h"

//Set settings to default values
Settings::Settings() {
  x0 = 16.0;
  x1 = 30.0;
  xinf = 1000;
  xmin = x0;

  dxfine = 1E-4;
  dxrough = 1E-2;

  mchi = 0.1;
  gchi = 4.0;
  IsFermion_chi = 1;

  mmed = 0.01;
  wmed = 10.0;
  gmed = 3;
  IsFermion_med = 0;
  cchi = 1.0;
  csm = 1.0;

  IsAntiSame = 0.0;
  NAngularIntegral = 10;
  NBetaIntegral = 100;
  eps_ode = 0.01;
  eps_xf = 0.01;
  FastModePrecision = 1E-8;
}

void Settings::PrintSettings() {
  printf("mchi  %E\ngchi  %E\nmmed  %E\nwmed  %E\ncchi  %E\ncsm %E\n\n x0 "
         "%E\nx1  %E",
         mchi, gchi, mmed, wmed, cchi, csm, x0, x1);
}
