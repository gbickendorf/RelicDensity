#include "dofHelper.h"

//Effective Degrees Of Freedom by energy
double dofHelper::gEeff[] = {
    106.75, 106.75, 106.75, 106.75, 106.75, 106.75, 106.75, 106.75, 106.75,
    106.75, 106.75, 106.75, 106.47, 105.07, 104.96, 104.38, 103.7,  102.23,
    100.78, 98.6,   95.01,  92.95,  91.64,  89.66,  87.72,  87.23,  87.14,
    86.66,  86.09,  85.52,  84.4,   84.31,  82.58,  81.41,  80.26,  78.52,
    76.24,  74.58,  72.42,  69.78,  67.24,  65.28,  63.87,  61.07,  59.75,
    57.56,  23.58,  21.05,  18.65,  17.3,   16.29,  15.81,  15.24,  14.68,
    13.93,  13.12,  12.17,  11.81,  11.39,  11.14,  10.99,  10.98,  10.92,
    10.93,  10.94,  10.95,  10.97,  10.89,  10.91,  10.92,  10.93,  10.85,
    10.78,  10.8,   10.72,  10.65,  10.58,  10.35,  10.05,  9.76,   9.19,
    8.72,   7.85,   7.11,   6.45,   5.13,   4.44,   3.81,   3.43,   3.36,
    3.36,   3.36,   3.36,   3.36,   3.36,   3.36,   3.36,   3.36,   3.36,
    3.36};
//Temperature at which gEeff by energy was taken
double dofHelper::TBenchmark[] = {
    1000,     807.74,   670.87,   557.2,    462.78,   384.35,   319.23,
    265.14,   220.17,   182.85,   151.88,   126.11,   104.74,   86.96,
    72.22,    59.97,    49.8,     41.34,    34.32,    28.49,    23.64,
    19.62,    16.29,    13.52,    11.22,    9.32,     7.74,     6.43,
    5.34,     4.43,     3.68,     3.06,     2.54,     2.11,     1.75,
    1.45,     1.2,      1,        0.83,     0.69,     0.56,     0.47,
    0.39,     0.33,     0.27,     0.22,     0.18,     0.15,     0.12,
    0.1,      8.57E-02, 7.11E-02, 5.90E-02, 4.89E-02, 4.06E-02, 3.37E-02,
    2.79E-02, 2.31E-02, 1.92E-02, 1.59E-02, 1.32E-02, 1.10E-02, 9.13E-03,
    7.58E-03, 6.30E-03, 5.23E-03, 4.34E-03, 3.61E-03, 3.00E-03, 2.49E-03,
    2.07E-03, 1.72E-03, 1.42E-03, 1.18E-03, 9.83E-04, 8.16E-04, 6.78E-04,
    5.62E-04, 4.67E-04, 3.87E-04, 3.21E-04, 2.66E-04, 2.21E-04, 1.83E-04,
    1.51E-04, 1.25E-04, 1.03E-04, 8.54E-05, 7.07E-05, 5.87E-05, 4.87E-05,
    4.04E-05, 3.36E-05, 2.79E-05, 2.32E-05, 1.92E-05, 1.60E-05, 1.33E-05,
    1.10E-05, 9.85E-06};
//Temperature at which gEeff by entropy s was taken
double dofHelper::TBenchmarks[] = {
    1000,     807.74,   670.87,   557.2,    462.78,   384.35,   319.23,
    265.14,   220.17,   182.85,   151.88,   126.11,   104.74,   86.96,
    72.22,    59.97,    49.8,     41.34,    34.32,    28.49,    23.64,
    19.62,    16.29,    13.52,    11.22,    9.32,     7.74,     6.43,
    5.34,     4.43,     3.68,     3.06,     2.54,     2.11,     1.75,
    1.45,     1.2,      1,        0.83,     0.69,     0.56,     0.47,
    0.39,     0.33,     0.27,     0.22,     0.18,     0.15,     0.12,
    0.1,      8.57E-02, 7.11E-02, 5.90E-02, 4.89E-02, 4.06E-02, 3.37E-02,
    2.79E-02, 2.31E-02, 1.92E-02, 1.59E-02, 1.32E-02, 1.10E-02, 9.13E-03,
    7.58E-03, 6.30E-03, 5.23E-03, 4.34E-03, 3.61E-03, 3.00E-03, 2.49E-03,
    2.07E-03, 1.72E-03, 1.42E-03, 1.18E-03, 9.83E-04, 8.16E-04, 6.78E-04,
    5.62E-04, 4.67E-04, 3.87E-04, 3.21E-04, 2.66E-04, 2.21E-04, 1.83E-04,
    1.51E-04, 1.25E-04, 1.03E-04, 6.84E-05, 5.68E-05, 4.71E-05, 4.44E-05,
    3.91E-05, 3.25E-05, 2.70E-05, 2.24E-05, 1.86E-05, 1.55E-05, 1.28E-05,
    1.07E-05, 9.90E-06};
//Effective Degrees Of Freedom by entopie
double dofHelper::gSeff[] = {
    106.75, 106.75, 106.75, 106.75, 106.75, 106.75, 106.75, 106.75, 106.75,
    106.75, 106.75, 106.75, 106.75, 104.96, 105.07, 104.38, 103.7,  102.23,
    100.78, 98.6,   95.01,  92.95,  91.64,  89.66,  87.72,  87.14,  87.23,
    86.66,  86.09,  85.52,  84.31,  84.4,   82.58,  81.41,  80.26,  78.52,
    76.24,  74.58,  72.42,  69.78,  67.24,  65.28,  63.87,  61.07,  59.75,
    57.56,  23.58,  21.05,  18.65,  17.3,   16.29,  15.81,  15.24,  14.68,
    13.93,  13.12,  12.17,  11.81,  11.39,  11.14,  10.98,  10.99,  10.92,
    10.93,  10.94,  10.95,  10.97,  10.89,  10.91,  10.92,  10.93,  10.85,
    10.78,  10.8,   10.72,  10.65,  10.58,  10.35,  10.05,  9.76,   9.19,
    8.72,   7.85,   7.11,   6.45,   5.13,   4.44,   4,      3.91,   3.91,
    3.91,   3.91,   3.91,   3.91,   3.91,   3.91,   3.91,   3.91,   3.91,
    3.91};
// https://arxiv.org/pdf/1609.04979.pdf
double dofHelper::xconts[100] = {};
double dofHelper::contrFermion[100] = {};
double dofHelper::contrBoson[100] = {};
double dofHelper::scontrFermion[100] = {};
double dofHelper::scontrBoson[100] = {};

void dofHelper::Init() {
  int NOrder = 100;
  for (int i = 0; i < 100; i++) {
    xconts[i] = 10.0 - 0.1 * i;

    double *wb1 = new double[NOrder];
    double *xb1 = new double[NOrder];
    cgqf(NOrder, xconts[i], 20, xb1, wb1);

    //Extrapolate on fermions and bosons
    for (int j = 0; j < NOrder; j++) {
      contrFermion[i] += 15.0 / pow(M_PI, 4) * pow(xb1[j], 2) *
                         sqrt(pow(xb1[j], 2) - pow(xconts[i], 2)) * wb1[j] /
                         (exp(xb1[j]) + 1);
      contrBoson[i] += 15.0 / pow(M_PI, 4) * pow(xb1[j], 2) *
                       sqrt(pow(xb1[j], 2) - pow(xconts[i], 2)) * wb1[j] /
                       (exp(xb1[j]) - 1);
      scontrFermion[i] += 15.0 / pow(M_PI, 4) *
                          pow(pow(xb1[j], 2) - pow(xconts[i], 2), 1.5) *
                          wb1[j] / (exp(xb1[j]) + 1);
      scontrBoson[i] += 15.0 / pow(M_PI, 4) *
                        pow(pow(xb1[j], 2) - pow(xconts[i], 2), 1.5) * wb1[j] /
                        (exp(xb1[j]) - 1);
    }
    scontrFermion[i] = (3.0 * contrFermion[i] + scontrFermion[i]) / 4.0;
    scontrBoson[i] = (3.0 * contrBoson[i] + scontrBoson[i]) / 4.0;
    delete[] wb1;
    delete[] xb1;
  }
}

//Look up closest known temperature
int dofHelper::GetIndexClosest(double x, double arr[], int N) {
  int low = 0;
  int high = N;
  int current = (high - low) / 2;
  while ((high - low) / 2 > 0) {
    current = low + (high - low) / 2;
    if (arr[current] >= x)
      low = current;
    else
      high = current;
  }
  return current - 1;
}

//Interpolate to actual temperature x
double dofHelper::Interpolate(double x, double arrx[], double arry[], int N) {
  int index = GetIndexClosest(x, arrx, N);
  if (index == 0 || index >= N - 1)
    return arry[index];
  return arry[index] + (x - arrx[index]) * (arry[index + 1] - arry[index]) /
                           (arrx[index + 1] - arrx[index]);
}

//Get Degees of freedom from Standard Model alone by energy
double dofHelper::dofSMContribution(double T) {
  if (T > TBenchmark[0]) {
    printf("Temperature to high\n");
    return 999999;
  }
  int n = 0;
  for (size_t i = 0; i < 100; i++) {
    if (T <= TBenchmark[i])
      n = i;
  }
  return exp(log(gEeff[n]) + (log(gEeff[n + 1]) - log(gEeff[n])) *
                                 (T - TBenchmark[n]) /
                                 (TBenchmark[n + 1] - TBenchmark[n]));
}

//Get Degees of freedom from Standard Model alone by entropy
double dofHelper::dofsSMContribution(double T) {
  if (T > TBenchmarks[0]) {
    printf("Temperature to high\n");
    return 999999;
  }
  int n = 0;
  for (size_t i = 0; i < 100; i++) {
    if (T <= TBenchmarks[i])
      n = i;
  }
  return exp(log(gSeff[n]) + (log(gSeff[n + 1]) - log(gSeff[n])) *
                                 (T - TBenchmarks[n]) /
                                 (TBenchmarks[n + 1] - TBenchmarks[n]));
}

//Get approximate Degree of Freedom by energy at temperature T
double dofHelper::DOF(double T, Settings *sett) {
  double total = Interpolate(T, TBenchmark, gEeff, 100);
  if (sett->IsFermion_chi)
    total +=
        sett->gchi * Interpolate(sett->mchi / T, xconts, contrFermion, 100);
  else
    total += sett->gchi * Interpolate(sett->mchi / T, xconts, contrBoson, 100);

  if (sett->IsFermion_med)
    total +=
        sett->gmed * Interpolate(sett->mmed / T, xconts, contrFermion, 100);
  else
    total += sett->gmed * Interpolate(sett->mmed / T, xconts, contrBoson, 100);

  return total;
}

//Get approximate Degree of Freedom by entropy at temperature T
double dofHelper::DOFs(double T, Settings *sett) {
  double total = Interpolate(T, TBenchmarks, gSeff, 100);
  if (sett->IsFermion_chi)
    total +=
        sett->gchi * Interpolate(sett->mchi / T, xconts, scontrFermion, 100);
  else
    total += sett->gchi * Interpolate(sett->mchi / T, xconts, scontrBoson, 100);

  if (sett->IsFermion_med)
    total +=
        sett->gmed * Interpolate(sett->mmed / T, xconts, scontrFermion, 100);
  else
    total += sett->gmed * Interpolate(sett->mmed / T, xconts, scontrBoson, 100);

  return total;
}
