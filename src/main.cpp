#include "legendre_rule.h"
#include "matrixelements.h"
#include "relic.h"
#include "thermalAverage.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <omp.h>
#include <random>
#include <stdio.h>
#include <vector>

//Scan parameterspace at random coupling constants and masses
void scanRadom() {
  dofHelper::Init();
  vector<mt19937_64> gens;
  for (int i = 0; i < omp_get_max_threads(); i++) {
    gens.push_back(mt19937_64(i + time(0)));
  }

	//Setup random point in parameterspace
  std::uniform_real_distribution<double> Massdistribution(-4.0, -0.0);
  std::uniform_real_distribution<double> CouplingDistribution(-6.0, 0.0);
#pragma omp parallel for
  for (size_t i = 0; i < 100000; i++) {
    Settings sett;
    sett.mchi = pow(10.0, Massdistribution(gens[omp_get_thread_num()]));
    sett.cchi = 1.0;
    sett.csm = pow(10.0, CouplingDistribution(gens[omp_get_thread_num()]));
    sett.mmed = sett.mchi * 3.0;
    relic rel(&sett);
    thermalAverage *thermAvg = new thermalAverage(&sett);
    double xf = rel.CalculateXf();
    sett.x0 = xf - 5.0;
    sett.x1 = xf + 15.0;
    double relicResult = rel.CalculateRelicDensity2(0, 0);
    //Write results to file
#pragma omp critical
    {
      FILE *pFileeq;
      printf("Thread: %d 	%E	%E	%E	%E	%E\n",
             omp_get_thread_num(), sett.mchi, sett.mmed, sett.csm, xf,
             relicResult);
      pFileeq = fopen("GLS_YukawaType.txt", "a");
      fprintf(pFileeq, "%E	%E	%E	%E	%E\n", sett.mchi,
              sett.mmed, sett.csm, xf, relicResult);
      fclose(pFileeq);
    }
    delete thermAvg;
  }
}

int main(int argc, char *argv[]) {
  scanRadom();
  return 0;
}
