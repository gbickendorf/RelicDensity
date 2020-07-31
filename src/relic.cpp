#include "relic.h"

#define VERBOSE_FILE 0


relic::relic(Settings *sett) {
  settings = sett;
  thermAvg = new thermalAverage(sett);
  sett->wmed = MatrixElements::Width(sett);
}

relic::~relic() { delete thermAvg; }

void relic::Init() {
	//Initialise rough and fine stepcount
  int Nfine = (int)(settings->x1 - settings->xmin) / settings->dxfine;
  int Nrough = (int)(settings->xinf - settings->x1) / settings->dxrough;

  x.resize(Nfine + Nrough);
  w.resize(Nfine + Nrough);
  y.resize(Nfine + Nrough);

	//Initialise walues only dependent on x
#pragma omp parallel for shared(x, w, y, settings)
  for (int i = 0; i < Nfine; i++) {
    x[i] = settings->xmin + settings->dxfine * i;
    w[i] = thermAvg->thermalAverageCrossection(settings->mchi / x[i]);
    y[i] = Yeq(x[i]);
  }
#pragma omp parallel for shared(x, w, y, settings)
  for (int i = Nfine; i < Nfine + Nrough; i++) {
    x[i] = settings->x1 + settings->dxrough * (i - Nfine);
    w[i] = thermAvg->thermalAverageCrossection(settings->mchi / x[i]);
    y[i] = 0.0;
  }
}

//Calculate the freezeout temperature X_f
double relic::CalculateXf() {
  int MaxGuard = 100;
  double xtrial = 100.0;
  double dx = 1.0;
  double currVal = -100.0;
  do {
    double thrmav =
        thermAvg->thermalAverageCrossection(settings->mchi / xtrial);
    currVal =
        log(settings->mchi / (2.0 * pow(M_PI, 3)) *
            sqrt(45.0 * mpl * mpl /
                 (2.0 * dofHelper::DOF(settings->mchi / xtrial, settings))) *
            thrmav / sqrt(xtrial)) -
        xtrial;

    if (currVal >= 0) {
      xtrial += dx;
      dx /= 10;
    }
    xtrial -= dx;
    MaxGuard--;
  } while (MaxGuard > 0 && dx > settings->eps_xf);
  return xtrial;
}

//Calculate equilibrium particel density normalised by entropy
double relic::Yeq(double x) {
  return 45.0 / (4 * pow(M_PI, 4)) * sqrt(M_PI / 2) * settings->gchi /
         dofHelper::DOFs(settings->mchi / x, settings) * pow(x, 1.5) * exp(-x);
}

//Calculate the derivative occuring in the diff. equation
double relic::fprime(int i, double f) {
  return -sqrt(M_PI /
               (45.0 * dofHelper::DOF(settings->mchi / x[i], settings))) *
         dofHelper::DOFs(settings->mchi / x[i], settings) * mpl *
         settings->mchi / (x[i] * x[i]) * 2.0 * w[i] *
         (pow(f, 2.0) - pow(y[i], 2.0));
}

//Calculate the relic densisty while setting the lowest temperature X_0 automatically
double relic::CalculateRelicDensityAutoX0(int writeFiles, int verbose) {
  double relicBefore = 0.0;
  double relicNow = CalculateRelicDensity(writeFiles, verbose);
  do {
    if (settings->x0 < settings->xmin) {
      printf("Too close to x0 = 0 \n");
      break;
    }
    settings->x0 -= 1;
    relicBefore = relicNow;
    relicNow = CalculateRelicDensity(writeFiles, verbose);
    if (verbose)
      printf("%f 		%E\n", settings->x0, relicNow);
  } while (abs((relicBefore - relicNow) / relicNow) > settings->eps_ode ||
           isinf(relicNow));

  return relicNow;
}

//Calculate the relic densisty
double relic::CalculateRelicDensity(int writeFiles, int verbose) {
#if VERBOSE_FILE > 0
  FILE *pFileeq;
  pFileeq = fopen("RelicEq.csv", "w");
  FILE *pFile;
  pFile = fopen("RelicTrial.csv", "w");
#endif
  double f = Yeq(x[0]);
  int Nfine = (int)(settings->x1 - settings->xmin) / settings->dxfine;
  int Nrough = (int)(settings->xinf - settings->x1) / settings->dxrough;
  for (int i = 1; i < Nfine + Nrough; i++) {
#if VERBOSE_FILE > 0
    fprintf(pFile, "%E,%E\n", x[i], f);
    fprintf(pFileeq, "%E,%E\n", x[i], Yeq(x[i]));
#endif
    if (isinf(f)) {
      break;
    }
    if (i % 100 == 0 && verbose)
      printf("%.3f %%\n", (double)i / (Nfine + Nrough) * 100.0);
    if (fprime(i, f) == 0.0) {
      printf("No changes anymore\n");
      break;
    }
    f += fprime(i, f) * (x[i] - x[i - 1]);
  }
#if VERBOSE_FILE > 0
  fclose(pFile);
#endif

  return f * settings->mchi * 2889.2 / 1.05E-5;
}

//Derivative for the gsl-routine for diff. equations
double relic::fprime_gsl(double xrun, double f) {
  return -sqrt(M_PI /
               (45.0 * dofHelper::DOF(settings->mchi / xrun, settings))) *
         dofHelper::DOFs(settings->mchi / xrun, settings) * mpl *
         settings->mchi / (xrun * xrun) * 2.0 * thermAvg->thermalAverageCrossection(settings->mchi / xrun) *
         (pow(f, 2.0) - pow(Yeq(xrun), 2.0));
  ;
}

int odefunc1(double xrun, const double y[], double f[], void *params) {
  relic *rel = static_cast<relic *>(params);
  f[0] = rel->fprime_gsl(xrun, y[0]);
  return GSL_SUCCESS;
}

//Calculate relic denisity this time with gsl solvers
double relic::CalculateRelicDensity2(int writeFiles, int verbose) {
  size_t dim = 1;
  gsl_odeiv2_system sys = {odefunc1, NULL, dim, this};

  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
      &sys, gsl_odeiv2_step_rk8pd, 1e-10, 1e-10, 0.0);
  int i;
  double x0 = settings->x0,
         xf = 1000;
  double xrun = x0;
  double f[1] = {Yeq(x0)}; /* initial value */
  int N = 100;
  for (i = 1; i <= N; i++) {
    double xi = x0 + 1.0 * i * (xf - x0) / N;
    int status = gsl_odeiv2_driver_apply(d, &xrun, xi, f);
    if (status != GSL_SUCCESS) {
      return INFINITY;
      break;
    }
  }
  gsl_odeiv2_driver_free(d);

  return f[0] * settings->mchi * 2889.2 / 1.05E-5;
}
