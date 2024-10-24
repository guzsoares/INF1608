#include "ode.h"

#include <stdio.h>
#include <math.h>

static int N;

#define Sf1 10.0542731796122
#define Sf2(t) (1.0 - exp(-10.0*t) / 2.0)

double f1(double t, double y) {
    N++;
    return t * y + pow(t, 3);
}

double f2(double t, double y) {
    N++;
    return 10 * (1 - y);
}

void Avalia (double t0, double t1, double y0, double h, double tol, double (*f) (double t, double y), double sol)
{
  printf("Solucao: %.16f\n", sol);
    
    N = 0;
    double rk = RungeKutta(t0, t1, y0, f, h);
    printf("RungeKutta   %d   %.16f   %g\n", N, rk, fabs(sol - rk));
    
    N = 0;
    double ad = RungeKuttaAdaptativo(t0, t1, y0, f, tol);
    printf("Adaptativo   %d   %.16f   %g\n", N, ad, fabs(sol - ad));
}

int main (void)
{
  double t0 = 0.0;
    double t1_f1 = 2.4;
    double y0_f1 = -1.0;
    double h = 0.001;
    double tol = 1e-12;

    // Evaluate for f1
    Avalia(t0, t1_f1, y0_f1, h, tol, f1, Sf1);

    double t1_f2 = 3.0;
    double y0_f2 = 0.0;

    // Evaluate for f2
    Avalia(t0, t1_f2, y0_f2, h, tol,f2, Sf2(t1_f2));

    return 0;
}
