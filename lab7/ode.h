#ifndef _RUNGE_KUTTA_H_
#define _RUNGE_KUTTA_H_

#ifndef MIN
#define MIN(a, b) (a < b ? a : b)
#endif

double RungeKutta (double t0, double t1, double y0, double (*f) (double t, double y), double h);
double RungeKuttaAdaptativo (double t0, double t1, double y0, double (*f) (double t, double y), double tol);

#endif