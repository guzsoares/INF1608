#include <math.h>
#include "otimizacao.h"

#define MAX_ITER 50
#define GOLDEN_RATIO ((sqrt(5) - 1) / 2)
#define TOLERANCE 1e-10


int parar(double x1, double x2, double tol) {
    return ((x2 - x1) / 2.0) <= tol;
}

int msa(double a, double b, double (*f)(double x), double tol, double *xmin) {
    int fevals_msa = 0;
    double x1, x2, fx1, fx2;
    int iter = 0;

    x1 = a + (1 - GOLDEN_RATIO) * (b - a);
    x2 = a + GOLDEN_RATIO * (b - a);
    fx1 = f(x1);
    fx2 = f(x2);
    fevals_msa += 2; 

    while (!parar(a, b, tol) && iter < MAX_ITER) {
        iter++;
        if (fx1 < fx2) {
            b = x2;
            x2 = x1;
            fx2 = fx1;
            x1 = a + (1 - GOLDEN_RATIO) * (b - a);
            fx1 = f(x1);
            fevals_msa++;
        } else {
            a = x1;
            x1 = x2;
            fx1 = fx2;
            x2 = a + GOLDEN_RATIO * (b - a);
            fx2 = f(x2);
            fevals_msa++;
        }
    }

    *xmin = (a + b) / 2.0;
    return iter;
}

int mips(double r, double delta, double (*f)(double x), double tol, double *xmin) {
    double s = r - delta, t = r + delta;
    double x, fx, fr, fs, ft, denom;
    int iter = 0;

    fr = f(r);
    fs = f(s);
    ft = f(t);

    while (iter < MAX_ITER) {
        if (iter >= 3 && fabs(fs - ft) <= tol) {
            *xmin = (s + t) / 2.0;
            return iter;
        }

        denom = 2 * ((s - r) * (ft - fs) - (fs - fr) * (t - s));
        if (fabs(denom) < 1e-10) {
            x = (r + s + t) / 3.0;
        } else {
            x = ((r + s) / 2.0) - ((fs - fr) * (t - r) * (t - s) / denom);
        }

        fx = f(x);
        r = s;
        fr = fs;
        s = t;
        fs = ft;
        t = x;
        ft = fx;
        iter++;
    }

    return 0;
}