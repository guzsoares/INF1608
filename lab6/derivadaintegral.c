#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "derivadaintegral.h"

double derivada(double (*f)(double x), double x, double h) {
    return (f(x + h) - f(x - h)) / (2 * h);
}

double simpson(double (*f)(double), double a, double b, int n) {
    if (n == 1) {
        double h = (b - a) /  n;
        double x0 = a;
        double x1 = (a + h/2.0);
        double x2 = a + h;
        double s = (h / 6.0) * (f(x0) + 4.0 * f(x1) + f(x2));
        return s;
    } else {
        double h = (b - a) / n;
        double s = f(a) + f(b);
        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            if (i % 2 == 0) {
                s += 2.0 * f(x);
            } else {
                s += 4.0 * f(x);
            }
        }
        s *= h / 3.0;
        return s;
    }
}

double simpsoncomposite(double (*f)(double), double a, double b) {
    double h = (b - a) / 2.0;
    return (h / 3.0) * (f(a) + 4 * f(a + h) + f(b));
}

double simpsonadaptativo(double (*f)(double), double a, double b, double tol) {
    double m = (a + b) / 2.0;
    double S = simpsoncomposite(f, a, b);
    double S1 = simpsoncomposite(f, a, m);
    double S2 = simpsoncomposite(f, m, b);
    
    if (fabs(S1 + S2 - S) < 15 * tol) {
        return S1 + S2 + (S1 + S2 - S) / 15.0;
    } else {
        return simpsonadaptativo(f, a, m, tol / 2.0) + simpsonadaptativo(f, m, b, tol / 2.0);
    }
}

double quadraturagauss2(double (*f)(double), double a, double b) {
    double m = (a + b) / 2.0;
    double h = (b - a) / 2.0;
    double x0 = -0.5773502692; 
    double x1 = 0.5773502692; 
    double c0 = 1.0;            
    double c1 = 1.0;            
    return h * (c0 * f(h * x0 + m) + c1 * f(h * x1 + m));
}
