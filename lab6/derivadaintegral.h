#ifndef DERIVADAINTEGRAL_H
#define DERIVADAINTEGRAL_H

double derivada(double (*f)(double x), double x, double h);
double simpson(double (*f)(double), double a, double b, int n);
double simpsonadaptativo(double (*f)(double), double a, double b, double tol);
double quadraturagauss2(double (*f)(double), double a, double b);

#endif // DERIVADAINTEGRAL_H
