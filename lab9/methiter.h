#ifndef METITER_H
#define METITER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int gaussseidel(int n, double** A, double* b, double* x, double tol);

int sor_gaussseidel(int n, double** A, double* b, double* x, double tol, double w);

#endif 
