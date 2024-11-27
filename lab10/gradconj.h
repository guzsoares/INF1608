// gradconj.h
#ifndef GRADCONJ_H
#define GRADCONJ_H

int gradconj(int n, double** A, double* b, double* x, double tol);
int gradconj_jacobi(int n, double** A, double* b, double* x, double tol);

#endif
