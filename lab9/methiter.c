#include "methiter.h"
#include <math.h>
#include <stdio.h>

double calc_norma(int n, double* r) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += r[i] * r[i];
    }
    return sqrt(sum);
}

int gaussseidel(int n, double** A, double* b, double* x, double tol) {
    double* r = (double*)malloc(n * sizeof(double));
    int iter = 0;

    while (1) {
        iter++;
        for (int i = 0; i < n; i++) {
            double sigma = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sigma += A[i][j] * x[j];
                }
            }
            x[i] = (b[i] - sigma) / A[i][i];
        }

        for (int i = 0; i < n; i++) {
            r[i] = b[i];
            for (int j = 0; j < n; j++) {
                r[i] -= A[i][j] * x[j];
            }
        }

        if (calc_norma(n, r) <= tol) {
            free(r);
            return iter;
        }
    }
}

int sor_gaussseidel(int n, double** A, double* b, double* x, double tol, double w) {
    double* r = (double*)malloc(n * sizeof(double));
    double* x_old = (double*)malloc(n * sizeof(double));
    int iter = 0;

    while (1) {
        iter++;

        for (int i = 0; i < n; i++) {
            x_old[i] = x[i];
        }

        for (int i = 0; i < n; i++) {
            double sigma = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sigma += A[i][j] * x[j];
                }
            }
            double x_new = (b[i] - sigma) / A[i][i];
            x[i] = (1 - w) * x_old[i] + w * x_new;
        }

        for (int i = 0; i < n; i++) {
            r[i] = b[i];
            for (int j = 0; j < n; j++) {
                r[i] -= A[i][j] * x[j];
            }
        }

        if (calc_norma(n, r) <= tol) {
            free(r);
            free(x_old);
            return iter;
        }
    }
}