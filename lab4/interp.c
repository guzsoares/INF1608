#include <math.h>
#include "interp.h"
#include "matriz.h"

#define PI 3.141592653589793

void regular(int n, double a, double b, double* xi) {
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) {
        xi[i] = a + i * step;
    }
}

void chebyshev(int n, double a, double b, double* xi) {
    for (int i = 0; i < n; i++) {
        double beta = (2 * i + 1) * PI / (2 * n);
        xi[i] = 0.5 * (a + b) + 0.5 * (b - a) * cos(beta);
    }
}

void coeficientes(int n, double* xi, double (*f)(double), double* bi) {
    double **matriz = mat_cria(n,n);

    for (int i = 0; i < n; i++) {
        matriz[i][0] = f(xi[i]);
    }

    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            matriz[i][j] = (matriz[i + 1][j - 1] - matriz[i][j - 1]) / (xi[i + j] - xi[i]);
        }
    }

    for (int i = 0; i < n; i++) {
        bi[i] = matriz[0][i];
    }

    mat_libera(n, matriz);
}



double avalia(int n, double* xi, double* bi, double x) {
    double resultado = bi[n - 1]; 
    for (int i = n - 2; i >= 0; i--) {
        resultado = resultado * (x - xi[i]) + bi[i]; 
    }
    return resultado;
}


