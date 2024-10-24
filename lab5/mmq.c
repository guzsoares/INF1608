#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matriz.h"
#include "sistlinear.h"

double mmq(int m, int n, double** A, double* b, double* x) {
    double** AT = mat_cria(n, m); 
    double** ATA = mat_cria(n, n); 
    double* ATb = (double*) malloc(n * sizeof(double)); 

    mat_transposta(m, n, A, AT);

    mat_multm(n, m, n, AT, A, ATA);

    mat_multv(n, m, AT, b, ATb);

    gauss(n, ATA, ATb, x);

    double* r = (double*) malloc(m * sizeof(double));
    double* Ax = (double*) malloc(m * sizeof(double));
    mat_multv(m, n, A, x, Ax);
    for (int i = 0; i < m; i++) {
        r[i] = b[i] - Ax[i];
    }

    double norma2 = 0.0;
    for (int i = 0; i < m; i++) {
        norma2 += r[i] * r[i];
    }
    norma2 = sqrt(norma2);

    mat_libera(n, AT);
    mat_libera(n, ATA);
    free(ATb);
    free(r);
    free(Ax);

    return norma2;
}

double ajuste_parabola(int n, double* px, double* py, double* a, double* b, double* c) {
    double** A = mat_cria(n, 3);
    for (int i = 0; i < n; i++) {
        A[i][0] = 1.0;
        A[i][1] = px[i]; 
        A[i][2] = px[i] * px[i]; 
    }

    double* coeficientes = (double*) malloc(3 * sizeof(double));
    double norma2 = mmq(n, 3, A, py, coeficientes);

    *a = coeficientes[0];
    *b = coeficientes[1];
    *c = coeficientes[2];

    mat_libera(n, A);
    free(coeficientes);

    return norma2; 
}

double ajuste_cubica(int n, double* px, double* py, double* a, double* b, double* c, double* d) {
    double** A = mat_cria(n, 4);
    for (int i = 0; i < n; i++) {
        A[i][0] = 1.0;  
        A[i][1] = px[i]; 
        A[i][2] = px[i] * px[i];  
        A[i][3] = px[i] * px[i] * px[i]; 
    }

    double* coeficientes = (double*) malloc(4 * sizeof(double));
    double norma2 = mmq(n, 4, A, py, coeficientes);

    *a = coeficientes[0];
    *b = coeficientes[1];
    *c = coeficientes[2];
    *d = coeficientes[3];

    mat_libera(n, A);
    free(coeficientes);

    return norma2; 
}

double ajuste_exponencial_exp(int n, double* px, double* py, double* a, double* b) {
    double* log_py = (double*) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        log_py[i] = log(py[i]);
    }

    double** A = mat_cria(n, 2);
    for (int i = 0; i < n; i++) {
        A[i][0] = 1.0;   
        A[i][1] = px[i];
    }

    double* coeficientes = (double*) malloc(2 * sizeof(double));
    double norma2 = mmq(n, 2, A, log_py, coeficientes);

    *a = exp(coeficientes[0]);
    *b = coeficientes[1];

    mat_libera(n, A);
    free(log_py);
    free(coeficientes);

    return norma2; 
}


