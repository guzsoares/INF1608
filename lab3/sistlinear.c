#include <stdio.h>
#include <math.h>

void troca_linhas(int n, double** A, double* b, int linha1, int linha2) {
    double* tempA = A[linha1];
    A[linha1] = A[linha2];
    A[linha2] = tempA;

    double tempB = b[linha1];
    b[linha1] = b[linha2];
    b[linha2] = tempB;
}

void gauss(int n, double** A, double* b, double* x) {

    for (int k = 0; k < n; k++) {
        int maxIndex = k;
        for (int i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > fabs(A[maxIndex][k])) {
                maxIndex = i;
            }
        }
        if (maxIndex != k) {
            troca_linhas(n, A, b, k, maxIndex);
        }

        for (int i = k + 1; i < n; i++) {
            double fator = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= fator * A[k][j];
            }
            b[i] -= fator * b[k];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}