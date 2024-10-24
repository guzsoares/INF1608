#include "taylor.h"
#include <stdlib.h>

double avalia_taylor (int n, double *c, double x0, double x){
    double result = 0.0;
    double pow = 1.0;
    double factorial = 1.0;

    for (int i = 0; i < n; i++) {
        result += c[i] * pow / factorial;
        pow *= (x - x0);
        factorial *= (i + 1);
    }

    return result;
}

double avalia_seno (int n, double x){
    const int variation[] = {0, 1, 0, -1};
    double coefficients[n];

    for (int i = 0; i < n; i++){
        int pos = i % 4;
        coefficients[i] = variation[pos];
    }

    return avalia_taylor(n, coefficients, 0.0, x);
}