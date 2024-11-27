#include <stdio.h>
#include <math.h>
#include "otimizacao.h"

double func1(double x) {
    return x * x + sin(x);
}

double func2(double x) {
    return pow(x, 6) - 11 * pow(x, 3) + 17 * x * x - 7 * x + 1;
}

int fevals;

double counted_func1(double x) {
    fevals++;
    return func1(x);
}

double counted_func2(double x) {
    fevals++;
    return func2(x);
}

int main() {
    double xmin;
    int iter;

    printf("MSA:\n");

    printf("x*x + sin(x)\n");
    fevals = 0;
    iter = msa(-0.5, 0.5, counted_func1, 1e-5, &xmin);
    printf("msa  [a =  -0.5, b =   0.5]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    fevals = 0;
    iter = msa(-1.5, 1.5, counted_func1, 1e-5, &xmin);
    printf("msa  [a =  -1.5, b =   1.5]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    // Testes do Método da Seção Áurea (MSA) para func2
    printf("x^6 - 11 x^3 + 17 x^2 - 7 x + 1\n");
    fevals = 0;
    iter = msa(-0.5, 0.5, counted_func2, 1e-5, &xmin);
    printf("msa  [a =  -0.5, b =   0.5]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    fevals = 0;
    iter = msa(-0.7, 0.7, counted_func2, 1e-5, &xmin);
    printf("msa  [a =  -0.7, b =   0.7]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    fevals = 0;
    iter = msa(1, 2, counted_func2, 1e-5, &xmin);
    printf("msa  [a =     1, b =     2]: it = %2d, xmin = %10.5f, fevals: %2d\n", iter, xmin, fevals);

    fevals = 0;
    iter = msa(1, 8, counted_func2, 1e-5, &xmin);
    printf("msa  [a =     1, b =     8]: it = %2d, xmin = %10.5f, fevals: %2d\n\n", iter, xmin, fevals);

    printf("MIPS:\n");

    printf("x*x + sin(x)\n");
    fevals = 0;
    iter = mips(0, 0.5, counted_func1, 1e-6, &xmin);
    printf("mips [r =     0, d =   0.5]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    fevals = 0;
    iter = mips(1.5, 1.0, counted_func1, 1e-6, &xmin);
    printf("mips [r =   1.5, d =     1]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    fevals = 0;
    iter = mips(1.2, 0.2, counted_func1, 1e-6, &xmin);
    printf("mips [r =   1.2, d =   0.2]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    printf("x^6 - 11 x^3 + 17 x^2 - 7 x + 1\n");
    fevals = 0;
    iter = mips(0, 0.5, counted_func2, 1e-6, &xmin);
    printf("mips [r =     0, d =   0.5]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    fevals = 0;
    iter = mips(0.6, 0.01, counted_func2, 1e-6, &xmin);
    printf("mips [r =   0.6, d =  0.01]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    fevals = 0;
    iter = mips(1, 0.5, counted_func2, 1e-6, &xmin);
    printf("mips [r =     1, d =   0.5]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    fevals = 0;
    iter = mips(1.4, 0.1, counted_func2, 1e-6, &xmin);
    printf("mips [r =   1.4, d =   0.1]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    fevals = 0;
    iter = mips(1, 1e-7, counted_func2, 1e-6, &xmin);
    printf("mips [r =     1, d = 1e-07]: it = %2d, xmin = %10.6f, fevals: %2d\n", iter, xmin, fevals);

    return 0;
}
