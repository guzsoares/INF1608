#include "raiz.h"
#include <math.h>
#include <stdio.h>

int bissecao(double a, double b, double (*f)(double x), double* r) {
    double c;
    int iteracoes = 0;
    
    if (f(a) * f(b) >= 0) {
        return 0;
    }
    
    while ((b - a) / 2.0 > 0.5e-8) {
        iteracoes++;
        c = (a + b) / 2.0;
        
        if (fabs(f(c)) < 1e-12) {
            *r = c;
            return iteracoes;
        } else if (f(a) * f(c) < 0) {
            b = c;
        } else {
            a = c;
        }
    }
    
    *r = (a + b) / 2.0;
    return iteracoes;
}
