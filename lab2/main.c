
#include "raiz.h"
#include <math.h>
#include <stdio.h>

double func1(double x) {
    return x * x * x + 2 * x * x - 2;
}

double func2(double x) {
    return x * x * x * x + x + 1;
}

double func3(double x) {
    double g = 9.8;
    double c = 15.0;
    double t = 9.0;
    double v = 35.0;
    return (g * x / c) * (1 - exp(-c * t / x)) - v;
}

int main(void) {
    double raiz;
    int iteracoes;
    int avaliacoes;
    printf("Q1 e Q2 de -3 a +3\n");
    iteracoes = bissecao(-3, 3, func1, &raiz);
    avaliacoes = iteracoes + 2;
    if (iteracoes > 0) {
        printf("Q1: %.9f it: %d av: %d\n", raiz, iteracoes, avaliacoes);
    } else {
        printf("Q1: Função não convergiu it: %d av: %d\n", iteracoes, avaliacoes);
    }

    iteracoes = bissecao(-3, 3, func2, &raiz);
    avaliacoes = iteracoes + 2;
    if (iteracoes > 0) {
        printf("Q2: %.9f it: %d av: %d\n", raiz, iteracoes, avaliacoes);
    } else {
        printf("Q2: Função não convergiu it: %d av: %d\n", iteracoes, avaliacoes);
    }

    iteracoes = bissecao(-5, 5, func1, &raiz);
    avaliacoes = iteracoes + 2;
    if (iteracoes > 0) {
        printf("Q1 de -5 a +5: %.9f it: %d av: %d\n", raiz, iteracoes, avaliacoes);
    } else {
        printf("Q1 de -5 a +5: Função não convergiu it: %d av: %d\n", iteracoes, avaliacoes);
    }

    iteracoes = bissecao(0, 150, func3, &raiz);
    avaliacoes = iteracoes + 2;
    if (iteracoes > 0) {
        printf("Q3 de 0 a 150: %.9f it: %d av: %d\n", raiz, iteracoes, avaliacoes);
    } else {
        printf("Q3 de 0 a 150: Função não convergiu it: %d av: %d\n", iteracoes, avaliacoes);
    }

    iteracoes = bissecao(0, 250, func3, &raiz);
    avaliacoes = iteracoes + 2;
    if (iteracoes > 0) {
        printf("Q3 de 0 a 250: %.9f it: %d av: %d\n", raiz, iteracoes, avaliacoes);
    } else {
        printf("Q3 de 0 a 250: Função não convergiu it: %d av: %d\n", iteracoes, avaliacoes);
    }

    iteracoes = bissecao(0, 50, func3, &raiz);
    avaliacoes = iteracoes + 2;
    if (iteracoes > 0) {
        printf("Q3 de 0 a 50: %.9f it: %d av: %d\n", raiz, iteracoes, avaliacoes);
    } else {
        printf("Q3 de 0 a 50: Função não convergiu it: %d av: %d\n", iteracoes, avaliacoes);
    }

    return 0;
}
