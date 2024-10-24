#include "taylor.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.141592653589793

int main (void)
{
  printf("Teste Q1 1: f(x)=e^x, x0=0, x=1\n");
  {
	double c[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
    double v = exp(1.0);
    for (int n=1; n<=5; ++n) { // teste com 1 até 5 termos
      double v1 = avalia_taylor(n,c,0.0,1.0);
      double d = v-v1;
      printf("%d: %.9f   %.9f\n",n,v1,d);
    }
  }
  printf("Teste Q1 2: f(x)=sin x, x0=0, x=PI/2\n");
  {
	double c[10] = {0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0};
    double v = sin(PI/2.0);
    for (int n=1; n<=10; ++n) { // teste com 1 até 10 termos
      double v1 = avalia_taylor(n,c,0.0,PI/2.0);
      double d = v-v1;
      printf("%d: %.9f   %.9f\n",n,v1,d);
    }
  }
  printf("Teste Q2 1: x = 0.4 * PI\n");
  {
    double x = 0.4 * PI;
    printf("Seno:\n");
    double v = sin(x);
    printf("x = %.16g, f(x) = %.16g\n",x,v);
    for (int n=1; n<=20; ++n) { // aproximando 1 até 20 termos
      double v1 = avalia_seno(n,x);
      double d = v-v1;
      printf("%d: %.16g   %.16g\n",n,v1,d);
    }
  }
  printf("Teste Q2 2: x = [PI/12, ... , PI/3]\n");
  {
    double x[5] = {PI/12.0, PI/9.0, PI/7.0, PI/4.0, PI/3.0};
    	  
    printf("Seno:\n");
    for (int i=0; i<5; ++i) {
      double v = sin(x[i]);
      printf("x = %.9f,  f(x) = %.9f\n",x[i],v);
      for (int n=1; n<=10; ++n) { // aproximando 1 até 10 termos
        double v1 = avalia_seno(n,x[i]);
        double d = v-v1;
        printf("%d: %.9f   %.9f\n",n,v1,d);
      }
    }
  }
  return 0;
}
