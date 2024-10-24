#include "interp.h"

#include <math.h>
#include <stdio.h>

#define PI 3.141592653589793

static double f (double x)
{
  return 1.0 / (1.0 + 12.0 * x * x);
}

#define N 10
int main (void)
{
  double xi[N];
  double bi[N];
  
  printf("Coeficientes espacados regularmente:\n");
  for (int i = 2; i <= N; i++) {
    regular(i, -2, 2, xi);
	printf("%d: ", i);
	for (int j = 0; j < i; j++)
	  printf("%.4f ", xi[j]);
	printf("\n");
  } 

  printf("Coeficientes de chebyshev:\n");
  for (int i = 2; i <= N; i++) {
    chebyshev(i, -2, 2, xi);
	printf("%d: ", i);
	for (int j = 0; j < i; j++)
	  printf("%.4f ", xi[j]);
	printf("\n");
  } 
  
  printf("Diff. Divididas de Newton:\n");
  printf("Regulares: ");
  regular(N, -2, 2, xi);
  coeficientes(N,xi,f,bi);
  for (int i=0; i < N; ++i)
    printf("%g ",bi[i]);
  printf("\n");
  
  printf("Chebyshev: ");
  chebyshev(N, -2, 2, xi);
  coeficientes(N,xi,f,bi);
  for (int i=0; i < N; ++i)
    printf("%g ",bi[i]);
  printf("\n");

  printf("Avalia funcao 1/(1+12x^2): \n");
  printf("Com coeficientes regulares: \n");
  regular(N, -2, 2, xi);
  coeficientes(N,xi,f,bi);
  for (float x=-2.0; x < 2.0001; x+=0.1) {
    double y = avalia(N,xi,bi,x);
    double yf = f(x);
    printf("%.1f %f %f %g\n",x,y,yf,fabs(y-yf));
  }
  printf("Com coeficientes chebyshev: \n");
  chebyshev(N, -2, 2, xi);
  coeficientes(N,xi,f,bi);
  for (float x=-2.0; x < 2.0001; x+=0.1) {
    double y = avalia(N,xi,bi,x);
    double yf = f(x);
    printf("%.1f %f %f %g\n",x,y,yf,fabs(y-yf));
  }
  
  return 0;
}
