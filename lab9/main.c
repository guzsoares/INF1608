#include "methiter.h"
#include "matriz.h"
#include "vetor.h"

#include <stdio.h>

static void printv (int n, double* v)
{
  for (int i=0; i<n; ++i)
    printf("%.6g\n",v[i]);
}
static void zerov (int n, double* v)
{
  for (int i=0; i<n; ++i)
    v[i] = 0.0;
}

int main (void)
{
  // Teste 1
  {
  double** A = mat_cria(2,2);
  A[0][0] = 3;
  A[0][1] = 1;
  A[1][0] = 1;
  A[1][1] = 2;
  double b[] = {5,5};
  
  double x[2];
  zerov(2,x);
  int n = gaussseidel(2,A,b,x,1e-7);
  printf("iter: %d\n",n);
  printv(2,x);

  zerov(2,x);
  n = sor_gaussseidel(2,A,b,x,1e-7,1.1);
  printf("iter: %d\n",n);
  printv(2,x);

  mat_libera(2,A);
  }
  // Teste 2
  {
  double** A = mat_cria(3,3);
  A[0][0] =  3;
  A[0][1] =  1;
  A[0][2] = -1;
  A[1][0] =  2;
  A[1][1] =  4;
  A[1][2] =  1;
  A[2][0] = -1;
  A[2][1] =  2;
  A[2][2] =  5;
  double b[] = {4, 1, 1};
  
  double x[3];
  zerov(3,x);
  int n = gaussseidel(3,A,b,x,1e-7);
  printf("iter: %d\n",n);
  printv(3,x);

  zerov(3,x);
  n = sor_gaussseidel(3,A,b,x,1e-7,1.1);
  printf("iter: %d\n",n);
  printv(3,x);

  mat_libera(3,A);
  }
  // Teste 3
  {
  int N = 6;
  double** A = mat_cria(N,N);
  for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
      if (i==j) 
        A[i][j] = 3.0;
      else if (i+j == N-1)
        A[i][j] = 0.5;
      else
        A[i][j] = 0.0;
	
  for (int i=0; i<N-1; ++i)
    A[i][i+1] = -1;
  for (int i=1; i<N; ++i)
    A[i][i-1] = -1;

  double b[N];
  b[0] = b[N-1] = 2.5;
  b[1] = b[N-2] = 1.5;
  for (int i=2; i<N-2; ++i)
    b[i] = 1.0;

  double x[N];
  zerov(N,x);
  int n = gaussseidel(N,A,b,x,1e-7);
  printf("iter: %d\n",n);
  printv(N,x);

  zerov(N,x);
  n = sor_gaussseidel(N,A,b,x,1e-7,1.1);
  printf("iter: %d\n",n);
  printv(N,x);

  mat_libera(N,A);
  }
  return 0;
}
