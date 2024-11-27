#include "gradconj.h"
#include "matriz.h"
#include "vetor.h"

#include <stdio.h>

static void printm (int n, double** A)
{
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j)
      printf("%10.6g  ",A[i][j]);
    printf("\n");
  }
  printf("\n");
}

static void printv (int n, double* v)
{
  for (int i=0; i<n; ++i)
    printf("%.6g  ",v[i]);
  printf("\n");
}

static void zerov (int n, double* v)
{
  for (int i=0; i<n; ++i)
    v[i] = 0.0;
}

#define GRAD_CONJ_TOL 1e-7
static void GCTest (int n, double** A, double* b)
{
  int it = 0;
  double* x = vet_cria(n);

  // descomentar para imprimir matriz A e vetor b
  //printf("A:\n");
  //mat_imprime(n, n, A);
  //printf("b: ");
  //printv(n, b);

  zerov(n, x);
  it = gradconj(n, A, b, x, GRAD_CONJ_TOL);
  printf("GC iter: %d\n", it);
  printv(n, x);

  zerov(n, x);
  it = gradconj_jacobi(n, A, b, x, GRAD_CONJ_TOL);
  printf("GCJC iter: %d\n", it);
  printv(n, x);

  vet_libera(x);
  printf("--------------------------------\n");
}

int main (void)
{
  { // result x: [1 ... 1]
    int n = 10;
    double** A = mat_cria(n, n);
	double* b = vet_cria(n);

    // initialize matrix	
    for (int i=0; i<n; ++i)
      for (int j=0; j<n; ++j)
        A[i][j] = 0;
	
	for (int i=0; i<n; ++i) {
	  A[i][i] = i + 2;
	  if (i+1 < n) {
		A[i][i + 1] = 0.5;
        A[i + 1][i] = 0.5;	  
	  }
	  if (i+2 < n) {
		A[i][i + 2] = 0.5;
        A[i + 2][i] = 0.5;	  
	  }
    }
	
	// create b vector from x
	double* xe = vet_cria(n);
    for (int i=0; i<n; ++i) xe[i] = 1.0;
    mat_multv(n, n, A, xe, b);

	GCTest(n, A, b);
	vet_libera(xe);
	vet_libera(b);
	mat_libera(n, A);
  }
  { // result x: [1 ... 1]
    int n = 10;
    double** A = mat_cria(n, n);
	double* b = vet_cria(n);

    // initialize matrix	
    for (int i=0; i<n; ++i)
      for (int j=0; j<n; ++j)
        A[i][j] = 0;
	
	for (int i=0; i<n; ++i) {
	  A[i][i] = i + 2;
	  if (i+1 < n) {
		A[i][i + 1] = 0.5;
        A[i + 1][i] = 0.5;	  
	  }
	  if (i+2 < n) {
		A[i][i + 2] = 0.5;
        A[i + 2][i] = 0.5;	  
	  }
	  if (2*i < n) {
		A[i][2*i] = 0.5;
        A[2*i][i] = 0.5;	  
	  }
    }
	
	// create b vector from x
	double* xe = vet_cria(n);
    for (int i=0; i<n; ++i) xe[i] = 1.0;
    mat_multv(n, n, A, xe, b);

	GCTest(n, A, b);
	vet_libera(xe);
	vet_libera(b);
	mat_libera(n, A);
  }
  { // result x: [1 ... 1]
    int n = 10;
    double** A = mat_cria(n, n);
	double* b = vet_cria(n);

    // initialize matrix	
    for (int i=0; i<n; ++i)
      for (int j=0; j<n; ++j)
        A[i][j] = 0;
	
	for (int i=0; i<n; ++i) {
	  A[i][i] = i + 2;
	  if (i+2 < n) {
		A[i][i + 2] = 0.5;
        A[i + 2][i] = 0.5;	  
	  }
	  if (i+4 < n) {
		A[i][i + 4] = 0.5;
        A[i + 4][i] = 0.5;	  
	  }
	  if (i < n - n/100) { 
	    A[i][n-1] = -0.1;
	    A[n-1][i] = -0.1;
      }
    }
	
	// create b vector from x
	double* xe = vet_cria(n);
    for (int i=0; i<n; ++i) xe[i] = 1.0;
    mat_multv(n, n, A, xe, b);

	GCTest(n, A, b);
	vet_libera(xe);
	vet_libera(b);
	mat_libera(n, A);
  }
  { // result x: [1 ... 1]
    // completar com matriz da questao 2
  }
    
  return 0;
}
