#include "matriz.h"
#include "sistlinear.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
void fatoracao (int n, double** a, int* p)
{
  for (int i=0; i<n; ++i)
    p[i] = i;
  for (int j=0; j<n-1; ++j) {
    // find piv�
    int m = j;
    for (int i=j+1; i<n; ++i)
      if (fabs(a[i][j]) > fabs(a[m][j]))
        m = i;
    // swap lines: j <-> m
    for (int k=0; k<n; ++k) {
      double t = a[j][k];
      a[j][k] = a[m][k];
      a[m][k] = t;
    }
    // register permutation
    int t = p[j];
    p[j] = p[m];
    p[m] = t;
    // elimination
    for (int i=j+1; i<n; ++i) {
      double f = a[i][j]/a[j][j];
      for (int k=j+1; k<n; ++k)
        a[i][k] -= f*a[j][k];
      a[i][j] = f;
    }
  }
}

void substituicao (int n, double** a, int* p, double* b, double* x)
{
  // forward substitution
  for (int i=0; i<n; ++i) {
    double s = 0;
    for (int j=0; j<i; ++j) 
      s += a[i][j]*x[j];
    x[i] = b[p[i]] - s;
  }
  // backward substitution
  for (int i=n-1; i>=0; --i) {
    double s = 0;
    for (int j=i+1; j<n; ++j) 
      s += a[i][j]*x[j];
    x[i] = (x[i] - s) / a[i][i];
  }
}

void gauss (int n, double** a, double* b, double* x)
{
  for (int j=0; j<n-1; ++j) {
    // find piv�
    int m = j;
    for (int i=j+1; i<n; ++i)
      if (fabs(a[i][j]) > fabs(a[m][j]))
        m = i;
    // swap lines: j <-> m
    double t;
    for (int k=j; k<n; ++k) {
      t = a[j][k];
      a[j][k] = a[m][k];
      a[m][k] = t;
    }
    t = b[j];
    b[j] = b[m];
    b[m] = t;

    // elimination
    for (int i=j+1; i<n; ++i) {
      double f = a[i][j]/a[j][j];
      for (int k=j; k<n; ++k)
        a[i][k] -= f*a[j][k];
      b[i] -= f*b[j];
    }
  }
  // backward substitution
  for (int i=n-1; i>=0; --i) {
    double s = 0;
    for (int j=i+1; j<n; ++j) 
      s += a[i][j]*x[j];
    x[i] = (b[i] - s) / a[i][i];
  }
}

void cholesky (int n, double** A)  // A is symmetric!
{
  for (int k=0; k<n; ++k) {
    A[k][k] = sqrt(A[k][k]);
    for (int i=k+1; i<n; ++i) { 
      A[i][k] /= A[k][k];
    }
    for (int i=k+1; i<n; ++i)
      for (int j=k+1; j<=i; ++j)
        A[i][j] -= A[i][k]*A[j][k];
  }
}

void substituicoes (int n, double** a, double* b, double* x)
{
  // forward substitution
  for (int i=0; i<n; ++i) {
    double s = 0;
    for (int j=0; j<i; ++j) 
      s += a[i][j]*x[j];
    x[i] = (b[i] - s) / a[i][i];
  }
  // backward substitution
  for (int i=n-1; i>=0; --i) {
    double s = 0;
    for (int j=i+1; j<n; ++j) 
      s += a[j][i]*x[j];
    x[i] = (x[i] - s) / a[i][i];
  }
}