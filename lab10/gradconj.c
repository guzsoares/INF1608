#include "gradconj.h"
#include "matriz.h"
#include "vetor.h"
#include <stdlib.h>

int gradconj(int n, double** A, double* b, double* x, double tol) {
    double* r = vet_cria(n);
    double* d = vet_cria(n);
    double* Ad = vet_cria(n);
    
    mat_multv(n, n, A, x, r);  
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - r[i];     
        d[i] = r[i];             
    }
    
    int k;
    for (k = 0; k < n; k++) {
        if (vet_norma2(n, r) < tol) break;
        
        mat_multv(n, n, A, d, Ad);
        double alpha = vet_escalar(n, r, r) / vet_escalar(n, d, Ad);

        for (int i = 0; i < n; i++) {
            x[i] += alpha * d[i];
            r[i] -= alpha * Ad[i];
        }
        
        double beta = vet_escalar(n, r, r) / (alpha * vet_escalar(n, d, Ad));
        
        for (int i = 0; i < n; i++) {
            d[i] = r[i] + beta * d[i];
        }
    }
    
    vet_libera(r);
    vet_libera(d);
    vet_libera(Ad);
    return k; 
}

int gradconj_jacobi(int n, double** A, double* b, double* x, double tol) {
    double* r = vet_cria(n);
    double* d = vet_cria(n);
    double* Ad = vet_cria(n);
    double* M_inv_r = vet_cria(n);

    double* M_inv = vet_cria(n);
    for (int i = 0; i < n; i++) {
        M_inv[i] = 1.0 / A[i][i]; 
    }
    
    mat_multv(n, n, A, x, r); 
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - r[i]; 
        M_inv_r[i] = M_inv[i] * r[i]; 
        d[i] = M_inv_r[i]; 
    }

    int k;
    double rkz_old = vet_escalar(n, r, M_inv_r); 
    for (k = 0; k < n; k++) {
  
        if (vet_norma2(n, r) < tol) break;


        mat_multv(n, n, A, d, Ad); 
        double alpha = rkz_old / vet_escalar(n, d, Ad);


        for (int i = 0; i < n; i++) {
            x[i] += alpha * d[i];      
            r[i] -= alpha * Ad[i];    
            M_inv_r[i] = M_inv[i] * r[i];  
        }


        double rkz_new = vet_escalar(n, r, M_inv_r); 
        double beta = rkz_new / rkz_old;
        rkz_old = rkz_new; 

        for (int i = 0; i < n; i++) {
            d[i] = M_inv_r[i] + beta * d[i];
        }
    }


    vet_libera(r);
    vet_libera(d);
    vet_libera(Ad);
    vet_libera(M_inv_r);
    vet_libera(M_inv);
    return k;  
}