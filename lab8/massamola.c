#include <stdio.h>
#include <math.h>

void forca(double t, double x, double y, double* fx, double* fy) {
    double m = 2.0;
    double k = 0.3; 
    double r = 200; 
    double g = 9.8; 
    double w0 = 10.0; 

    *fx = 0;
    *fy = m * g;

    *fx += w0 * exp(-t / 20);
    *fy += 0;

    double dist = sqrt(x * x + y * y);
    double fx_mola = -k * (dist - r) * (x / dist);
    double fy_mola = -k * (dist - r) * (y / dist);
    
    *fx += fx_mola;
    *fy += fy_mola;
}


double evolui(double t, double h, double x, double y, double xp, double yp, double* xn, double* yn) {
    double fx, fy;
    forca(t, x, y, &fx, &fy);

    double delta = 0.002;
    double m = 2.0;

    *xn = x + (1 - delta) * (x - xp) + (h * h * fx / m);
    *yn = y + (1 - delta) * (y - yp) + (h * h * fy / m);

    return t + h;
}

void simula(double x0, double y0, double t_final, int n, double* x, double* y) {
    double h = t_final / n;
    double xp = x0, yp = y0;
    x[0] = x0;
    y[0] = y0;

    double t = 0;
    for (int i = 1; i < n; i++) {
        t = evolui(t, h, x[i - 1], y[i - 1], xp, yp, &x[i], &y[i]);
        xp = x[i - 1];
        yp = y[i - 1];
    }
}
