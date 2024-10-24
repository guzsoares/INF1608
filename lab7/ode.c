#include <stdio.h>
#include <math.h>

#define MIN_H 1e-8  
#define MAX_FACTOR 1.2 

double RungeKutta(double t0, double t1, double y0, double (*f)(double t, double y), double h) {
    double t = t0;
    double y = y0;
    while (t < t1) {
        if (t + h > t1) {
            h = t1 - t;
        }
        
        double k1 = h * f(t, y);
        double k2 = h * f(t + h / 2, y + k1 / 2);
        double k3 = h * f(t + h / 2, y + k2 / 2);
        double k4 = h * f(t + h, y + k3);
        
        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        
        t += h;
    }
    return y;
}

double RungeKuttaAdaptativo(double t0, double t1, double y0, double (*f)(double, double), double tol) {
    double h = 1e-7;
    double t = t0;
    double y = y0;

    while (t < t1) {
        if (t + h > t1) {
            h = t1 - t; 
        }

        double k1 = h * f(t, y);
        double k2 = h * f(t + h / 2.0, y + k1 / 2.0);
        double k3 = h * f(t + h / 2.0, y + k2 / 2.0);
        double k4 = h * f(t + h, y + k3);
        double y1 = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

        double h2 = h / 2.0;
        double k1_2 = h2 * f(t, y);
        double k2_2 = h2 * f(t + h2 / 2.0, y + k1_2 / 2.0);
        double k3_2 = h2 * f(t + h2 / 2.0, y + k2_2 / 2.0);
        double k4_2 = h2 * f(t + h2, y + k3_2);
        double y_temp = y + (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2) / 6.0;

        double k1_3 = h2 * f(t + h2, y_temp);
        double k2_3 = h2 * f(t + h2 + h2 / 2.0, y_temp + k1_3 / 2.0);
        double k3_3 = h2 * f(t + h2 + h2 / 2.0, y_temp + k2_3 / 2.0);
        double k4_3 = h2 * f(t + h, y_temp + k3_3);
        double y2 = y_temp + (k1_3 + 2 * k2_3 + 2 * k3_3 + k4_3) / 6.0;

        double delta = fabs(y2 - y1) / 15.0;

        double factor = pow(tol / fabs(delta), 1.0 / 5.0);

        if (delta <= tol) {
            t += h;
            y = y2 + delta;
            
            h *= fmin(MAX_FACTOR, factor);
        } else {
            h *= fmax(0.8 * factor, MIN_H / h);
        }
        if (h < MIN_H) {
            h = MIN_H;
        }
    }

    return y;
}

