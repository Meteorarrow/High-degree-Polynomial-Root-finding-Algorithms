#include <stdio.h>
#include <math.h>
#include <stdbool.h>


#define MAX_DEGREE 100
#define EPSILON 1e-7
#define MAX_ITER 1000


double evaluate_polynomial(double coeffs[], int degree, double x) {
    double result = coeffs[degree];
    for (int i = degree - 1; i >= 0; i--) {
        result = result * x + coeffs[i];
    }
    return result;
}

void derive_polynomial(double poly[], int degree, double deriv[]) {
    for (int i = 1; i <= degree; i++) {
        deriv[i - 1] = poly[i] * i;
    }
}

// Bisection Method
double bisection(double coeffs[], int degree, double a, double b) {
    double fa = evaluate_polynomial(coeffs, degree, a);
    double fb = evaluate_polynomial(coeffs, degree, b);
    
    if (fa * fb >= 0) {
        printf("Bisection failed: No root or even number of roots in [a, b]\n");
        return NAN;
    }

    double mid;
    for (int i = 0; i < MAX_ITER; i++) {
        mid = (a + b) / 2.0;
        double fmid = evaluate_polynomial(coeffs, degree, mid);
        
        if (fabs(fmid) < EPSILON) return mid;
        
        if (fa * fmid < 0) b = mid;
        else {
            a = mid;
            fa = fmid;
        }
    }
    return mid;
}

// Newton's Method
double newton(double poly[], double deriv[], int degree, double x0) {
    double x = x0;
    for (int i = 0; i < MAX_ITER; i++) {
        double fx = evaluate_polynomial(poly, degree, x);
        double dfx = evaluate_polynomial(deriv, degree - 1, x);
        
        if (fabs(dfx) < 1e-10) {
            printf("Newton failed: Derivative too small at x = %f\n", x);
            return NAN;
        }
        
        double x_next = x - fx / dfx;
        if (fabs(x_next - x) < EPSILON) return x_next;
        x = x_next;
    }
    return x;
}

int main() {
    // 範例：f(x) = x^5 - x - 1 (五次多項式)
    double poly[] = {-1, -1, 0, 0, 0, 1};
    int degree = 5;
    double deriv[MAX_DEGREE];
    
    derive_polynomial(poly, degree, deriv);
    
    printf("--- Root Finding for High-Degree Polynomial ---\n");
    
    // 測試二分法
    double root_bis = bisection(poly, degree, 1.0, 2.0);
    printf("Bisection root: %f\n", root_bis);
    
    // 測試牛頓法
    double root_newton = newton(poly, deriv, degree, 1.5);
    printf("Newton root: %f\n", root_newton);
    
    return 0;
}
