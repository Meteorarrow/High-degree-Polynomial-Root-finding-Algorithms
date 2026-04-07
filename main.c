#include "solver.h"

int main() {
    int d, steps_b, steps_n;
    double poly[MAX_DEGREE];
    double eps = 1e-7; // 設定精細度 eps

    printf("=== High-Degree Polynomial Root Finder ===\n");
    printf("Enter the degree of polynomial (d): ");
    scanf("%d", &d);

    printf("Enter coefficients from a0 to ad (e.g., for x^2 - 2, enter -2 0 1):\n");
    for (int i = 0; i <= d; i++) {
        printf("a%d: ", i);
        scanf("%lf", &poly[i]);
    }

    // --- 實驗 1: Bisection Method ---
    double a, b;
    printf("\n[Bisection] Enter search interval [a, b]: ");
    scanf("%lf %lf", &a, &b);
    double res_b = bisection_method(poly, d, a, b, eps, &steps_b);

    // --- 實驗 2: Newton's Method ---
    double x0;
    printf("[Newton] Enter initial guess x0: ");
    scanf("%lf", &x0);
    double res_n = newton_method(poly, d, x0, eps, &steps_n);

    // --- 結果對比分析 ---
    printf("\n-------------------------------------------\n");
    printf("Results with eps = %.1e:\n", eps);
    
    if (!isnan(res_b)) 
        printf("Bisection Root: %f (Steps: %d)\n", res_b, steps_b);
    else 
        printf("Bisection failed (Check your interval).\n");

    if (!isnan(res_n)) 
        printf("Newton Root:    %f (Steps: %d)\n", res_n, steps_n);
    else 
        printf("Newton failed (Check your initial guess).\n");
    
    printf("-------------------------------------------\n");
    printf("Efficiency Gain: Newton was %.1fx fewer steps than Bisection.\n", 
           (double)steps_b / steps_n);

    return 0;
}
