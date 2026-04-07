#include "solver.h"

// Horner's Method: O(d) 的高效求值
double evaluate_horner(double poly[], int d, double x) {
    double result = poly[d];
    for (int i = d - 1; i >= 0; i--) {
        result = result * x + poly[i];
    }
    return result;
}

// 自動求導：將 nx^(n-1) 的邏輯應用到係數陣列
void get_derivative(double poly[], int d, double deriv[]) {
    for (int i = 1; i <= d; i++) {
        deriv[i - 1] = poly[i] * i;
    }
}

// Bisection Method 實作
double bisection_method(double poly[], int d, double a, double b, double eps, int *steps) {
    double fa = evaluate_horner(poly, d, a);
    double fb = evaluate_horner(poly, d, b);
    *steps = 0;

    // 檢查區間端點是否異號 (勘根定理)
    if (fa * fb >= 0) return NAN;

    double mid;
    while ((b - a) > eps && *steps < MAX_ITER) {
        (*steps)++;
        mid = (a + b) / 2.0;
        double fmid = evaluate_horner(poly, d, mid);

        if (fabs(fmid) < 1e-15) break; // 剛好踩到根

        if (fa * fmid < 0) b = mid;
        else {
            a = mid;
            fa = fmid;
        }
    }
    return mid;
}

// Newton's Method 實作
double newton_method(double poly[], int d, double x0, double eps, int *steps) {
    double deriv[MAX_DEGREE];
    get_derivative(poly, d, deriv); // 自動獲取 f'(x) 的係數
    
    double x = x0;
    *steps = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        (*steps)++;
        double fx = evaluate_horner(poly, d, x);
        double dfx = evaluate_horner(deriv, d - 1, x);

        // 防止斜率為 0 (切線水平) 導致除以零錯誤
        if (fabs(dfx) < 1e-12) return NAN;

        double x_next = x - fx / dfx;

        // 判斷是否達到 eps 精細度
        if (fabs(x_next - x) < eps) return x_next;
        x = x_next;
    }
    return x;
}
