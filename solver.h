#ifndef SOLVER_H
#define SOLVER_H

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

// 定義最大多項式次數與最大疊代次數
#define MAX_DEGREE 100
#define MAX_ITER 1000

// --- 函式宣告 ---

// 使用 Horner's Method 計算 f(x)
double evaluate_horner(double poly[], int d, double x);

// 自動計算多項式的導數係數
void get_derivative(double poly[], int d, double deriv[]);

// 二分搜尋法 (Bisection Method)
// [a, b] 為搜尋區間，eps 為精細度
double bisection_method(double poly[], int d, double a, double b, double eps, int *steps);

// 牛頓疊代法 (Newton's Method)
// x0 為初始猜測值，eps 為精細度
double newton_method(double poly[], int d, double x0, double eps, int *steps);

#endif
