#include <stdio.h>
#include <math.h>

#define MAX_DEGREE 100
#define MAX_ITER 1000

// Bolzano's Theorem
double evaluate(double poly[],int d,double x) {
    double result=poly[0];
    for(int i=1;i<=d;i++) {
        result=result*x+poly[i];
    }
    return result;
}

// Bisection Method
double bisection(double poly[],int d,double a,double b,double eps,int *steps) {
    double fa=evaluate(poly,d,a);
    double fb=evaluate(poly,d,b);
    *steps=0;

    if (fa*fb>=0) return NAN; 

    double mid;
    while ((b-a)>eps&&*steps<MAX_ITER) {
        (*steps)++;
        mid=(a+b)/2.0;
        double fmid=evaluate(poly,d,mid);
        if(fabs(fmid)<1e-15) break;

        if(fa*fmid<0) b=mid;
        else{a=mid;fa=fmid;}
    }
    return mid;
}

// Newton's Method
double newton(double poly[],int d,double x0,double eps,int*steps) {
    double deriv[MAX_DEGREE];
    for (int i=0;i<d;i++) {
        deriv[i]=poly[i]*(d-i);
    }
    
    double x=x0;
    *steps=0;
    for (int i=0;i<MAX_ITER;i++) {
        (*steps)++;
        double fx=evaluate(poly,d,x);
        double dfx=evaluate(deriv,d-1,x); 

        if(fabs(dfx)<1e-12) return NAN;

        double x_next=x-fx/dfx;
        if(fabs(x_next-x)<eps) return x_next;
        x=x_next;
    }
    return x;
}

int main() {
    int d,s_b,s_n;
    double poly[MAX_DEGREE+1];
    double eps=1e-7;

    printf("=== Polynomial Root Finder ===\n");
	printf("Please set your function\n");
    printf("What is the degree of function? ");
    scanf("%d",&d);

    printf("Please enter the coefficients in descending order:\n");
    for (int i=0;i<=d;i++) {
        printf("x^%d: ",d-i);
        scanf("%lf", &poly[i]);
    }
    printf("\nYour function is: f(x) = ");
    bool first=true;
    for (int i=0;i<=d;i++) {
        int power=d-i;
        if (poly[i]==0) continue;

        if (!first&&poly[i]>0) printf(" + ");
        if (poly[i]<0) printf(" - ");
        
        double abs_val=fabs(poly[i]);
        if(abs_val!=1.0||power==0) printf("%.2f",abs_val);
        
        if(power>0) printf("x");
        if(power>1) printf("^%d",power);
        first=false;
    }
    printf("\n");

    double a,b,x0;
    printf("\n[Bisection] Enter search interval [a, b]: ");
    scanf("%lf %lf",&a,&b);
    printf("[Newton] Enter initial guess x0: ");
    scanf("%lf",&x0);

    double res_b=bisection(poly,d,a,b,eps,&s_b);
    double res_n=newton(poly,d,x0,eps,&s_n);

    printf("\n-------------------------------------------");
    printf("\nAnalysis Results (eps = %.1e):", eps);
    if (!isnan(res_b)) printf("\nBisection Root: %.7f (Steps: %d)",res_b,s_b);
    else printf("\nBisection: Failed (Check Interval)");

    if (!isnan(res_n)) printf("\nNewton Root: %.7f (Steps: %d)",res_n,s_n);
    else printf("\nNewton: Failed (Check Initial Guess)");

    printf("\n");

    return 0;
}
