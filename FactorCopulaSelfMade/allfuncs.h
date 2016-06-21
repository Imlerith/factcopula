//
//  GaussLegendre.h
//  FactorCopulaSelfMade
//
//  Created by Sergey Nasekin on 16/06/16.
//  Copyright (c) 2016 Sergey Nasekin. All rights reserved.
//

#ifndef __FactorCopulaSelfMade__allfuncs__
#define __FactorCopulaSelfMade__allfuncs__

#include <stdio.h>
#include <vector>
using namespace std;

struct Vecs{
    
    vector<double> index;
    vector<double> vals;
    
};

#ifdef __cplusplus
extern "C"
{
#endif
    
    /* Numerical computation of int(f(x),x=a..b) by Gauss-Legendre n-th order high precision quadrature
     [in]n     - quadrature order
     [in]f     - integrand
     [in]data  - pointer on user-defined data which will
     be passed to f every time it called (as second parameter).
     [in][a,b] - interval of integration
     
     return:
     -computed integral value or -1.0 if n order quadrature is not supported
     */
    double gauss_legendre(int n, double (*f)(double,void*), void* data, double a, double b);
    
    /* 2D Numerical computation of int(f(x,y),x=a..b,y=c..d) by Gauss-Legendre n-th order high precision quadrature
     [in]n     - quadrature order
     [in]f     - integrand
     [in]data  - pointer on user-defined data which will
					be passed to f every time it called (as third parameter).
     [in][a,b]x[c,d] - interval of integration
     
     return:
     -computed integral value or -1.0 if n order quadrature is not supported
     */
    double gauss_legendre_2D_cube(int n, double (*f)(double,double,void*), void* data, double a, double b, double c, double d);
    
    /* Computing of abscissas and weights for Gauss-Legendre quadrature for any(reasonable) order n
     [in] n   - order of quadrature
     [in] eps - required precision (must be eps>=macheps(double), usually eps = 1e-10 is ok)
     [out]x   - abscisass, size = (n+1)>>1
     [out]w   - weights, size = (n+1)>>1
     */
    void gauss_legendre_tbl(int n, double* x, double* w, double eps);
    
    
#ifdef __cplusplus
}
#endif

//ADDITIONAL FUNCTIONS IN MY PROJECT

//Test functions for integration
double f(double x, void* data);

double ff(double x);

double fff(double x);

//Trapezoidal numeric integration
//double trapz(double a, double (*func)(double), double b, int n);
double trapz(vector<double> xrange, vector<double> yrange);

//Gauss-Hermite integration
double Gauss_Hermite_Integration_40pts( double (*f)(double,double,double,double,double), double xx, double alpha, double nu, double nuCF );

//Convolution function
double Conv(double t, double x, double alpha, double nu, double nuCF);

//Sorting
Vecs sort_indexes(vector<double> &v);

//Cumulative sum
vector <double> cumsum(vector <double> input);

//Beta function
double beta(double a, double b);

//Student-t pdf
double tpdf(double x, double nu);

//Linspace making
vector<double> linspace(double a, double b, int n);

vector<double> invcdf(double u, double nuCF, double alpha, double nu);


#endif /* defined(__FactorCopulaSelfMade__allfuncs__) */




