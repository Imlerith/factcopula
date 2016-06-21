//
//  tCopLogLik.cpp
//  FactorCopulaSelfMade
//
//  Created by Sergey Nasekin on 16/06/16.
//  Copyright (c) 2016 Sergey Nasekin. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "allfuncs.h"
#include "InterpSpline.h"


double f(double x, void* data)
{
    return sin(x);
}

double ff(double x)
{
    return sin(x);
}

double fff(double x)
{
    //Multiply by the exponential factor for Gauss-Hermite integration
    return exp(x)*sin(x);
}


vector<double> invcdf(double u, double nuCF, double alpha, double nu)
{
    double C;
    vector<double> xrange;
    vector<double> GHintegral;
    vector<double> convcdf;
    vector<double> cumcdf;
//    double funcval;
    Vecs sortedcdf;
    vector<double> Y;
    vector<double> X;
    
    C = 1/( (alpha*sqrt(1-pow(alpha,2)))*sqrt(nuCF*nu)*beta(0.5*nuCF,0.5)*beta(0.5*nu,0.5) );
    xrange = linspace(-15.0,10.0,151);
    
    for(int i = 0; i < xrange.size(); i++){
        
        double loopint = C*Gauss_Hermite_Integration_40pts( Conv, xrange[i], alpha, nu, nuCF );
    
        GHintegral.push_back(loopint);
    
    }
    
//    for(int i = 0; i < xrange.size(); i++){
//        
//        double inttr = trapz(xrange,GHintegral);
//        convcdf.push_back(inttr);
//        
//    }
    
    int lag = 50;
    vector<double> xvec;
    vector<double> yvec;
    
    for (int i = 0; i < xrange.size() - lag; i++){

        for (int j = 0; j < lag+i; j++){
        
            xvec.push_back(xrange[j]);
            yvec.push_back(GHintegral[j]);
        
        }
        
        convcdf.push_back(trapz(xvec,yvec));
        xvec.clear();
        yvec.clear();
    
    }
    
    
//    vector<double>::iterator i;
//    vector<double>::iterator j;
//    for( i = xrange.begin(),j = GHintegral.begin(); i != xrange.end()-49; ++i, ++j ){
//        
//        vector<double> xvec (0,i+49);
//        vector<double> yvec (0,j+49);
//        
//        convcdf.push_back(trapz(xvec,yvec));
//        
//    }
    
    
//    cumcdf = cumsum(convcdf);
//    
//    sortedcdf = sort_indexes(cumcdf);
//    
//    for(int i = 0; i < sortedcdf.index.size(); i++){
//        
//        double sorti = sortedcdf.index[i];
//        
//        Y[i] = xrange[sorti];
//        
//    }
//    
//    X = sortedcdf.vals;
//    
//    tk::spline s;
//    s.set_points(X,Y);
//    
//    
//    funcval = s(u);
    
    
    return convcdf;
    
}

vector<double> linspace(double a, double b, int n) {
    vector<double> array;
    double step = (b-a) / (n-1);
    
    while(a <= b) {
        array.push_back(a);
        a += step;           // could recode to better handle rounding errors
    }
    return array;
}






