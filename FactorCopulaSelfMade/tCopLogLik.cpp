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


double invcdf(double u, double nuCF, double alpha, double nu)
{
    double C;
    int points = 120;
    vector<double> xrange;
    vector<double> GHintegral;
    vector<double> cumcdf;
    double funcval;

    //Calculate the constant for the pdf and the abscissa space for the integral
    C = 1/( (alpha*sqrt(1-pow(alpha,2)))*sqrt(nuCF*nu)*beta(0.5*nuCF,0.5)*beta(0.5*nu,0.5) );
    xrange = linspace(-15.0,10.0,points);
    
    
    //Perform numerical integration for the pdf
    for(int i = 0; i < xrange.size(); i++){
        
        double loopint = C*Gauss_Hermite_Integration_40pts( Conv, xrange[i], alpha, nu, nuCF );
    
        GHintegral.push_back(loopint);
    
    }
    
    
    //Compute the cumulative pdf (cdf)
    int lag = 30;
    vector<double> xvec;
    vector<double> yvec;
    
    for (int i = 0; i < xrange.size() - lag; i++){

        for (int j = 0; j < lag+i; j++){
        
            xvec.push_back(xrange[j]);
            yvec.push_back(GHintegral[j]);
        
        }
        
        cumcdf.push_back(trapz(xvec,yvec));
        xvec.clear();
        yvec.clear();
    
    }
    
    //Cut out the corresponding vector of abscissas
    vector<double>::const_iterator first = xrange.begin() + lag;
    vector<double>::const_iterator last = xrange.end();
    vector<double> Y(first,last);
    
    
    //Initialize an object from namespace "tk", class "spline"
    tk::spline s;
    s.set_points(cumcdf,Y);
    
    //Find the inverse at u
    funcval = s(u);
    
    //Return the answer
    return funcval;
    
}


//Help function to generate the abscissa space (with equal-sized intervals)
vector<double> linspace(double a, double b, int n) {
    vector<double> array;
    double step = (b-a) / (n-1);
    
    while(a <= b) {
        array.push_back(a);
        a += step;           // could recode to better handle rounding errors
    }
    return array;
}






