//
//  TrapzInt.cpp
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
#include "fastgl.h"

double trapz(vector<double> xrange, vector<double> yrange)
{
    long n = yrange.size();
    double a = xrange[0];
    double b = xrange[1];
    
    //Evaluate delta
    double h = (b-a);
    
    //Evaluate endpoints
    double value = 0.5*h*( yrange[0] + yrange[n-1]);
    
    //Evaluate midpoints
    for(int k = 1; k < n-1; k++){
        value += h*yrange[k];
    }

    return value;
}


//double trapz(double a, double (*func)(double), double b, int n)
//{
//    double h = (b-a)/(n-1);
//    //  Evaluate endpoints
//    double value = 0.5*( (*func)(a)+ (*func)(b));
//    //  Now the midpoints
//    for(int k=2; k < n; k++){
//        value += (*func)( a + h*(k-1) );
//    }
//    value*=h;
//    return value;
//}