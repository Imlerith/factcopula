//
//  ConvFunc.cpp
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


double Conv(double t, double x, double alpha, double nu, double nuCF){

    double funcval;
    
    funcval = exp(pow(t,2))*( ( pow((1+  (pow(t,2))  /((pow(alpha,2))*nuCF)   ),(-0.5*(1+nuCF))  ) )*\
                              ( (pow(   (1 + (pow(x-t,2))/(((1-pow(alpha,2)))*nu)   ),(-0.5*(1+nu))   )) ) );
    
    return funcval;

}