//
//  TDist_Gamma_Beta.cpp
//  FactorCopulaSelfMade
//
//  Created by Sergey Nasekin on 17/06/16.
//  Copyright (c) 2016 Sergey Nasekin. All rights reserved.
//

#include <stdio.h>
#include <cmath>

//Define Pi
#ifndef PI
    #define PI 3.1415926535897932384626433832795028841971693993751
#endif

double tpdf(double x, double nu){

    double pdfval = ( (tgamma((nu+1)/2))/(tgamma(nu/2)) )*(1/sqrt(nu*PI))*( 1/( pow((1+((pow(x,2))/nu)),((nu+1)/2)) ) );
    
    return pdfval;

}

double beta(double a, double b){
    
    double betafunc = ( tgamma(a)*tgamma(b) )/tgamma(a+b);
    
    return betafunc;
    
}