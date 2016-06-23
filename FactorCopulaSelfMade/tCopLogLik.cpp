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
#include "fastgl.h"



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
    try{
        funcval = s(u);
    }
    catch(...){
        funcval = numeric_limits<double>::quiet_NaN();
    }
    
    
    //Return the answer
    return funcval;
    
}

//PAIR COPULA ESTIMATION
double paircop(double v, double u, double alpha, double nu, double margpar){

    double h = 0.001;
    double funcl;
    double funcr;
    double approxder;
    double funcval;
    
    //Compute approximate derivatives
    funcl = invcdf(u+h,margpar,alpha,nu);
    funcr = invcdf(u-h,margpar,alpha,nu);
    approxder = (funcl-funcr)/(2*h);
    
    //Compute the pair copula
    funcval = tpdf(( invcdf(u,margpar,alpha,nu) - alpha*v )/( sqrt(1 - pow(alpha,2)) ),nu)*(1/sqrt(1-pow(alpha,2)))*approxder;
    
    return funcval;

}


double tCopLogLik(double udata[255][27], double param[27][2], double margpar){
    
    int nodes = 10;
    double loglik = 0.0;

    for(int i = 255; i--; ){
        
        //Initializing the integral (G.-L.)
        double Intgl = 0.0;
        for(int k = 1 ; k <= nodes ; ++k){
            
            fastgl::QuadPair p = fastgl::GLPair(nodes, k); //initializing G.-L. the nodes-weights object
            //Initialize the product inside the integral
            
            double copprod = paircop( 0.5*(p.x()+1.0), udata[i][0], param[0][0], param[0][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][1], param[1][0], param[1][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][2], param[2][0], param[2][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][3], param[3][0], param[3][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][4], param[4][0], param[4][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][5], param[5][0], param[5][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][6], param[6][0], param[6][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][7], param[7][0], param[7][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][8], param[8][0], param[8][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][9], param[9][0], param[9][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][10], param[10][0], param[10][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][11], param[11][0], param[11][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][12], param[12][0], param[12][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][13], param[13][0], param[13][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][14], param[14][0], param[14][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][15], param[15][0], param[15][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][16], param[16][0], param[16][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][17], param[17][0], param[17][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][18], param[18][0], param[18][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][19], param[19][0], param[19][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][20], param[20][0], param[20][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][21], param[21][0], param[21][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][22], param[22][0], param[22][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][23], param[23][0], param[23][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][24], param[24][0], param[24][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][25], param[25][0], param[25][1], margpar)*\
            paircop( 0.5*(p.x()+1.0), udata[i][26], param[26][0], param[26][1], margpar);
            
            
//            double copprod = 1.0;
//            for(int j = 0; j < 27; j++ ){
//                copprod *= paircop( 0.5*(p.x()+1.0), udata[i][j], param[j][0], param[j][1], margpar);
//            }//end of the product sum
            
            Intgl += 0.5*p.weight*copprod;
            

        }//end of the quadrature sum

        if (Intgl != 0.0){
            loglik += log(Intgl);
        }
        

    }//end of the sample sum
    
    return -loglik;

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






