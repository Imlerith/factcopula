//
//  main.cpp
//  FactorCopulaSelfMade
//
//  Created by Sergey Nasekin on 16/06/16.
//  Copyright (c) 2016 Sergey Nasekin. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <numeric>      // for: std::partial_sum


#include <time.h> 
#include "allfuncs.h"
#include "InterpSpline.h"
#include "fastgl.h"


static double time_consumed = 0; //global variable to count time


using namespace std;

//Define Pi
#ifndef PI
    #define PI 3.1415926535897932384626433832795028841971693993751
#endif


//MAIN FUNCTION

int main () {
    clock_t start, end;
    
    double DataArray[255][27];
    double ParamArray[27][2];
    double margpar = 5.49;
    
    
    ifstream datafile ("/Users/nasekins/Desktop/ufile.txt");
    ifstream paramfile ("/Users/nasekins/Desktop/parfile.txt");
//    ofstream newfile ("/Users/nasekins/Desktop/ofile.txt");

    
    //IMPORT DATA INTO THE DATA ARRAY
    for(int i = 0; i < 255; i++){
        for(int j = 0; j < 27; j++){
            datafile >> DataArray[i][j];
            /*
             newfile << myArray[i][j];
             newfile << ";";
             */
            cout << "Data entry is " << DataArray[i][j] << endl;
        }
        /*
         newfile << "\n";
         */
    }
    
    //IMPORT DATA INTO THE PARAMETER ARRAY
    for(int i = 0; i < 27; i++){
        for(int j = 0; j < 2; j++){
            paramfile >> ParamArray[i][j];
            cout << "Parameter entry is " << ParamArray[i][j] << endl;
        }
    }
    

//    // Numerical approximation of an integral
//    double approx1;
//    double approx2;
//    double approx3;
//    
//    //Order of integration
//    int n = 20;
//    
//    //Limits of integration
//    double a = 0;
//    double b = 1;
//    double xx = -9.3418;
//    double alpha = 0.5;
//    double nu = 4.0;
//    double nuCF = 6.0;
    
//    vector<double> X(5), Y(5);
//    X[0]=0.1;
//    X[1]=0.4;
//    X[2]=1.2;
//    X[3]=1.8;
//    X[4]=2.0;
//    Y[0]=0.1;
//    Y[1]=0.7;
//    Y[2]=0.6;
//    Y[3]=1.1;
//    Y[4]=0.9;
//
//    
//    
//    //Compute the integral
//    
//    approx1 = gauss_legendre(n,f,NULL,a,b);
////    approx2 = trapz(a, ff, b, n);
//    approx2 = trapz(X,Y);
//    approx3 = Gauss_Hermite_Integration_40pts( Conv, xx, alpha, nu, nuCF );
//    
//    printf("n = %4d: GL value = %.15g\n",n,approx1); //this is C-style output! 
//    printf("n = %4d: TR value = %.15g\n",n,approx2);
//    printf("n = %4d: GH value = %.15g\n",40,approx3);
    
//
//    tk::spline s;
//    s.set_points(X,Y);
//    
//    double x=1.5;
//    
//    printf("spline at %f is %f\n", x, s(x));
//    
//    //Compute the convolution integral
//    vector<double> vec(10);
//    
//    vec[0]=0.1;
//    vec[1]=0.4;
//    vec[2]=1.2;
//    vec[3]=1.8;
//    vec[4]=2.0;
//    vec[5]=0.3;
//    vec[6]=0.7;
//    vec[7]=0.6;
//    vec[8]=1.1;
//    vec[9]=0.9;
//    
//    Vecs output;
//    vector<double> cumulsum;
//    
//    //Sorted vector and indices in the original array
//    output = sort_indexes(vec);
//    
//    //Computing the cumulative sum
//    cumulsum = cumsum(output.vals);
//    
//    //Printing out the corresponding indices
//    for (vector<double>::const_iterator i = output.index.begin(); i != output.index.end(); ++i)
//        cout << "The corresponding index is: " << *i << "\n" << endl;
//    
//    //Printing out the values of the sorted array
//    for (vector<double>::const_iterator i = output.vals.begin(); i != output.vals.end(); ++i){
//        
//        cout << "The sorted vector is: " << *i << "\n" << endl;
//        
//    }
//    
//    //Printing out cumulative sum
//    for(int i = 0; i < cumulsum.size(); i++){
//        
//        cout << "cumsum: " << cumulsum[i] << "\n";
////        newfile << cumulsum[i];
////        newfile << "\n";
//    }
//    
//    cout << "Beta function: " << beta(nu,nuCF) << "\n";
//    cout << "tpdf: " << beta(0.5,nuCF) << "\n";
//    
//    vector<double> linsp;
//    linsp = linspace(-10, 10, 101);
//    
//    for(int i = 0; i < linsp.size(); i++){
//        
//        cout << "Linspace: " << linsp[i] << "\n";
//    }

    
//    double inv;
//    inv = invcdf(0.1, nuCF, alpha, nu);
//    cout << "InvCdf: " << inv << "\n";
    
//    double pcop;
//    pcop = paircop(0.1, 0.2, alpha, nu, nuCF);
//    cout << "Pair copula: " << pcop << "\n";
    
    
    double fcopll;
    start = clock();
    fcopll = tCopLogLik(DataArray, ParamArray, margpar);
    cout << "Log-likelihood: " << fcopll << "\n";
    end = clock();
    
//
//
//    //double val = invcdf(0.1, nuCF, alpha, nu);
    

    time_consumed += (double)(end - start) / CLOCKS_PER_SEC;
    
    cout << "Time consumed is: " <<time_consumed << endl;
    
    return 0;
}








