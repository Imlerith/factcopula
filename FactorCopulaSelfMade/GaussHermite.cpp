//
//  GaussHermite.cpp
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

static const double x[] = {
    1.74537214597582383493e-01,    5.23874713832277192629e-01,
    8.74006612357088077427e-01,    1.22548010904628903095e+00,
    1.57886989493161388625e+00,    1.93479147228229579332e+00,
    2.29391714187508342199e+00,    2.65699599844289579501e+00,
    3.02487988390128443770e+00,    3.39855826585962834626e+00,
    3.77920675343522349311e+00,    4.16825706683250020142e+00,
    4.56750207284439485502e+00,    4.97926097854525587152e+00,
    5.40665424797012760839e+00,    5.85409505603040010826e+00,
    6.32825535122008195560e+00,    6.84023730524935541768e+00,
    7.41158253148546880941e+00,    8.09876113925085005206e+00
};

static const double A[] = {
    3.38643277425589218194e-01,    2.65728251877377076134e-01,
    1.63378732713271457086e-01,    7.84746058654043913071e-02,
    2.93125655361723698462e-02,    8.46088800825813244030e-03,
    1.87149682959795277949e-03,    3.13853594541331475638e-04,
    3.93693398109249277033e-05,    3.63157615069302351182e-06,
    2.41114416367052344169e-07,    1.12123608322758101747e-08,
    3.52562079136541190293e-10,    7.15652805269031870828e-12,
    8.80570764521613225655e-14,    6.00835878949081669008e-16,
    1.98918101211650248554e-18,    2.56759336541166966046e-21,
    8.54405696377551077351e-25,    2.59104371384708147343e-29
};

#define NUM_OF_POSITIVE_ZEROS  sizeof(x) / sizeof(double)
#define NUM_OF_ZEROS           NUM_OF_POSITIVE_ZEROS+NUM_OF_POSITIVE_ZEROS



double Gauss_Hermite_Integration_40pts( double (*f)(double,double,double,double,double), double xx, double alpha, double nu, double nuCF ) {
    
    //All variables are constant apart from the integration variable (t)
    double integral = 0.0;
    const double *pA = &A[NUM_OF_POSITIVE_ZEROS];
    const double *px;
    
    for (px = &x[NUM_OF_POSITIVE_ZEROS - 1]; px >= x; px--)
        integral += *(--pA) * ( (*f)(*px, xx, alpha, nu, nuCF) + (*f)(- *px, xx, alpha, nu, nuCF) );
    
    return integral;
}




