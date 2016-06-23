//
//  fastgl.h
//  FactorCopulaSelfMade
//
//  Created by Sergey Nasekin on 22/06/16.
//  Copyright (c) 2016 Sergey Nasekin. All rights reserved.
//

#ifndef __FactorCopulaSelfMade__fastgl__
#define __FactorCopulaSelfMade__fastgl__

#include <stdio.h>
#include <stddef.h>
#include <cmath>

// Functions for fastgl in double precision
namespace fastgl {
    // A struct for containing a Node-Weight pair
    struct QuadPair {
        double theta, weight;
        
        // A function for getting the node in x-space
        double x() {return cos(theta);}
        
        // A constructor
        QuadPair(double t, double w) : theta(t), weight(w) {}
        QuadPair() {}
    };
    
    // Function for getting Gauss-Legendre nodes & weights
    // Theta values of the zeros are in [0,pi], and monotonically increasing.
    // The index of the zero k should always be in [1,n].
    // Compute a node-weight pair:
    QuadPair GLPair(size_t n, size_t k);
}

#endif /* defined(__FactorCopulaSelfMade__fastgl__) */
