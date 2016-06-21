//
//  CumSum.cpp
//  FactorCopulaSelfMade
//
//  Created by Sergey Nasekin on 17/06/16.
//  Copyright (c) 2016 Sergey Nasekin. All rights reserved.
//

#include <stdio.h>
#include <vector>
using namespace std;


vector <double> cumsum(vector <double> input){
    
    vector<double> cumul;

    double sum = 0;
    
    for (vector<double>::const_iterator i = input.begin(); i != input.end(); ++i){
        
        sum += *i;
        cumul.push_back(sum);
        
    }
    
    return cumul;

}

