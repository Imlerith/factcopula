//
//  SortValues.cpp
//  FactorCopulaSelfMade
//
//  Created by Sergey Nasekin on 17/06/16.
//  Copyright (c) 2016 Sergey Nasekin. All rights reserved.
//

#include <stdio.h>
#include <vector>
using namespace std;

//Returns both the sorted vector and the corresponding indices
struct Vecs{

    vector<double> index;
    vector<double> vals;

};


Vecs sort_indexes(vector<double> &v) {
    
    // initialize original index locations
    vector<double> idx(v.size());
    vector<double> vsorted(v.size());
    for (int i = 0; i != idx.size(); ++i) idx[i] = i;
    
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](int i1, int i2) {return v[i1] < v[i2];}); //using C++11 lambdas here!
    
    for (vector<double>::size_type i = 0; i != idx.size(); i++){ //using indices, NOT iterators here!
    
        vsorted[i] = v[idx[i]];
        
    }
    
    Vecs result = {idx,vsorted};
    
    return result;
}