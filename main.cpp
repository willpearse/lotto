//
//  main.cpp
//  lottery
//   - a simple lottery model of species turnover
//
//  Individuals can continue to the next time-step, turn into something else, or die.
//  Every time-step, species are randomly thrown into the community (the number needed!)
//   - and these probabilities are calculated too.
//  The event matrix is calculated first, then the transition matrix optimised around that (iteratively).
//
//  Created by Will Pearse on 14/03/2012.
//  Copyright 2012 Imperial College London. All rights reserved.
//

#include <iostream>
#include "data.h"
using namespace std;

int main (int argc, const char * argv[])
{
    //Welcome header
    cout << "Lottery simulation model" << endl;
    cout << "Will Pearse - 2011 (will.pearse@gmail.com)" << endl;
    
    if(argc == 1)
    {
        cout << "To run with real data, enter two arguments: the community file, then the name of an ouput file" << endl;
        cout << "To simulate, enter four arguments: the number of communities, number of years in each, number of species, and number of individuals." << endl;
        cout << "*Anything* else will either do nothing, or cause a confusing-looking error." << endl;
    }
    
    if(argc == 3)
    {
        Data data(argv[2]);
        
    }
    return 0;
}

