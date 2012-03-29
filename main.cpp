//
//  main.cpp
//  lottery
//   - a simple lottery model of species turnover
//
//  Individuals can continue to the next time-step, turn into something else, or die.
//  Every time-step, species are randomly thrown into the community (the number needed!)
//   - and these probabilities are calculated too.
//  The event matrix is calculated first, then the transition matrix optimised around that (iteratively).
//  *Makes it very difficult to detect death (in my opinion)*
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
        cout << "To run with real data, enter two arguments: the community file, then the number of overall iterations" << endl;
        cout << "To simulate, enter six arguments: the number of communities, number of years in each, number of species, number of starting individuals, the number of additions per year, a stable transition rate, and a random seed." << endl;
        cout << "*Anything* else will either do nothing, or cause a confusing-looking error." << endl;
    }
    
    if(argc == 3)
    {
        //Loading data from file
        cout << "...loading from " << argv[1] << endl;
        Data data(argv[1]);
        
        //Guess parameters
        int runs(atoi(argv[2]));
        data.optimise(0,0,runs);
        
        //Ouput
        data.print_parameters();
	data.print_event_matrix(-1,0);
    }
    
    if(argc == 8)
    {
        //RANDOMISATIONS
        //Setup
        int n_communities = atoi(argv[1]);
        int n_years = atoi(argv[2]);
        int n_additions = atoi(argv[5]);
        cout << n_additions << endl;
        vector<string> sp_names(atoi(argv[3]));
        int n_individuals(atoi(argv[4]));
        double static_freq = atof(argv[6]);
        cout << static_freq << endl;
        int rnd_seed = atoi(argv[7]);
        char letter = 'a';
        for(int i=0; i<sp_names.size(); ++i)
            sp_names[i] = letter++;
        vector<string> community_names(n_communities);
        char big_letter = 'A';
        for(int i=0; i<community_names.size(); ++i)
            community_names[i] = big_letter++;
        double turnover_freq = (1.0 - static_freq)/(sp_names.size()+2);
        boost::numeric::ublas::matrix<double> transition_matrix(sp_names.size(), sp_names.size()+2);
        for(int i=0; i<transition_matrix.size1(); ++i)
            for(int j=0; j<transition_matrix.size2(); ++j)
                if(i==j)
                    transition_matrix(i,j) = static_freq;
                else
                    transition_matrix(i,j) = turnover_freq;
        
        //Randomisations
        Data data(n_communities, n_years, n_individuals, n_additions, sp_names, transition_matrix, community_names, rnd_seed);
        data.print_event_matrix(-1,0,8,0);
        data.print_event_matrix(-1,0,8,1);
        
        //Guess parameters
        cout << endl << "..'optimising'..." << endl << endl;
        data.set_transitions();
        
        data.print_event_matrix(-1,0);
        data.print_parameters();
        
        data.optimise(0,0);
        //Output
        cout << "Log-likelihood: " << data.likelihood() << endl << "Parameters:" << endl;
        data.print_parameters();
        cout << endl << "Estimated events:";
        data.print_event_matrix(-1,0);
    }
    return 0;
}

