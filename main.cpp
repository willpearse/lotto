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
        cout << "To simulate, enter six arguments: the number of communities, number of years in each, number of species, number of starting individuals, the number of additions per year, and a random seed." << endl;
        cout << "*Anything* else will either do nothing, or cause a confusing-looking error." << endl;
        
        vector<string> sp_names(5);
        char letter = 'a';
        for(int i=0; i<5; ++i)
            sp_names[i] = letter++;
        vector<string> community_names(5);
        char big_letter = 'A';
        for(int i=0; i<5; ++i)
            community_names[i] = big_letter++;
        double static_freq = 0.8;
        double turnover_freq = (1.0 - static_freq)/(sp_names.size()+2);
        boost::numeric::ublas::matrix<double> transition_matrix(sp_names.size(), sp_names.size()+2);
        for(int i=0; i<transition_matrix.size1(); ++i)
            for(int j=0; j<transition_matrix.size2(); ++j)
                if(i==j)
                    transition_matrix(i,j) = static_freq;
                else
                    if (j > sp_names.size())
                        transition_matrix(i,j) = 1.0 / sp_names.size();
                    else
                        transition_matrix(i,j) = turnover_freq;
        Data data(5, 10, 100, 5, sp_names, transition_matrix, community_names, 123456);
        data.set_transitions();
        data.optimise(0,0);
        data.print_parameters();
        //data.print_event_matrix(0,0);
        //data.print_event_matrix(0,1);
        //data.print_event_matrix(2,1);
        //data.print_event_matrix(2,5);
    }
    
    if(argc == 3)
    {
        //Loading data from file
        Data data(argv[2]);
        
    }
    
    if(argc == 7)
    {
        //RANDOMISATIONS
        //Setup
        int n_communities = atoi(argv[1]);
        int n_years = atoi(argv[2]);
        int n_additions = atoi(argv[5]);
        vector<string> sp_names(atoi(argv[3]));
        int n_individuals(atoi(argv[4]));
        int rnd_seed = atoi(argv[6]);
        char letter = 'a';
        for(int i=0; i<sp_names.size(); ++i)
            sp_names[i] = letter++;
        vector<string> community_names(n_communities);
        char big_letter = 'A';
        for(int i=0; i<community_names.size(); ++i)
            community_names[i] = big_letter++;
        double static_freq = 0.8;
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
        
        //Guess parameters
        cout << "Log-likelihood: " << data.likelihood() << endl;
        data.set_transitions();
        data.optimise(0,0);
        //Output
        data.print_parameters();
    }
    return 0;
}

