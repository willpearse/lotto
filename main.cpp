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
#include <fstream>
#include <sstream>
#include "data.h"
using namespace std;

static void print_iteration(vector<Data> results, int n_communities, int n_years, int total_individuals, int total_additions, boost::numeric::ublas::matrix<double> transition_matrix, vector<string> species_names)
{
    //Setup
    stringstream ss;
    ss << "results_" << n_communities << "-" << species_names.size() << "-" << n_years << "-" << total_individuals << "-" << total_additions << ".csv";
    string filename(ss.str());
    ofstream file;
    file.open(filename.c_str());
    
    //Header
    file << "iter,row,column,estimated,real" << endl;
    //Loop through
    for(int i=0; i<results.size(); ++i)
    {
        for(int j=0; j<results[i].transition_matrices[0].size1(); ++j)
        {
            for(int k=0; k<results[i].transition_matrices[0].size2(); ++k)
            {
                file << i << "," << j << "," << k << "," << results[i].transition_matrices[0](j,k) << "," << transition_matrix(j,k) << endl;
            }
        }
    }
}


int main (int argc, const char * argv[])
{
    //Welcome header
    cout << "Lottery simulation model" << endl;
    cout << "Will Pearse - 2011 (will.pearse@gmail.com)" << endl;
    
    if(argc == 1)
    {
        cout << "To run with real data, enter two arguments: the community file, then the number of overall iterations" << endl;
        cout << "To run the simulations I present in the paper (the PDF), enter one argument: a random seed" << endl;
        cout << "To do some simulations, enter six arguments: the number of communities, number of years in each, number of species, number of starting individuals, the number of additions per year, a stable transition rate, and a random seed." << endl;
        cout << endl << "*Anything* else will either do nothing, or cause a confusing-looking error." << endl;
    }
    
    if(argc == 2)
    {
        //////////////
        //Setup
        //////////////
        int nComm[] = {5, 10, 20, 50};
        int nYears[] = {5, 15, 30};
        int nIndiv[] = {100, 200};
        int nAdds[] = {5, 20};
        double nFreqs[] = {0.2, 0.5, 0.8};
        int nIter = 30;
        int rnd_seed = atoi(argv[1]);
        boost::mt19937 rnd_generator(rnd_seed);
        vector<Data> iter_results;
        
        //////////////
        //Tests
        //////////////
        //Create species' names
        vector<string> sp_names(5);
        char letter = 'a';
        for(int i=0; i<5; ++i)
            sp_names[i] = letter++;
        
        for(int com=0; com<4; ++com)
        {
            //Create community names
            vector<string> community_names(nComm[com]);
            char LETTER = 'A';
            for(int i=0; i<nComm[com]; ++i)
                community_names[i] = LETTER++;
            
            for(int yr=0; yr<3; ++yr)
                for(int ind=0; ind<2; ++ind)
                    for(int add=0; add<2; ++add)
                        for(int frq=0; frq<3; ++frq)
                        {
                            //Create transition matrix
                            double turnover_freq = (1.0 - nFreqs[frq]) / sp_names.size();
                            double add_freq = 1.0 / sp_names.size();
                            boost::numeric::ublas::matrix<double> transition_matrix(sp_names.size(), sp_names.size()+2);
                            for(int i=0; i<transition_matrix.size1(); ++i)
                            {
                                int j=0;
                                for(j=0; j<(transition_matrix.size2()-1); ++j)
                                    if(i==j)
                                        transition_matrix(i,j) = nFreqs[frq];
                                    else
                                        transition_matrix(i,j) = turnover_freq;
                                transition_matrix(i,j) = add_freq;
                            }
                            
                            for(int iter=0; iter<nIter; ++iter)
                            {
                                //Random number
                                int local_seed = rnd_generator();
                                
                                //Create data...
                                Data data(nComm[com], nYears[yr], nIndiv[ind], nAdds[add], sp_names, transition_matrix, community_names, local_seed);
                                
                                //Optimise
                                data.optimise(0,0);
                                
                                //Store
                                iter_results.push_back(data);
                            }
                            print_iteration(iter_results, nComm[com], nYears[yr], nIndiv[ind], nAdds[add], transition_matrix, sp_names);
                            iter_results.clear();
                        }
            }
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
        
        //Guess parameters
        cout << endl << "..'optimising'..." << endl << endl;
        
        data.optimise(0,0);
        //Output
        cout << "Log-likelihood: " << data.likelihood() << endl << "Parameters:" << endl;
        data.print_parameters();
        cout << endl << "Estimated events:";
        data.print_event_matrix(-1,0);
    }
    return 0;
}

