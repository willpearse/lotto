//
//  communityData.h
//  transitionMatrix
//
//  Created by Will Pearse on 18/11/2011.
//  Copyright 2011 Imperial College London. All rights reserved.
//
//  Data-holding class
//  Two constructors: randomisations, real data from file
//  NOTE: 'Reproduction' and 'death' are never kept in species_name variables - they must be added each time for forwards-compatibility
//

#ifndef lottery_data_h
#define lottery_data_h

#include <functional>
#include <algorithm>
#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <numeric>
#include <boost/bind.hpp>
#include <math.h>
#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/tools/minima.hpp>
#include "community.h"

class Data{
private:
    int n_communities;
    std::vector<std::string> community_names;
    int n_species;
    std::vector< std::string > species_names;
    std::vector<int> n_sites;
    std::vector<std::vector<int> > total_individuals;
    std::vector<Community> communities;
    boost::numeric::ublas::matrix<double> real_transition_matrix;

public:
    std::vector<boost::numeric::ublas::matrix<double> > transition_matrices;
    //Real data from file
    Data(const char *file);
    //Randomisations
    Data(int n_communities, int n_years, int total_individuals, int total_additions, std::vector<std::string> sp_names, boost::numeric::ublas::matrix<double> transition_matrix, std::vector<std::string> community_names, int rnd_seed);
    //Print summary to screen
    void summary(void);
    //Ouput details to file
    void output_file(std::string filename);
    //Get transitions for each community
    void set_transitions(void);
    //Calculate likelihood of all communities
    double likelihood(void);
    //Optimise transition matrix, with a maximum number of subsets for communities and years
    void optimise(int max_communities, int max_years, int mat_iter=5, int param_iter=100);
    //Print a community out
    void print_community(int community_index, int year_index, int width=8);
    //Print out parameters of all matrices
    void print_parameters(int width=8, int real=0);
    //Print out event matrix for a community's years
    void print_event_matrix(int community_index, int transition_index, int width=8, int real=0);
};

double inverse_integ_log_lik_trans(double param, std::vector<Community> communities, boost::numeric::ublas::matrix<double> transition_matrix, int t_m_index, int row, int column);
double inverse_integ_log_lik_add(double param, std::vector<Community> communities, boost::numeric::ublas::matrix<double> transition_matrix, int t_m_index, int sp);

#endif
