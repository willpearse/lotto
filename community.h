//
//  community.h
//  lottery
//
//  Created by Will Pearse on 14/03/2012.
//  Copyright 2012 Imperial College London. All rights reserved.
//

#ifndef lottery_community_h
#define lottery_community_h

#include <functional>
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

class Community{
private:
    int n_years, n_species;
    std::vector<int> years;
    std::vector<std::string> species_names;
    std::vector<int> n_individuals;
    std::vector<int> transition_matrix_index;
    boost::numeric::ublas::matrix<int> real_t_m;
    std::vector<boost::numeric::ublas::matrix<int> > real_e_m;
    std::vector<std::vector<std::string> > communities;
    std::string community_name;

public:
    Community(std::string species, std::string abundance, std::string year, std::string name);
    Community(int n_years, int total_individuals, int total_additions, std::vector<std::string> sp_names, boost::numeric::ublas::matrix<double> transition_matrix, std::string name, int rnd_seed);
    void initialise(std::vector<boost::numeric::ublas::matrix<double> > transition_matrices);
    std::vector<boost::numeric::ublas::matrix<int> > event_matrices;
    double calc_likelihood(void);
    void add_species(std::string species, std::string abundance, std::string year);
    void print_year(int index, int width);
    void print_event_matrix(int transition_index, int width, int real=0);
    //Find likely transitions
    boost::numeric::ublas::matrix<int> set_transitions(boost::numeric::ublas::matrix<double> transition_matrix, int community_transition);
    
    friend class Data;
    friend double inverse_integ_log_lik_trans(double param, std::vector<Community> communities, boost::numeric::ublas::matrix<double> transition_matrix, int t_m_index, int row, int column);
    friend double inverse_integ_log_lik_add(double param, std::vector<Community> communities, boost::numeric::ublas::matrix<double> transition_matrix, int t_m_index, int sp);
};

#endif
