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
    std::vector<boost::numeric::ublas::matrix<int> > event_matrices;
    std::vector<std::vector<std::string> > communities;
    std::string community_name;
    
public:
    Community(std::string species, std::string abundance, std::string year, std::string name);
    double calc_likelihood(void);
    void add_species(std::string species, std::string abundance, std::string year);
    
    friend class Data;
};

#endif
