//
//  communityData.h
//  transitionMatrix
//
//  Created by Will Pearse on 18/11/2011.
//  Copyright 2011 Imperial College London. All rights reserved.
//
//  Data-holding class
//  Two constructors: randomisations, real data from file
//

#ifndef lottery_data_h
#define lottery_data_h

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
#include "community.h"

class Data{
private:
    int n_communities;
    std::vector<std::string> community_names;
    int n_species;
    std::vector< std::string > species_names;
    std::vector<int> n_sites;
    std::vector<std::vector<int> > total_individuals;
    std::vector<boost::numeric::ublas::matrix<double> > transition_matrices;
    std::vector<Community> communities;
public:
    //Real data from file
    Data(const char *file);
    //Randomisations
    Data(int n_communities, int n_years, int total_individuals[]);
    //Print summary to screen
    void summary(void);
    //Ouput details to file
    void output_file(std::string filename);
    //Optimise transition matrix, with a maximum number of subsets for communities and years
    void optimise(int max_communities, int max_years);
    
    friend class Community;
};

#endif
