//
//  community.cpp
//  lottery
//
//  Created by Will Pearse on 17/03/2012.
//  Copyright 2012 Imperial College London. All rights reserved.
//

#include <iostream>
#include "community.h"
using namespace std;

//////////////////
//HOUSEKEEPING////
//////////////////

static int string_to_int(string var)
{
    int var_i;
    stringstream convert;
    convert << var;
    convert >> var_i;
    return var_i;
}
static vector<string> fill_community(string species, int abundance)
{
    vector<string> output;
    while(abundance > 0)
    {
        output.push_back(species);
        --abundance;
    }
    return output;
}
static double log_add(double x, double y)
{
    return(x + log(y));
}

//////////////////
//CONSTRUCTION////
//////////////////

Community::Community(string species, string abundance, string year, string name)
{
    //Setup
    int abundance_int=string_to_int(abundance), year_int=string_to_int(year);
    vector<string> community=fill_community(species, abundance_int);
    
    //Make the starting values
    communities.push_back(community);
    n_species = 1;
    years.push_back(year_int);
    species_names.push_back(species);
    community_name = name;
}

void Community::add_species(string species, string abundance, string year)
{
    //Setup
    int abundance_int=string_to_int(abundance), year_int=string_to_int(year),i=0;
    vector<string> community=fill_community(species, abundance_int);
    vector<string>::iterator iter;
    
    //Do we already have a community for this year?
    // - can't use find without an iterator, can't use interators as indices --> write own search function
    for(; i<years.size(); ++i)
    {
        if(years[i] == year_int)
        {
            communities[i].insert(communities[i].end(), community.begin(), community.end());
            break;
        }
    }
    
    //If we didn't find a match, we should make a new year
    if(i==years.size())
    {
        communities.push_back(community);
        years.push_back(year_int);
    }
    
    //Insert any new names
    iter = find(species_names.begin(), species_names.end(), species);
    if(iter == species_names.end())
        species_names.push_back(species);
}

Community::Community(int n_years, int total_individuals, int total_additions, std::vector<std::string> sp_names, boost::numeric::ublas::matrix<double> transition_matrix, std::vector<double> addition_rates, std::string name)
{
    //Setup
    assert(total_individuals >= 1);
    vector<string> community_vec;
    int current_vec = 0;
    
    //Make first community from uniform distribution
    boost::mt19937 rnd_generator(NULL);
    vector<double> uniform_weights(sp_names.size());
    for(int i=0; i<sp_names.size(); ++i)
        uniform_weights[i] = 1.0 / sp_names.size();
    boost::random::discrete_distribution<> uniform(uniform_weights.begin(), uniform_weights.end());
    for(int i=0; i<total_individuals; ++i){
        int x = uniform(rnd_generator);
        community_vec.push_back(sp_names[x]);
    }
    
    sort(community_vec.begin(), community_vec.end());
    communities.push_back(community_vec);
    
    //Make generators for all the species
    vector<boost::random::discrete_distribution<> > species_distributions;
    vector<double> current_rates(sp_names.size()+1);
    for(int i=0; i<sp_names.size(); ++i)
    {
        for(int j=0; j<current_rates.size(); ++j)
            current_rates[j] = transition_matrix(i,j);
        boost::random::discrete_distribution<> curr_dist(current_rates.begin(), current_rates.end());
        species_distributions.push_back(curr_dist);
    }
    
    //Make addition generators for all the species
    boost::random::discrete_distribution<> sp_addition_distribution(addition_rates.begin(), addition_rates.end());
    
    //Add in death for the species
    species_names = sp_names;
    n_species = sp_names.size();
    sp_names.push_back("DEATH");
    assert(sp_names.size() == transition_matrix.size2());
    
    //Make real transition matrix
    real_t_m = transition_matrix;
    real_a_r = addition_rates;
    
    //Make next load of communities
    while(++current_vec < n_years)
    {
        //Setup
        boost::numeric::ublas::matrix<int> event_matrix = boost::numeric::ublas::zero_matrix<int>(species_names.size(), species_names.size()+2);
        vector<string> current_com;
        //Go along all species
        for(int i=0; i<communities[current_vec-1].size(); ++i)
        {
            //Find current species
            int j=0;
            for(; j<sp_names.size(); ++j)
                if(sp_names[j]==communities[current_vec-1][i])
                    break;
            
            //Choose next step and record it in the event matrix
            int next_step = species_distributions[j](rnd_generator);
            ++event_matrix(j,next_step);
            //Add the new species in the list, two if it's a reproduction, and nothing if it's death
            if(sp_names[next_step] != "DEATH")
            {
                current_com.push_back(sp_names[next_step]);
            }
        }
        real_e_m.push_back(event_matrix);
        
        //Do additions
        for(int i=0; i<total_additions; ++i)
            current_com.push_back(sp_names[sp_addition_distribution(rnd_generator)]);
        
        //Book-keeping            
        sort(current_com.begin(), current_com.end());
        communities.push_back(current_com);
    }
    //Construct a community from all this data
    // - there's a lot of commonality between the constructors...
    for(int i=0; i<n_years; ++i)
    {
        years.push_back(i);
        communities[i].erase(communities[i].begin());
    }
    community_name = name;
}
