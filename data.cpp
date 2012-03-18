//
//  data.cpp
//  lottery
//
//  Created by Will Pearse on 14/03/2012.
//  Copyright 2012 Imperial College London. All rights reserved.
//

#include "data.h"
using namespace std;


///////////////////
//Constructors/////
///////////////////

Data::Data(const char *file)
{
    //Setup
    ifstream str(file);
    string line,cell;
    vector<string> splits;
    string species,abundance,year,community_name;
    int i;
    
    //For each line
    while(getline(str, line))
    {   
        //Get abundance, species, community name and year
        stringstream stream(line);
        getline(stream, abundance, ',');
        getline(stream, species, ',');
        getline(stream, community_name, ',');
        getline(stream, year, ',');
        
        //Do we already have this community?
        for(i=0; i<community_names.size(); ++i)
            if(community_name == community_names[i])
                break;
        if(i == community_names.size())
        {
            //No, so make a new one and store its attributes
            Community temp(species, abundance, year, community_name);
            communities.push_back(temp);
            community_names.push_back(community_name);
        }
        else //Yes, so add to that particular community
            //communities[i].add_species(species, abundance, year);
        
        //Add species to list if it's new
        for(i=0; i<species_names.size(); ++i)
            if(species == species_names[i])
                break;
        if(i == species_names.size())
            species_names.push_back(species);
    }
    
    n_species = species_names.size();
    
    //Initialise all communities
    //for(i=0; i<communities.size(); ++i)
    //    communities[i].initialise();
    
    //Make transition matrix
    transition_matrices.push_back(boost::numeric::ublas::matrix<int>(n_species, n_species));
    double null_freq = 1.0 / n_species;
    for(int i=0; i<n_species; ++i)
        for(int j=0; j<n_species; ++j)
            transition_matrices[0](i,j) = null_freq;
    
}

Data::Data(int no_communities, int n_years, int total_individuals, int total_additions, std::vector<std::string> sp_names, boost::numeric::ublas::matrix<double> transition_matrix, std::vector<std::string> com_names, int rnd_seed)
{
    boost::mt19937 rnd_generator(rnd_seed);
    transition_matrices.push_back(transition_matrix);
    species_names = sp_names;
    n_communities = no_communities;
    community_names = com_names;
    n_species = sp_names.size();
    for(int i=0; i<no_communities; ++i)
    {
        int new_seed = rnd_generator();
        Community temp(n_years, total_individuals, total_additions, sp_names, transition_matrix, com_names[i], new_seed);
        communities.push_back(temp);
    }
}

//////////////
//LIKELIHOOD//
//////////////
//Optimise this function later (pointer the fuck out of it)
double Data::likelihood(void)
{
    vector<double> all_likelihoods;
    //Go through each community
    for(int i=0; i<communities.size(); ++i)
    {
        for(int j=0; j<(communities[i].n_years-1); ++j)
        {
            //Pull out the correct t_m
            boost::numeric::ublas::matrix<double> curr_t_m = transition_matrices[communities[i].transition_matrix_index[j]];
            //Make copy of event_matrix to use as a counter
            boost::numeric::ublas::matrix<int> curr_e_m = communities[i].event_matrices[j];
            //Prepare to hold likelihoods
            vector<double> likelihoods;
            //Acucmulate event probabilities
            for(int k=0; k<curr_e_m.size1(); ++k)
                for(int l=0; l<curr_e_m.size2(); ++l)
                    while(curr_e_m(k,l) > 0)
                    {
                        likelihoods.push_back(log(curr_t_m(k,l)));
                        --curr_e_m(k,l);
                    }
            
            all_likelihoods.push_back(accumulate(likelihoods.begin(), likelihoods.end(), 0.0));
        }
    }
    return accumulate(all_likelihoods.begin(), all_likelihoods.end(), 0.0);
}

///////////
//DISPLAY//
///////////
void Data::print_community(int community_index, int year_index, int width)
{
    communities[community_index].print_year(year_index, width);
}