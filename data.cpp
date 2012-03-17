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

Data::Data(int n_communities, int n_years, int total_individuals, int total_additions, std::vector<std::string> sp_names, boost::numeric::ublas::matrix<double> transition_matrix, std::vector<double> addition_rates, std::vector<std::string> community_names)
{
    for(int i=0; i<n_communities; ++i)
    {
        
        Community temp(n_years, total_individuals, total_additions, sp_names, transition_matrix, addition_rates, community_names[i]);
        communities.push_back(temp);
    }
}