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
