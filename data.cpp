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
            communities[i].add_species(species, abundance, year);
        
        //Add species to list if it's new
        for(i=0; i<species_names.size(); ++i)
            if(species == species_names[i])
                break;
        if(i == species_names.size())
            species_names.push_back(species);
    }
    
    n_species = species_names.size();
    n_communities = communities.size();
    for(int i=0; i<communities.size(); ++i)
        n_sites.push_back(communities[i].communities.size());
    
    //Make transition matrix
    transition_matrices.push_back(boost::numeric::ublas::matrix<int>(n_species, n_species+2));
    double null_stationary_freq = 1.0 / (n_species-1);
    double null_background_freq = (1.0 - null_stationary_freq) / (n_species-1);
    double null_add_freq = 1.0 / n_species;
    for(int i=0; i<n_species; ++i)
    {
        int j=0;
        for(; j<(n_species+1); ++j)
            transition_matrices[0](i,j) = null_background_freq;
        transition_matrices[0](i,j) = null_add_freq;
        transition_matrices[0](i,i) = null_stationary_freq;
    }
    
    //Initialise all communities
    for(i=0; i<communities.size(); ++i)
        communities[i].initialise(transition_matrices);
    
}

Data::Data(int no_communities, int n_years, int total_individuals, int total_additions, std::vector<std::string> sp_names, boost::numeric::ublas::matrix<double> transition_matrix, std::vector<std::string> com_names, int rnd_seed)
{
    boost::mt19937 rnd_generator(rnd_seed);

    species_names = sp_names;
    n_communities = no_communities;
    community_names = com_names;
    n_species = sp_names.size();
    
    //Make transition matrix
    transition_matrices.push_back(boost::numeric::ublas::matrix<int>(n_species, n_species+2));
    double null_stationary_freq = 1.0 / (n_species-1);
    double null_background_freq = (1.0 - null_stationary_freq) / (n_species-1);
    double null_add_freq = 1.0 / n_species;
    for(int i=0; i<n_species; ++i)
    {
        int j=0;
        for(; j<(n_species+1); ++j)
            transition_matrices[0](i,j) = null_background_freq;
        transition_matrices[0](i,j) = null_add_freq;
        transition_matrices[0](i,i) = null_stationary_freq;
    }
    
    for(int i=0; i<no_communities; ++i)
    {
        int new_seed = rnd_generator();
        Community temp(n_years, total_individuals, total_additions, sp_names, transition_matrix, com_names[i], new_seed);
        communities.push_back(temp);
    }
}

///////////////
//TRANSITIONS//
///////////////
//Equally, pointer the fuck out of this
void Data::set_transitions(void)
{
    //Go through each community
    for(int i=0; i<communities.size(); ++i)
    {
        communities[i].event_matrices.clear();
        //Go through each community->community transition
        for(int j=0; j<(communities[i].n_years-1); ++j)
        {
            //Pull out the correct t_m
            boost::numeric::ublas::matrix<double> curr_t_m = transition_matrices[communities[i].transition_matrix_index[j]];
            //Use it
            communities[i].event_matrices.push_back(communities[i].set_transitions(curr_t_m, j));
        }
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

//Re-calculate transitions and likelihood with one changed parameter
double inverse_integ_log_lik_trans(double param, vector<Community> communities, boost::numeric::ublas::matrix<double> transition_matrix, int t_m_index, int row, int column)
{
    //Change parameter values in rest of transition_matrix
    double leftover = 1 - param;
    double fudge_factor = 1 - transition_matrix(row,column);
    for(int i=0; i<transition_matrix.size2()-1; ++i)
        transition_matrix(row,i) = (transition_matrix(row,i) / fudge_factor) * leftover;
    transition_matrix(row,column) = param;
	
    vector<double> all_likelihoods;
    //Go through each community
    for(int i=0; i<communities.size(); ++i)
    {
        for(int j=0; j<(communities[i].n_years-1); ++j)
        {
            //Only do this if we're dealing with the transition matrix of interest
            if (j==t_m_index)
            {
                //Re-calculate transitions
                communities[i].set_transitions(transition_matrix, j);
                //Make copy of event_matrix to use as a counter
                boost::numeric::ublas::matrix<int> curr_e_m = communities[i].event_matrices[j];
                //Prepare to hold likelihoods
                vector<double> likelihoods;
                //Acucmulate transition event probabilities
                for(int k=0; k<curr_e_m.size1(); ++k)
                    for(int l=0; l<curr_e_m.size2()-1; ++l)
                        while(curr_e_m(k,l) > 0)
                        {
                            likelihoods.push_back(log(transition_matrix(k,l)));
                            --curr_e_m(k,l);
                        }
                
                all_likelihoods.push_back(accumulate(likelihoods.begin(), likelihoods.end(), 0.0));
            }
        }
        
    }
    return (0.0 - accumulate(all_likelihoods.begin(), all_likelihoods.end(), 0.0));
}

double inverse_integ_log_lik_add(double param, vector<Community> communities, boost::numeric::ublas::matrix<double> transition_matrix, int t_m_index, int sp)
{
    //Change parameter values in rest of transition_matrix
    double leftover = 1 - param;
    double fudge_factor = 1 - transition_matrix(sp,transition_matrix.size2()-1);
    for(int i=0; i<transition_matrix.size1(); ++i)
        transition_matrix(i,transition_matrix.size2()-1) = (transition_matrix(i,transition_matrix.size2()-1) / fudge_factor) * leftover;
    transition_matrix(sp,transition_matrix.size2()-1) = param;
    
    vector<double> all_likelihoods;
    //Go through each community
    for(int i=0; i<communities.size(); ++i)
    {
        for(int j=0; j<(communities[i].n_years-1); ++j)
        {
            //Only do this if we're dealing with the transition matrix of interest
            if (j==t_m_index)
            {
                //Re-calculate transitions
                communities[i].set_transitions(transition_matrix, j);
                //Make copy of event_matrix to use as a counter
                boost::numeric::ublas::matrix<int> curr_e_m = communities[i].event_matrices[j];
                //Prepare to hold likelihoods
                vector<double> likelihoods;
                //Acucmulate addition event probabilities
                for(int k=0; k<curr_e_m.size1(); ++k)
                    while(curr_e_m(k,transition_matrix.size2()-1) > 0)
                    {
                        likelihoods.push_back(log(transition_matrix(k,transition_matrix.size2()-1)));
                        --curr_e_m(k,transition_matrix.size2()-1);
                    }
                all_likelihoods.push_back(accumulate(likelihoods.begin(), likelihoods.end(), 0.0));
            }
        }
        
    }
    return (0.0 - accumulate(all_likelihoods.begin(), all_likelihoods.end(), 0.0));
}

////////////
//OPTIMISE//
////////////

void Data::optimise(int max_communities, int max_years, int mat_iter, int param_iter)
{
    //Setup
    for(int m_iter=0; m_iter<mat_iter; ++m_iter)
    {
        //Randomise the order in which we estimate parameters
        vector<int> order(transition_matrices[0].size1() * transition_matrices[0].size2());
        vector<int> row_order(order.size()), col_order(order.size());
        int row_track=0, col_track=0;
        for(int x=0; x<order.size(); ++x)
        {
            order[x] = x;
            if(col_track < transition_matrices[0].size2())
            {
                col_order[x] = col_track++;
                row_order[x] = row_track;
            }
            else
            {
                col_order[x] = 0;
                col_track = 1;
                ++row_track;
            }
            row_order[x] = row_track;
        }
        vector<int> rnd_order(order);
        random_shuffle(rnd_order.begin(), rnd_order.end());
        //Loop over transition matrices
        for(int i=0; i<transition_matrices.size(); ++i)
        {
            //Randomly loop through each parameter
            for(int j=0; j<rnd_order.size(); ++j)
            {
                //Row and column indices
                int row = row_order[rnd_order[j]];
                int col = col_order[rnd_order[j]];
                
                //Set new parameter
                if(row==0)
                    cout << "Before (" << col << "): " << transition_matrices[i](row, col);
                double fudge_factor = 1.0 - transition_matrices[i](row, col);
                transition_matrices[i](row,col) = boost::math::tools::brent_find_minima(boost::bind(inverse_integ_log_lik_trans, _1, communities, transition_matrices[i], i, row,col), 0.01, 0.99, param_iter).first;
                if(row==0)
                    cout << "," << transition_matrices[i](row,col) << endl;
                //Shift around the others - differently if dealing with an addition
                double leftover = 1.0 - transition_matrices[i](row,col);
                if(col!=n_species+1)
                {
                    for(int l=0; l<(transition_matrices[i].size2()-1); ++l)
                        if(col!=l)
                            transition_matrices[i](row,l) = (transition_matrices[i](row,l) / fudge_factor) * leftover;
                }
                else
                    for(int l=0; l<transition_matrices[i].size1(); ++l)
                        if(l!=row)
                            transition_matrices[i](l,n_species+1) = (transition_matrices[i](l,n_species+1) / fudge_factor) * leftover;
                double total = 0.0;
                for(int x=0; x<transition_matrices[i].size2(); ++x)
                    total += transition_matrices[i](row,x);
                if(total > 1.0){
                    print_parameters();
                    print_event_matrix(-1,0);}
                assert(total <= 1.0);
                
            }
        print_parameters();
        }
    }
}
///////////
//DISPLAY//
///////////
void Data::print_community(int community_index, int year_index, int width)
{
    communities[community_index].print_year(year_index, width);
}

void Data::print_event_matrix(int community_index, int transition_index, int width, int real)
{
    //Are we printing a total event matrix
    if(community_index == -1)
    {
        //Assign the total event matrix
        boost::numeric::ublas::matrix<int> total_events(n_species, n_species+2);
        for(int i=0; i<n_species; ++i)
            for(int j=0; j<(n_species+2); ++j)
                total_events(i,j) = 0;
        if(real)
        {
            for(int i=0; i<n_communities; ++i)
                for(int j=0; j<communities[i].real_e_m.size(); ++j)
                    total_events += communities[i].real_e_m[j];
        }
        else
        {
            for(int i=0; i<n_communities; ++i)
                for(int j=0; j<communities[i].event_matrices.size(); ++j)
                    total_events += communities[i].event_matrices[j];
        }
        //Print it out
        //Header
        cout << endl << setw(width) << "" ;
        for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
            cout << setw(width) << *iter;
        cout << setw(width) << "Death" << setw(width) << "Add." << endl;
        
        //Looping through
        for(int i = 0; i<total_events.size1(); ++i)
        {
            cout << setw(width) << species_names[i];
            for(int j=0; j<total_events.size2(); ++j)
                cout << setw(width) << total_events(i,j);
            cout << endl;
        }
    }
    else
        communities[community_index].print_event_matrix(transition_index, width, real);
}

void Data::print_parameters(int width)
{
    for(int i=0; i<transition_matrices.size(); ++i)
    {
        //Header
        cout << endl << setw(width) << "" ;
        for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
            cout << setw(width) << *iter;
        cout << setw(width) << "Death" << setw(width) << "Add." << endl;
        
        //Looping through
        for(int j=0; j<transition_matrices[i].size1(); ++j)
        {
            cout << setw(width) << species_names[j];
            for(int k=0; k<transition_matrices[i].size2(); ++k)
                if(transition_matrices[i](j,k)>0.0001)
                    cout << setw(width) << setprecision(4) << transition_matrices[i](j,k);
                else
                    cout << setw(width) << setprecision(4) << 0;
            cout << endl;
        }
    }
}
