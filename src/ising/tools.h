/* This module defines a number of auxilary tools for the generalized Ising solver.
 *
 * Author: Wenxuan Huang
 * Maintainer: Wenxuan Huang, Daniil Kitchaev
 * Date: 15 March, 2015
 *
 * Copyright: Wenxuan Huang (C) 2014, All rights reserved.
 */

#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <cerrno>
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"
#include <boost/lexical_cast.hpp>
#include "boost/regex.hpp"
#include "boost/algorithm/string/regex.hpp"

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cmath>
#include "math.h"
//#include "tools.cpp"



using namespace std;
using namespace ::boost::tuples;
using namespace ::boost;

void matchnumber(map<tuple<int,int,int,int,int,int>, double >& dictionary,double number,set<tuple<int,int,int,int,int,int> > &list);

inline double round(double d)
{
    return floor(d + 0.5);
}

long long myPow(long long x, long long p);

inline int positive_modulo(int i, int n)
{
    return (i % n + n) % n;
}

template <typename Iterator>
inline bool next_combination(const Iterator first, Iterator k, const Iterator last)
{
    /* Credits: Thomas Draper */
    if ((first == last) || (first == k) || (last == k))
        return false;
    Iterator itr1 = first;
    Iterator itr2 = last;
    ++itr1;
    if (last == itr1)
        return false;
    itr1 = last;
    --itr1;
    itr1 = k;
    --itr2;
    while (first != itr1)
    {
        if (*--itr1 < *itr2)
        {
            Iterator j = k;
            while (!(*itr1 < *j)) ++j;
            std::iter_swap(itr1,j);
            ++itr1;
            ++j;
            itr2 = k;
            std::rotate(itr1,j,last);
            while (last != j)
            {
                ++j;
                ++itr2;
            }
            std::rotate(k,itr2,last);
            return true;
        }
    }
    std::rotate(first,k,last);
    return false;
}

template <typename type1>
void convertJ_neg1_pos1toJ(map<set<type1>,double> &J_negative1_positive1, map<set<type1>,double> &J,double &constant_term)
{
    for (typename map<set<type1>,double>::iterator it1=J_negative1_positive1.begin(); it1!=J_negative1_positive1.end(); it1++) {

        set<type1> tuple_set=it1->first;
        vector<type1> tuple_vector;

        for (typename set<type1>::iterator it=tuple_set.begin();it!=tuple_set.end();it++)
        {
            tuple_vector.push_back(*it);
        }

        double value=it1->second;
        int size_of_set=tuple_set.size();
//        cout<<"\n what is the tuple set";
//        printvector(tuple_set);

//        cout<<"\n start enumeration";
        for (size_t k=size_of_set; k>0.1; k--) {

//            cout<<"\n now k is "<<k;
            int k_int=k;
            do
            {
                set<type1> subset;
                for(typename vector<type1>::iterator it2=tuple_vector.begin(); it2!=tuple_vector.begin()+k;it2++)
                {
                    subset.insert(*it2);
                }
//                cout<<"\n what is the subset generated";
//                printvector(subset);
                J[subset]+=value*myPow(2,k_int)*myPow(-1,size_of_set-k_int);
            }
            while(next_combination(tuple_vector.begin(),tuple_vector.begin() + k,tuple_vector.end()));
        }

        constant_term+=value*myPow(-1,size_of_set);

    }
}

template<typename type1>
void findminmax(map<type1, double > &dictionary,double &min,double &max )
{
    min=1e100;
    max=-1e100;
    for (typename map<type1, double >::iterator ait=dictionary.begin(); ait!=dictionary.end(); ait++)
    {

        if (ait->second<min)
        {
            min=ait->second;
        }

        if (ait->second>max)
            max=ait->second;

    }
}

template<typename type1,typename type2>
void printmap(map<type1,type2>&thismap)
{
    cout<<"\n";
    for(typename map<type1,type2>::iterator it = thismap.begin(); it != thismap.end(); ++it)
    {
        cout<<(*it).first<<":"<<(*it).second<<endl;
        //do something
    }
}

template<typename type1,typename type2,typename type3>
void printmapofmap(map<type1,map<type2,type3> >&thismap)
{
    cout<<"\n";
    for(typename map<type1,map<type2,type3> >::iterator it = thismap.begin(); it != thismap.end(); ++it)
    {
        cout<<(*it).first<<":"<<endl;
        printmap((*it).second);
        //do something
    }
}

template<typename type1>
void printvector(type1 &thisvector)
{
    //    cout<<"\n";
    for(typename type1::iterator it = thisvector.begin(); it != thisvector.end(); ++it)
    {
        cout<<(*it)<<", ";
        //do something
    }
}

template<typename type1>
void printvectorofvector(type1 &thisvector)
{
    //    cout<<"\n";
    for(typename type1::iterator it = thisvector.begin(); it != thisvector.end(); ++it)
    {
        cout<<"\n";
        printvector(*it);
        //do something
    }
}

template<typename type1>
void printvectorofmap(type1 &thisvector)
{
    //    cout<<"\n";
    for(typename type1::iterator it = thisvector.begin(); it != thisvector.end(); ++it)
    {
        cout<<"\n ";
        printmap(*it);
        //do something
    }
}

template<typename type1,typename vectortype>
void printmapofvector(map<type1,vectortype>&thismap)
{
    cout<<"\n";
    for(typename map<type1,vectortype >::iterator it = thismap.begin(); it != thismap.end(); ++it)
    {

        cout<<(*it).first<<": ";
        printvector((*it).second);
        cout<<"\n";
    }

}

template<typename type1,typename type2>
void printmapfromsets(map<type1,type2>&thismap)
{
    cout<<"\n";
    for(typename map<type1,type2>::iterator it = thismap.begin(); it != thismap.end(); ++it)
    {

        printvector((*it).first);
        cout<<": ";
        cout<<(*it).second;
        cout<<"\n";
    }

}

template<typename type1,typename type2,typename type3>
void getkeystoset(map<type1,type2> &thismap, set<type1> &thisset)
{
    for(typename map<type1,type2>::iterator it = thismap.begin(); it != thismap.end(); ++it)
    {
        thisset.insert(*(it).first);
    }
}

template<typename type1,typename type2,typename type3>
void getkeystovector(map<type1,type2> &thismap, vector<type1> &thisset)
{
    for(typename map<type1,type2>::iterator it = thismap.begin(); it != thismap.end(); ++it)
    {
        thisset.push_back(*(it).first);
    }
}

template<int k>
int getimax(map<tuple<int,int>,int> & spin)
{
    vector<int> thisvector;
    for(typename map<tuple<int,int>,int>::iterator it = spin.begin(); it != spin.end(); ++it)
    {
        thisvector.push_back((*it).first.get<k>());
    }

    vector<int>::iterator it=max_element(thisvector.begin(), thisvector.end());

    return *it;

}

template<int k>
int getimin(map<tuple<int,int>,int> & spin)
{
    vector<int> thisvector;
    for(typename map<tuple<int,int>,int>::iterator it = spin.begin(); it != spin.end(); ++it)
    {
        thisvector.push_back((*it).first.get<k>());
    }

    vector<int>::iterator it=min_element(thisvector.begin(), thisvector.end());

    return *it;

}

template<int k>
int getimax(map<tuple<int,int,int,int>,int> & spin)
{
    vector<int> thisvector;
    for(typename map<tuple<int,int,int,int>,int>::iterator it = spin.begin(); it != spin.end(); ++it)
    {
        thisvector.push_back((*it).first.get<k>());
    }

    vector<int>::iterator it=max_element(thisvector.begin(), thisvector.end());

    return *it;

}

template<int k>
int getimin(map<tuple<int,int,int,int>,int> & spin)
{
    vector<int> thisvector;
    for(typename map<tuple<int,int,int,int>,int>::iterator it = spin.begin(); it != spin.end(); ++it)
    {
        thisvector.push_back((*it).first.get<k>());
    }

    vector<int>::iterator it=min_element(thisvector.begin(), thisvector.end());

    return *it;

}

std::string exec(std::string command);

std::vector<std::string> split(const std::string &s, char delim, std::vector<std::string> &elems);

void split_is(const std::string &s, char delim, std::vector<std::string> &elems);

void split_is(const std::string &s, std::string delim_regex, std::vector<std::string> &elems);

std::string to_string_is(int n);

std::string to_string_is(double d);

int stoi_is(std::string s);

double stod_is(std::string s);

void convertspintostring(map<tuple<int,int,int,int>,int>& spin, string &spinstring);

void printblock(map<tuple<int,int,int,int>,int> &thisblock);

void printblocklist(vector<map<tuple<int,int,int,int>,int> > &thisblocklist);

void calculate_formation_energy(map< set<tuple<int,int,int,int,int> >, double> J,map< set<tuple<int,int,int,int,int> >, double> mu, map< set<tuple<int,int,int,int,int> >, double>cluster_type_temp, double constant, double mu_constant, double &formation_energy_temp);

void print_poscar_out(tuple<int,int,int,int,int,int> periodicity,map< tuple<int,int,int,int>,int> unitcell,map<int,int> components,bool tern_alg=false);


void convert_to_starting_point_cluster(set<tuple<int,int,int,int,int> > prototype_set,set<tuple<int,int,int,int,int> > &returned_equivalent_set);

void get_concentration_formation_energy_from_spin_and_J_in(tuple<int,int,int,int,int,int> periodicity_now,map< set<tuple<int,int,int,int,int> >, double > mu,map< set<tuple<int,int,int,int,int> >, double > J_fixed_part, map<tuple<int,int,int,int>,int>  spin_now, double constant, double mu_constant,double&formation_energy_out, double &concentration_out);

void evaluate_energy_at_hull(double concentration_now, map<double,double> hull_map, double &energy_at_hull);

void print_poscar_no_vacancy(tuple<int,int,int,int,int,int> periodicity,map< tuple<int,int,int,int>,int> unitcell,map<int,int> components,bool tern_alg=false);

std::string get_file_contents(const char *filename);

struct solver_variable {
    bool use_gradual_introduction_new_variables;
    bool do_monte_carlo;
    tuple<int,int,int,int,int,int> monte_carlo_periodicity;
    double initial_T;
    double final_T;
    double T_decrement;
    long long each_T_how_many_runs;
    bool dedicated1D;
    bool ternary_alg;
    map<set< tuple<int,int,int,int,int,int> >,double > J_input_tern_unaltered;
    map<set< tuple<int,int,int,int,int> >,double > mu_tern_first_group;
    map<set< tuple<int,int,int,int,int> >,double > mu_tern_second_group;
    int model_numbers;
    int thread_number;
    bool ternary_output_states_above_hull;
    double ternary_x_min;
    double ternary_x_max;
    double ternary_y_min;
    double ternary_y_max;
    double ternary_z_min;
    double ternary_z_max;
    double ternary_output_states_above_hull_gap;
    double ternary_output_states_above_hull_ceil;
    bool ternary_debug;
    double mu1;
    double mu2;
    bool use_incomplete_solver;
    string which_incomplete_solver;
    double incomplete_time_cap;
    int N_to_start_incomplete;
    map<tuple<int,int,int,int,int,int>, double > formation_periodic_to_energy_transfer;
    map<tuple<int,int,int,int,int,int>, double > formation_periodic_to_concentration_transfer;
    bool Major_periodicity;
    double prec;
    bool activate_large_cell_algo_binary;
    int large_cell_binary_cubic_periodicity;
    double large_cell_binary_mu_init;
    double large_cell_binary_mu_step;
    double large_cell_binary_mu_final;
    double large_cell_binary_mu_init_todo;
    double large_cell_binary_mu_step_todo;
    double large_cell_binary_mu_final_todo;
    map<double, double > formation_mu_to_energy_transfer;
    map<double, double > formation_mu_to_concentration_transfer;
    map<double, map<tuple<int,int,int,int>,int> > map_mu_to_spin_transfer;
    double binary_x_min;
    double binary_x_max;
    bool Use_simulated_annealing_as_incomplete_solver;

    
    
    
    
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & use_gradual_introduction_new_variables;
        ar & do_monte_carlo;
        ar & monte_carlo_periodicity;
        ar & initial_T;
        ar & final_T;
        ar & T_decrement;
        ar & each_T_how_many_runs;
        ar & dedicated1D;
        ar & ternary_alg;
        ar & J_input_tern_unaltered;
        ar & mu_tern_first_group;
        ar & mu_tern_second_group;
        ar & model_numbers;
        ar & thread_number;
        ar & ternary_output_states_above_hull;
        ar & ternary_x_min;
        ar & ternary_x_max;
        ar & ternary_y_min;
        ar & ternary_y_max;
        ar & ternary_z_min;
        ar & ternary_z_max;
        ar & ternary_output_states_above_hull_gap;
        ar & ternary_output_states_above_hull_ceil;
        ar & ternary_debug;
        ar & mu1;
        ar & mu2;
        ar & use_incomplete_solver;
        ar & which_incomplete_solver;
        ar & incomplete_time_cap;
        ar & N_to_start_incomplete;
        ar & formation_periodic_to_energy_transfer;
        ar & formation_periodic_to_concentration_transfer;
        ar & Major_periodicity;
        ar & prec;
        ar & map_mu_to_spin_transfer;
        ar & activate_large_cell_algo_binary;
        ar & large_cell_binary_cubic_periodicity;
        ar & large_cell_binary_mu_init;
        ar & large_cell_binary_mu_step;
        ar & large_cell_binary_mu_final;
        //note mu here is 1 0 based, will change it when necessaary later ...
        ar &  large_cell_binary_mu_init_todo;
        ar &  large_cell_binary_mu_step_todo;
        ar &  large_cell_binary_mu_final_todo;
        ar &  formation_mu_to_energy_transfer;
        ar &  formation_mu_to_concentration_transfer;
        ar &  binary_x_min;
        ar &  binary_x_max;
        ar &  Use_simulated_annealing_as_incomplete_solver;

        
    }

    
};

struct spin_no_periodic_struct{
public:
    map<tuple<int,int,int,int>,int> spin;
    
    bool operator<(const spin_no_periodic_struct& rhs) const {
        string str1,str2;
        map<tuple<int,int,int,int>,int> spinleft=spin;
        map<tuple<int,int,int,int>,int> spinright=rhs.spin;
        convertspintostring(spinleft, str1);
        convertspintostring(spinright, str2);
        return str1<str2;

    }
    bool operator==(const spin_no_periodic_struct& rhs) const {
        map<tuple<int,int,int,int>,int> spinleft=spin;
        map<tuple<int,int,int,int>,int> spinright=rhs.spin;
        return spinleft==spinright;
    }
    spin_no_periodic_struct() {};
    spin_no_periodic_struct(map<tuple<int,int,int,int>,int> spin_temp)
    {
        spin=spin_temp;
    }
    void print_spin()
    {
        string str1;
        map<tuple<int,int,int,int>,int> spinnow=spin;
        convertspintostring(spinnow, str1);
        cout<<str1;
    }
    
    spin_no_periodic_struct remove_x_neg();
    spin_no_periodic_struct remove_x_pos();
    spin_no_periodic_struct remove_y_neg();
    spin_no_periodic_struct remove_y_pos();
    spin_no_periodic_struct remove_z_neg();
    spin_no_periodic_struct remove_z_pos();

};

void calculate_cluster_type_and_energy_no_periodic(map<set<tuple<int,int,int,int,int> >, double> J_for_proof,double &energy,map<set<tuple<int,int,int,int,int> >, double> &cluster_type_here,map<tuple<int,int,int,int>,int>spin_now)
;

void calculate_cluster_type_and_energy_periodic(map<set<tuple<int,int,int,int,int> >, double> J_for_proof,tuple<int,int,int,int,int,int>periodicity ,double &energy,map<set<tuple<int,int,int,int,int> >, double> &cluster_type_here,map<tuple<int,int,int,int>,int>spin_now);

map<tuple<int,int,int,int>,int> extend_periodic_spin(tuple<int,int,int,int,int,int>periodicity,map<tuple<int,int,int,int>,int>spin_now ,map<set<tuple<int,int,int,int,int> >, double> J);

//template<typename type1>

void construct_set_of_equivalent_from_prototype_set(set<set<tuple<int,int,int,int,int> > > &set_of_equivalent,set<tuple<int,int,int,int,int> > &prototype_set,int x_range,int y_range,int z_range);


struct spin_periodic_struct{
public:
    map<tuple<int,int,int,int>,int> spin;
    tuple<int,int,int,int,int,int> periodicity;
    
    bool operator<(const spin_no_periodic_struct& rhs) const {
        string str1,str2;
        map<tuple<int,int,int,int>,int> spinleft=spin;
        map<tuple<int,int,int,int>,int> spinright=rhs.spin;
        convertspintostring(spinleft, str1);
        convertspintostring(spinright, str2);
        return str1<str2;
        
    }
    bool operator==(const spin_no_periodic_struct& rhs) const {
        map<tuple<int,int,int,int>,int> spinleft=spin;
        map<tuple<int,int,int,int>,int> spinright=rhs.spin;
        return spinleft==spinright;
    }
    spin_periodic_struct() {};
    spin_periodic_struct(map<tuple<int,int,int,int>,int> spin_temp,tuple<int,int,int,int,int,int> input_periodicity)
    {
        spin=spin_temp;
        periodicity=input_periodicity;
    }
    void print_spin()
    {
        string str1;
        map<tuple<int,int,int,int>,int> spinnow=spin;
        convertspintostring(spinnow, str1);
        cout<<str1;
    }

    
};





void calculate_range_from_J(map<set<tuple<int,int,int,int,int> >, double> &J,
                            int &x_range,int &y_range,int &z_range,
                            map<int, int>&component);



void  convert_J_input_tern_unaltered_to_J_tern_in(map<set< tuple<int,int,int,int,int,int> >,double > J_input_tern_unaltered,map<set< tuple<int,int,int,int,int> >,double >& J_tern_in);


void get_concentration_formation_energy_from_spin_and_J_in_ternary(tuple<int,int,int,int,int,int> periodicity_now,map< set<tuple<int,int,int,int,int> >, double > mu_first, map< set<tuple<int,int,int,int,int> >, double > mu_second ,map< set<tuple<int,int,int,int,int> >, double > J_fixed_part, map<tuple<int,int,int,int>,int>  spin_now, double constant, double mu_constant,double&formation_energy_out, tuple<double,double> &concentration_out,map< set<tuple<int,int,int,int,int> >, double > &clustertype_out);

int floor_int_division(int up,int down);




#endif
