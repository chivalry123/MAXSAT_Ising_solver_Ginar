/* This module defines the lower bound solver for the generalized Ising model.
 *
 * Author: Wenxuan Huang
 * Maintainer: Wenxuan Huang, Daniil Kitchaev
 * Date: 15 March, 2015
 *
 * Copyright: Wenxuan Huang (C) 2014, All rights reserved.
 */

#ifndef TREEOFDEVIL_H
#define TREEOFDEVIL_H

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include "gurobi_c++.h"
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"
#include <boost/lexical_cast.hpp>
#include "tools.h"

using namespace std;
using namespace ::boost::tuples;
using namespace ::boost;




void get_lowerbound(int a0,int a1,int a2,int a3,int a4,int a5,
                 map<set<tuple<int,int,int,int,int> >, double> &J,
                 int xrange,int yrange,int zrange,
                 map<int,int> componentnumber,
                 map<tuple<int,int,int,int>,int> &spin,
                 double &averageenergy,
                 map<tuple<int,int,int,int,int,int>, map<set<tuple<int,int,int,int,int> >, double> > &clustertypeperiodic,
                 map<set<tuple<int,int,int,int,int> >, double> &J_for_proof,
                 std::string id,
                 bool pseudo_mode,
                 bool pseudo_mode_with_proof,
                 bool basic_exact_mode,
                 bool obscenely_verbose,
                 double limit_dimension, bool use_level_method, bool use_weighted_dual_average,solver_variable &global_parameters);

void get_lowerbound_sub(int a0,int a1,int a2,int a3,int a4,int a5,
                    map<set<tuple<int,int,int,int,int> >, double> &J,
                    int xrange,int yrange,int zrange,
                    map<int,int> componentnumber,
                    map<tuple<int,int,int,int>,int> &spin,
                    double &averageenergy,
                    map<tuple<int,int,int,int,int,int>, map<set<tuple<int,int,int,int,int> >, double> > &clustertypeperiodic,
                    double LB_so_far,
                    std::string id,
                    bool stop_at_first_iteration,
                    bool stop_when_upperbound_smaller_than_predicted_lower_bound,
                    bool stop_till_exact_result,
                    bool obscenely_verbose,solver_variable &global_parameters);

void graph_model_and_prototype_set_construction(bool obscenely_verbose,set<spin_no_periodic_struct> lowest_energy_spin_struct_set, int x_range, int y_range, int z_range, bool &shall_i_return,long long model_update_loops, bool &warm_restart,set<spin_no_periodic_struct> &warm_start_spin_structs,map<set<tuple<int,int,int,int,int> >, double> &warm_start_J, map<set<tuple<int,int,int,int,int> >, double> J_input_best_memory,map<long long, set<set<tuple<int,int,int,int,int> > > > index_to_equivalent_sets,map<spin_no_periodic_struct, double> PI_value,map<set<set<tuple<int,int,int,int,int> > >,long long> revert_index_to_equivalent_sets,  set<tuple<int,int,int,int,int> > &prototype_set_here,set<long long> significant_index);


void dedicated1D_alg(map<set<tuple<int,int,int,int,int> >, double> J, int x_range,map<int, int> component ,spin_periodic_struct &spin_struct_min,string id, solver_variable global_parameters, double &lowerbound_from_compat,bool obscenely_verbose);

#endif

