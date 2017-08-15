/* This module defines the general Ising solver procedure and provides
 * an interface for calling the solver with a defined Hamiltonian, in the
 * form of a set of ECI's, and a constraint on the maximum size of the
 * desired solution.
 *
 * Author: Wenxuan Huang
 * Maintainer: Wenxuan Huang, Daniil Kitchaev
 * Date: 15 March, 2015
 *
 * Copyright: Wenxuan Huang (C) 2014, All rights reserved.
 */

#ifndef SOLVER_H
#define SOLVER_H

#include <unistd.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <string>
#include <vector>
#include <map>
#include <set>
#include "gurobi_c++.h"
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"
#include <boost/lexical_cast.hpp>
#include "periodicfunction.h"
#include "treeofdevil.h"
#include "tools.h"
#include <iomanip>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "boost/thread.hpp"
#include "common.h"



using namespace std;
using namespace ::boost::tuples;
using namespace ::boost;

void run_solver(int max_sites,
                map< set< tuple<int,int,int,int,int> >, double> &J,
                map< set< tuple<int,int,int,int,int> >, double> &lowerboundclustertype,
                map< set< tuple<int,int,int,int,int> >, double> &upperboundclustertype,
                map< tuple<int,int,int,int>, int> &cellrepresentation,
                double &lower_bound,
                double &upper_bound,
                double &exact_lower_bound,
                map< tuple<int,int,int,int>, int> &unitcell,
                tuple <int,int,int,int,int,int> &periodicity,
                map< set< tuple<int,int,int,int,int> >, double> &J_for_proof,
                std::string id,
                map< set< tuple<int,int,int,int,int> >, double> &mu,
                double mu_constant,
                bool work_with_mu,
                double &formation_energy_UB,
                double &formation_energy_LB, bool new_cluster_algorithm,int use_new_pair_terms,int use_new_triplet_terms, map<tuple<int,int,int,int,int,int>, map<tuple<int,int,int,int>,int> > &map_periodicity_to_spin, bool use_level_method,bool use_weighted_dual_average,solver_variable &global_parameters,
                double prec = 1000000.0,
                int num_loops = 4,
                bool basic_exact_mode = false,
                bool pseudo_mode = true,
                bool pseudo_mode_with_proof = false,
                bool verbose = true,
                bool very_verbose = false,
                bool obscenely_verbose = false,
                bool input_PRIM_output_PRIMOUT_mode=false,
                double limit_dimension=1000000,
                double constant=0,map< set< tuple<int,int,int,int,int> >, double> J_fixed_part= map< set< tuple<int,int,int,int,int> >, double>()
                );

void corecode(map<int, int> component,
              int x_range, int y_range, int z_range,
              int max_sites,
              int loopnumber,
              map<set<tuple<int,int,int,int,int> >, double> &J,
              map< set< tuple<int,int,int,int,int> >, double> &J_fixed_part,
              map<set<tuple<int,int,int,int,int> >, double> &lowerboundclustertype,
              map<set<tuple<int,int,int,int,int> >, double> &upperboundclustertype,
              map<tuple<int,int,int,int>,int> &cellrepresentation,
              double &lower_bound, double &upper_bound,
              map<tuple<int,int,int,int>,int> &unitcell,
              tuple<int,int,int,int,int,int> &periodicity,
              map<set<tuple<int,int,int,int,int> >, double> &J_for_proof,
              std::string id,
              bool pseudo_mode,
              bool pseudo_mode_with_proof,
              bool basic_exact_mode,
              bool very_verbose,
              bool obscenely_verbose,
              double limit_dimension,bool new_cluster_algorithm,int use_new_pair_terms,int use_new_triplet_terms, map<tuple<int,int,int,int,int,int>, map<tuple<int,int,int,int>,int> > &map_periodicity_to_spin, bool use_level_method,bool use_weighted_dual_average ,solver_variable &global_parameters,
              double constant=0,double mu_constant=0);



void simulated_annealing(map< set<tuple<int,int,int,int,int> >, double> J,solver_variable global_parameters,double &monte_energy, spin_periodic_struct & monte_spin_struct,double constant);

#endif
