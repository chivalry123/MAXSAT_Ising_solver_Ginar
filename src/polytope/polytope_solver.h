#ifndef PTOPE_H
#define PTOPE_H

#include <iostream>
#include <iomanip>
#include <streambuf>

#include <fstream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <map>
#include <vector>
#include <set>
#include <sstream>
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <stdio.h>
#include "serialize_tuple.h"
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/set.hpp>
#include "../ising/tools.h"
#include "eval_queue.h"

using namespace std;
using namespace ::boost::mpi;
using namespace ::boost::algorithm;
using namespace ::boost::tuples;
using namespace ::boost;

void find_polytope(int max_sites, int num_loops, double prec, int mode, double ediff,
                   int max_iter, bool opt_init, int n_rand_init, int save_freq,
                   map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                   set< tuple<int,int,double> > &constraints,
                   std::string prefix, int n_procs,
                   bool load_save, std::string infile, bool include_nonconverged,
                   boost::mpi::environment &env,
                   boost::mpi::communicator &world,        int use_new_pair_terms,
                   int use_new_triplet_terms);

// Solving polytope
void initialize_polytope(int n_rand_init, bool opt_init,
                        map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                        set< map<int, double> > &Jqueued);

void evaluate_queue(int max_sites, int num_loops, double prec, int mode,
                    map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                    set< map<int, double> > &Jqueued, set< map<int, double> > &Jdone,
                    int &n_iter, int save_freq,
                    vector< tuple< map<int, double>, double > > &lowerbound_points,
                    vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                    int n_procs, boost::mpi::environment &env, boost::mpi::communicator &world,int use_new_pair_terms,
                    int use_new_triplet_terms);

void construct_polytope(map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                        vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                        set< int > &polytope_vertices,
                        set< set<int> > &polytope_facets,
                        set< tuple<map<int, double>, double> > &polytope_normals);

void generate_new_J_vectors(map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                            set< tuple<int,int,double> > &constraints,
                            vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                            set< tuple<map<int, double>, double> > &polytope_normals,
                            vector< tuple< map<int, double>, double > > &lowerbound_points,
                            set< map<int, double> > &Jqueued,
                            set< map<int, double> > &Jdone);

double get_lower_bound(map<int, double> Jin, vector< tuple< map<int, double>, double > > &lowerbound_points);

void save_calculations(map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                       vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                       int n_iter);

void save_polytope(std::string prefix, double ediff,
                   map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                   set< tuple<int,int,double> > &constraints,
                   vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                   vector< tuple< map<int, double>, double > > &lowerbound_points,
                   set< int > &polytope_vertices,
                   set< set<int> > &polytope_facets,
                   set< tuple<map<int, double>, double> > &polytope_normals);

void load_polytope(std::string infile, double ediff,
                   bool include_nonconverged,
                   map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                   vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                   vector< tuple< map<int, double>, double > > &lowerbound_points,
                   set< int > &polytope_vertices,
                   set< set<int> > &polytope_facets,
                   set< tuple<map<int, double>, double> > &polytope_normals,
                   set< map<int, double> > &Jdone);

map<int, double> read_J(std::string J_str);

bool check_J(map<int, double> J, set< tuple<int,int,double> > &constraints);

bool check_J_ptope(map<int, double> J, double Eu, set< tuple<int,int,double> > &constraints,
                   vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points);

bool is_converged(double ediff, map<int, map<int, set< tuple<int, int, int, int> > > > &bases,
                                set< tuple<int,int,double> > &constraints,
                               set< tuple<map<int, double>, double> > &polytope_normals,
                               vector< tuple< map<int, double>, map<int, double>, double, double,
                                    tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                               vector< tuple< map<int, double>, double> > &lowerbound_points);

// Printing
void print_basis_set(map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases);

void print_constraints(set< tuple<int,int,double> > &constraints);

void print_J_set(set< map<int, double> > &Jqueued);

void print_J(map<int, double> &J);

void print_Ubound(vector< tuple< map<int, double>, map<int, double>, double, double,
                        tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points);

void print_Lbound(vector< tuple< map<int, double>, double> > &lowerbound_points);

void print_vertices(vector< tuple< map<int, double>, map<int, double>, double, double,
                                   tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                    set< int > &polytope_vertices);


void print_normals_convergence(double ediff,
                               map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                               set< tuple<int,int,double> > &constraints,
                               set< tuple<map<int, double>, double> > &polytope_normals,
                               vector< tuple< map<int, double>, map<int, double>, double, double,
                                    tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                               vector< tuple< map<int, double>, double> > &lowerbound_points,
                               bool print_all,
                               bool print_non_conv,
                               bool print_conv);

// Utility
std::string to_string_is(map<tuple<int,int,int,int>,int> &spin);

std::string to_string_is(tuple<int,int,int,int,int,int> &pc);


void print_J_in_file_non_converged(vector< tuple< map<int, double>, map<int, double>, double, double,
                                   tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > upperbound_points,double ediff,map<int, map<int, set< tuple<int,int,int,int,int> > > > bases);

void update_lowerbound_points(vector< tuple< map<int, double>, map<int, double>, double, double,
                              tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > upperbound_points,vector< tuple< map<int, double>, double> > &lowerbound_points);

#endif
