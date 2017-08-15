#ifndef EVALQUEUE_H
#define EVALQUEUE_H

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <map>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/list.hpp>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "serialize_tuple.h"
#include "../ising/solver.h"
#include "../ising/tools.h"
#include "polytope_solver.h"

using namespace std;
using namespace ::boost::mpi;
using namespace ::boost::algorithm;
using namespace ::boost::tuples;
using namespace ::boost;

void eval_queue(int max_sites, int num_loops, double prec, int mode,
                map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                vector< map<int, double> > proc_J,
                vector< tuple< map<int, double>, double > > &lowerbound_points,
                vector< tuple< map<int, double>, map<int, double>, double, double,
                               tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>, int>, double> > &upperbound_points,
                int proc,        int use_new_pair_terms,
                int use_new_triplet_terms);

void eval_queue_slave(int proc_id, int n_procs, boost::mpi::environment &env, boost::mpi::communicator &world);
#endif
