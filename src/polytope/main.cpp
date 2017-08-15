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
#include <stdio.h>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/set.hpp>

#include "polytope_solver.h"
#include "serialize_tuple.h"
#include "eval_queue.h"
#include "../ising/tools.h"

using namespace std;
using namespace ::boost::mpi;
using namespace ::boost::algorithm;
using namespace ::boost::tuples;
using namespace ::boost;

map<int, map<int, set< tuple<int,int,int,int,int> > > > load_bases(std::string infile);
set< tuple<int, int, double> > load_constraints(std::string infile);

int main(int argc, char* argv[]){
    // Initialize MPI stuff - root process (rank 0) runs most of the polytope solver while the slaves hang around and
    // help with the queue evaluation until everything is done.
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(11);
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    int n_procs = world.size();
    int proc_id = world.rank();

    if (proc_id > 0){
        eval_queue_slave(proc_id, n_procs, env, world);
    }else{
        bool restart = false;
        bool load_nonconverged = false;
        std::string restart_file = "";
        if(argc > 1){
            if(std::string(argv[1]) == "restart"){
                cout << "Restarting from only converged portion of the polytope." << endl;
                restart = true;
            }else if(std::string(argv[1]) == "continue"){
                cout << "Restarting from previously saved polytope (including nonconverged parts)." << endl;
                restart = true;
                load_nonconverged = true;
            }
            restart_file = argv[2];
        }

        time_t start_time;
        time_t end_time;
        time(&start_time);
        cout << "Polytope solver running on " << n_procs << " processes." << endl;
        std::string config_file = "config.in";
        std::string basis_file = "basis.in";
        std::string constraint_file = "constraints.in";

        int max_sites = -1;
        int max_iter = -1;
        int num_loops = -1;
        int mode = -1;
        int save_freq = -1;
        double prec = -1.0;
        double ediff = -1.0;
        bool opt_init = true;
        int n_rand_init = -1;
        int use_new_pair_terms=100;
        int use_new_triplet_terms=100;
        
        
        std::string prefix = "";

        std::string line;
        ifstream infile(config_file.c_str());
        if (infile.is_open()){
            while ( getline (infile, line) ){
                if(line.find("MAX_SITES") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    max_sites = stoi_is(val[1]);
                }else if(line.find("NUM_LOOPS") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    num_loops = stoi_is(val[1]);
                }else if(line.find("MAX_POLYTOPE_ITER") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    max_iter = stoi_is(val[1]);
                }else if(line.find("PREC") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    prec = stod_is(val[1]);
                }else if(line.find("EDIFF") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    ediff = stod_is(val[1]);
                }else if(line.find("MODE") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    mode = stoi_is(val[1]);
                }else if(line.find("N_RAND_INIT") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    n_rand_init = stoi_is(val[1]);
                }else if(line.find("OPT_INIT") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    if(val[1].find("False") != string::npos){
                        opt_init = false;
                    }
                }else if(line.find("SAVE_FREQ") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    save_freq = stoi_is(val[1]);
                }else if(line.find("NAME") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    prefix = val[1];
                }else if(line.find("use_new_pair_terms") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    use_new_pair_terms = stoi_is(val[1]);
                }else if(line.find("use_new_triplet_terms") != string::npos){
                    vector<std::string> val;
                    split_is(line, '=', val);
                    use_new_triplet_terms = stoi_is(val[1]);
                }
                
            }
            infile.close();
        }
        if (max_sites < 0 || num_loops < 0 || prec < 0 || mode < 0 || save_freq < -1 || ediff < 0 ||
            max_iter < 0 || n_rand_init < 0 || prefix == ""){
            cout << "Error reading config file. Exiting." <<endl;
            exit(1);
        }

        map<int, map<int, set< tuple<int,int,int,int,int> > > > bases = load_bases(basis_file);
        set< tuple<int, int, double> > constraints = load_constraints(constraint_file);

        find_polytope(max_sites, num_loops, prec, mode, ediff,
                      max_iter, opt_init, n_rand_init, save_freq,
                      bases, constraints,
                      prefix, n_procs,
                      restart, restart_file, load_nonconverged,
                      env, world,  use_new_pair_terms,
                       use_new_triplet_terms);
        time(&end_time);
        cout << "\nExecution took " << (int)difftime(end_time, start_time) << " seconds on the master thread." << endl;
        cout << "Master Thread Done." << endl;
    }
}

map<int, map<int, set< tuple<int,int,int,int,int> > > > load_bases(std::string infile){
    map<int, map<int, set< tuple<int,int,int,int,int> > > > bases;
    map<int, set< tuple<int,int,int,int,int> > > basis_group;
    set< tuple<int,int,int,int,int> > basis;
    bool in_group = false;
    std::string line;
    ifstream basefile(infile.c_str());
    if (basefile.is_open()){
        while ( getline (basefile, line) ){
            if(line.find("BASES") != string::npos){
                continue;
            }else if(line.c_str()[0] == '#'){
                continue;
            }else if(line.find("SYMMETRY_GROUP") != string::npos){
                if(in_group){
                    bases[bases.size()] = basis_group;
                    basis_group.clear();
                }else{
                    in_group = true;
                }
            }else if(in_group){
                vector<std::string> elements;
                split_is(line, '/', elements);
                for(uint el = 0; el < elements.size(); el++){
                    string str_el = elements[el];
                    vector<std::string> coords;
                    split_is(str_el, ' ', coords);
                    basis.insert(make_tuple(stoi_is(coords[0]), stoi_is(coords[1]), stoi_is(coords[2]), stoi_is(coords[3]), stoi_is(coords[4]) ));
                }
                basis_group[basis_group.size()] = basis;
                basis.clear();
            }
        }
        if(in_group){
            bases[bases.size()] = basis_group;
            basis_group.clear();
            in_group = false;
        }
    }
    return bases;
}

set< tuple<int, int, double> > load_constraints(std::string infile){
    set< tuple<int, int, double> > constraints;
    std::string line;
    ifstream basefile(infile.c_str());
    if (basefile.is_open()){
        while ( getline (basefile, line) ){
            cout << line << endl;
            if(line.find("CONSTRAINTS") != string::npos){
                continue;
            }else if(line.c_str()[0] == '#'){
                continue;
            }else{
                cout <<"reading constraint"<<endl;
                vector<std::string> const_el;
                split_is(line, ' ', const_el);
                constraints.insert(make_tuple(stoi_is(const_el[0]), stoi_is(const_el[1]), stod_is(const_el[2])));
            }
        }
    }
    return constraints;
}
