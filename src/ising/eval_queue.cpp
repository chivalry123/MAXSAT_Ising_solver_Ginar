#include "eval_queue.h"

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
                int use_new_triplet_terms) {

    bool pseudo_mode = false;
    bool exact_mode = false;
    bool proof_mode = false;
    if (mode == 0){
        exact_mode = true;
    }else if (mode == 1){
        pseudo_mode = true;
    }else if (mode == 2){
        pseudo_mode = true;
        proof_mode = true;
    }

    for(uint Ji = 0; Ji < proc_J.size(); Ji++){
        // Initialize run
        double lower_bound_J, exact_lower_bound, upper_bound_J;
        map< set< tuple<int,int,int,int,int> >, double> H, lowerboundclustertype,
                                                        upperboundclustertype, J_for_proof;
        map< tuple<int,int,int,int>,int> unitcell;
        tuple<int,int,int,int,int,int> periodicity;
        map< tuple<int,int,int,int>, int> cellrepresentation;
        map<int, double> J = proc_J[Ji];

        for (uint b = 0; b < bases.size(); b++){
            for (uint bs = 0; bs < bases[b].size(); bs++){
                H[bases[b][bs]] = J[b];
            }
        }

        cout << "(PID " << proc <<")[" << (Ji+1) << "/" << proc_J.size() <<"] Jin: ";
        print_J(J);
        cout << "\t ... \t";
        int run_id = Ji;
        std::string id = to_string_is(proc) + "."+to_string_is(run_id);

        time_t start_time;
        time(&start_time);

//        run_solver(max_sites, H, lowerboundclustertype, upperboundclustertype, cellrepresentation,
//                    lower_bound_J, exact_lower_bound, upper_bound_J, unitcell, periodicity,
//                    J_for_proof, id, prec, num_loops,
//                    exact_mode, pseudo_mode, proof_mode, false, false, false);
//        
//        
//        run_solver(int max_sites,
//                   map< set< tuple<int,int,int,int,int> >, double> &J,
//                   map< set< tuple<int,int,int,int,int> >, double> &lowerboundclustertype,
//                   map< set< tuple<int,int,int,int,int> >, double> &upperboundclustertype,
//                   map< tuple<int,int,int,int>, int> &cellrepresentation,
//                   double &lower_bound,
//                   double &exact_lower_bound,
//                   double &upper_bound,
//                   map< tuple<int,int,int,int>, int> &unitcell,
//                   tuple<int,int,int,int,int,int> &periodicity,
//                   map< set< tuple<int,int,int,int,int> >, double> &J_for_proof,
//                   std::string id,
//                   double prec,
//                   int num_loops,
//                   bool basic_exact_mode,
//                   bool pseudo_mode,
//                   bool pseudo_mode_with_proof,
//                   bool verbose,
//                   bool very_verbose,
//                   bool obscenely_verbose);
        

        map< set< tuple<int,int,int,int,int> >, double> mu;
        double mu_constant=0;
        bool work_with_mu=0;
        double formation_energy_UB=0;
        double formation_energy_LB=0;
        bool new_cluster_algorithm=1;
        map<tuple<int,int,int,int,int,int>, map<tuple<int,int,int,int>,int> > map_periodicity_to_spin;
        bool use_level_method=0;
        bool use_weighted_dual_average=0;
        bool input_PRIM_output_PRIMOUT_mode=false;
        double limit_dimension=1000000;
        double constant=0;
        solver_variable global_parameters;
        global_parameters.use_gradual_introduction_new_variables=true;
        bool obscenely_verbose=false;
        
        run_solver(max_sites, H, lowerboundclustertype, upperboundclustertype, cellrepresentation,
                   lower_bound_J, upper_bound_J, exact_lower_bound, unitcell, periodicity,
                   J_for_proof, id, mu,
                    mu_constant,
                    work_with_mu,
                    formation_energy_UB,
                   formation_energy_LB,  new_cluster_algorithm, use_new_pair_terms, use_new_triplet_terms, map_periodicity_to_spin,  use_level_method, use_weighted_dual_average,global_parameters, prec, num_loops,
                   exact_mode, pseudo_mode, proof_mode, false, false, obscenely_verbose,input_PRIM_output_PRIMOUT_mode,limit_dimension,constant);
        
//        void run_solver(int max_sites,
//                        map< set< tuple<int,int,int,int,int> >, double> &J,
//                        map< set< tuple<int,int,int,int,int> >, double> &lowerboundclustertype,
//                        map< set< tuple<int,int,int,int,int> >, double> &upperboundclustertype,
//                        map< tuple<int,int,int,int>, int> &cellrepresentation,
//                        double &lower_bound,
//                        double &exact_lower_bound,
//                        double &upper_bound,
//                        map< tuple<int,int,int,int>, int> &unitcell,
//                        tuple <int,int,int,int,int,int> &periodicity,
//                        map< set< tuple<int,int,int,int,int> >, double> &J_for_proof,
//                        std::string id,
//                        map< set< tuple<int,int,int,int,int> >, double> &mu,
//                        double mu_constant,
//                        bool work_with_mu,
//                        double &formation_energy_UB,
//                        double &formation_energy_LB, bool new_cluster_algorithm,int use_new_pair_terms,int use_new_triplet_terms, map<tuple<int,int,int,int,int,int>, map<tuple<int,int,int,int>,int> > &map_periodicity_to_spin, bool use_level_method,bool use_weighted_dual_average,
//                        double prec = 1000000.0,
//                        int num_loops = 4,
//                        bool basic_exact_mode = false,
//                        bool pseudo_mode = true,
//                        bool pseudo_mode_with_proof = false,
//                        bool verbose = true,
//                        bool very_verbose = false,
//                        bool obscenely_verbose = false,
//                        bool input_PRIM_output_PRIMOUT_mode=false,
//                        double limit_dimension=1000000,
//                        double constant=0)
//        
//        
        

        time_t end_time;
        time(&end_time);

        // Process results
        map<int, double> sigma;
        for (uint b = 0; b < bases.size(); b++){
            sigma[b] = 0;
            for (uint bs = 0; bs < bases[b].size(); bs++){
                sigma[b] += upperboundclustertype[bases[b][bs]];
            }
        }

        lowerbound_points.push_back(make_tuple(J, lower_bound_J));
        upperbound_points.push_back(make_tuple(J, sigma, upper_bound_J, lower_bound_J,
                                               periodicity, unitcell, (double)difftime(end_time, start_time)));

        cout << "\n\t --> SIGMA: ";
        print_J(sigma);
        cout << "\t\t Eu: "<< upper_bound_J <<"\tEl: "<< lower_bound_J << "\tdE: "<< upper_bound_J - lower_bound_J << endl;
    }
}

void eval_queue_slave(int proc_id, int n_procs, boost::mpi::environment &env, boost::mpi::communicator &world){
    tuple< map<int, map<int, set< tuple<int,int,int,int,int> > > >, int, int, double, int,int,int> config;
    boost::mpi::broadcast(world, config, 0);
    map<int, map<int, set< tuple<int,int,int,int,int> > > > bases = boost::get<0>(config);
    int max_sites = boost::get<1>(config);
    int num_loops = boost::get<2>(config);
    double prec = boost::get<3>(config);
    int mode = boost::get<4>(config);
    int use_new_pair_terms= boost::get<5>(config);
    int use_new_triplet_terms=boost::get<6>(config);
    
    int MAX_LOOPS = 100000;
    int n_iter = 0;
    while(n_iter < MAX_LOOPS){
        vector< map<int, double> > proc_J;
        world.recv(0, n_iter, proc_J);
        unsigned int microseconds = proc_id * 10000;
        usleep(microseconds);

        //cout << "Process " << proc_id << " received " << proc_H.size() <(int)difftime(end_time, start_time)< " jobs to run." << endl;
        if (proc_J.size() > 0 && proc_J[0][0] == -10){
            cout << "(PID " << proc_id << ") Slave Thread Done." << endl;
            return;
        }

        vector< tuple< map<int, double>, double> > lowerbound_points;

        vector< tuple< map<int, double>, map<int, double>, double, double,
                       tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>, double> > upperbound_points;

        eval_queue(max_sites, num_loops, prec, mode,
                    bases, proc_J, lowerbound_points, upperbound_points, proc_id,use_new_pair_terms,use_new_triplet_terms);

        world.send(0, n_iter * n_procs * 2 + proc_id, lowerbound_points);
        world.send(0, n_iter * n_procs * 2 + proc_id + 1, upperbound_points);
        n_iter++;
    }
}

