#include "polytope_solver.h"

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
                   int use_new_triplet_terms){

    srand (time(NULL)); // Random seed
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(11);

    // ----------------------------------------------------------------------------------------------------------
    // Variables
    // ----------------------------------------------------------------------------------------------------------
    // Lower bound points - vector of (J_in, E_lower) points
    vector< tuple< map<int, double>,
                   double                               > > lowerbound_points;

    // Upper bound points - vector of (J_in, sigma_u_out, E_upper, E_lower, unit_cell, periodicity, time) points
    vector< tuple< map<int, double>,
                   map<int, double>,
                   double, double,
                   tuple<int,int,int,int,int,int>,
                   map< tuple<int,int,int,int>, int>,
                   double                               > > upperbound_points;

    // Polytope - vector of (J_in, sigma_u_out, unit_cell) points
    set< int > polytope_vertices;
    set< set<int> > polytope_facets;
    set< tuple< map<int, double>, double> > polytope_normals;

    // Inputs that have been run
    set< map<int, double> > Jdone;

    // Inputs that need to run
    set< map<int, double> > Jqueued;

    // Load save if needed
    if(load_save){
        
        load_polytope(infile, ediff, include_nonconverged, bases, upperbound_points, lowerbound_points,
                      polytope_vertices, polytope_facets, polytope_normals, Jdone);
        
        
        cout<<"\n From WX let's print the J_in files that are not converged\n";
        print_J_in_file_non_converged(upperbound_points, ediff, bases);
        
        cout<<"\n From WX let's update_lowerbound_points \n";
        update_lowerbound_points(upperbound_points,lowerbound_points);
        
// from Wenxuan add one more function in correcting the polytope normals
        
//       debug from Wenxuan load polytope and do analysis
//        cout<<"\ndebug from WX load polytope and do analysis";
//        cout<<"\n let's first print out the upper bound points";
//        for (  int up_count=0 ; up_count<upperbound_points.size(); up_count++) {
//            tuple< map<int, double>, map<int, double>,double, double,tuple<int,int,int,int,int,int>,map< tuple<int,int,int,int>, int>,double> this_point=upperbound_points[up_count];
//            map<int, double> J_now=this_point.get<0>();
//            cout << "\n\t --> J: ";
//            print_J(J_now);
//
//            map<int, double> sigma_now=this_point.get<1>();
//            cout << "\n\t --> SIGMA: ";
//            print_J(sigma_now);
//
//            cout<<"\nEupper: "<<this_point.get<2>()<<" Elower: "<<this_point.get<3>();
//            
//        }
//
//        cout<<"\n################################################\n Then let's look at what the hell is the initial convergence out put" ;
//        print_normals_convergence(ediff, bases, constraints, polytope_normals, upperbound_points, lowerbound_points, true, true, true);
//        
//        
//        cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n Now!! Then let's look at what the hell is the new convergence out put after reconstruct polytope" ;
//        
//        construct_polytope(bases, upperbound_points, polytope_vertices, polytope_facets, polytope_normals);
//        
//        print_normals_convergence(ediff, bases, constraints, polytope_normals, upperbound_points, lowerbound_points, true, true, true);
        
    }

    // Basis set
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYTOPE SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    cout << "Running polytope solver on basis set:" << endl;
    print_basis_set(bases);
    cout << "Under the following constraints:" << endl;
    print_constraints(constraints);
    cout << "\n\n";

    // ----------------------------------------------------------------------------------------------------------
    // Algorithm
    // ----------------------------------------------------------------------------------------------------------
    // 0. Broadcast some key information to all processes
    tuple< map<int, map<int, set< tuple<int,int,int,int,int> > > >, int, int, double, int,int,int> config =
    make_tuple(bases, max_sites, num_loops, prec, mode, use_new_pair_terms,
               use_new_triplet_terms);
    boost::mpi::broadcast(world, config, 0);

    // 1. Initialize first set of J vectors
    if(!load_save){
        cout << "(MASTER) Generating initial set of J-vectors..." << endl;
        initialize_polytope(n_rand_init, opt_init, bases, Jqueued);
    }else{
        cout << "(MASTER) Generating new J vectors ..." << endl;
        generate_new_J_vectors(bases, constraints, upperbound_points, polytope_normals, lowerbound_points, Jqueued, Jdone);
    }

    cout << "Starting iterative polytope solver..." << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    int n_iter = 0;
    do{
        
        // Run through queue of J vectors
        cout << "(MASTER) Evaluating queue ..." << endl;
        evaluate_queue(max_sites, num_loops, prec, mode,
                       bases, Jqueued, Jdone,
                       n_iter, save_freq,
                       lowerbound_points, upperbound_points,
                       n_procs, env, world,   use_new_pair_terms,
                        use_new_triplet_terms);

        cout << "(MASTER) Updating polytope ..." << endl;
        construct_polytope(bases, upperbound_points, polytope_vertices, polytope_facets, polytope_normals);
        
        if (n_iter < max_iter){
            // Update polytope and bounds
    
//            cout<<"\n debug from wenxuan here check polytope convergence";
//            print_normals_convergence(ediff, bases, constraints, polytope_normals, upperbound_points, lowerbound_points, true, true, true);

            // Generate new set of J vectors from non-converged facets
            
            cout << "(MASTER) Generating new J vectors ..." << endl;
            generate_new_J_vectors(bases, constraints, upperbound_points, polytope_normals, lowerbound_points, Jqueued, Jdone);
        }

        
        std::string s_prefix = prefix + +"."+to_string_is(n_iter);
        
        cout << "(MASTER) Saving polytope ..." << endl;
        
        update_lowerbound_points(upperbound_points,lowerbound_points);
        save_polytope(s_prefix, ediff, bases, constraints, upperbound_points, lowerbound_points, polytope_vertices,
                      polytope_facets, polytope_normals);
        
        cout<<"\n From WX let's print the J_in files that are not converged\n";
        print_J_in_file_non_converged(upperbound_points, ediff, bases);
        
        
    }while(!Jqueued.empty() and n_iter < max_iter);

    // Kill remaining processes
    vector< vector< map<int, double> > > J_distrib;
    for (int p = 0; p < n_procs; p++){
        vector< map<int, double> > proc_H;
        map<int, double> Jtemp;
        Jtemp[0] = -10;
        proc_H.push_back(Jtemp);
        J_distrib.push_back(proc_H);
    }

    for (int p = 1; p < n_procs; p++){
        world.send(p, n_iter, J_distrib[p]);
    }

//    cout<<"\nBefore polytope solver finished, let's: construct polytope again, and update lower bound points";
//    construct_polytope(bases, upperbound_points, polytope_vertices, polytope_facets, polytope_normals);
//    update_lowerbound_points(upperbound_points,lowerbound_points);
    
    cout << "--------------------- Polytope solver finished -----------------------" << endl;

    cout << "Polytope vertices: " << endl;
    print_vertices(upperbound_points, polytope_vertices);
    cout <<"\n\n";

    cout << "Polytope facets: " << endl;
    print_normals_convergence(ediff, bases, constraints, polytope_normals, upperbound_points, lowerbound_points, true, true, true);

    std::string s_prefix = prefix + ".final";
    cout << "Saving to " << s_prefix << " ... " << endl;
    save_polytope(s_prefix, ediff, bases, constraints, upperbound_points, lowerbound_points, polytope_vertices, polytope_facets, polytope_normals);
    
    
    
//    cout<<"\n debug from Wenxuan, I am looking at whether upper bound points missed some entry";
//
//    
//    for (  int up_count=0 ; up_count<upperbound_points.size(); up_count++) {
//        tuple< map<int, double>, map<int, double>,double, double,tuple<int,int,int,int,int,int>,map< tuple<int,int,int,int>, int>,double> this_point=upperbound_points[up_count];
//        map<int, double> J_now=this_point.get<0>();
//        cout << "\n\t --> J: ";
//        print_J(J_now);
//        
//        map<int, double> sigma_now=this_point.get<1>();
//        cout << "\n\t --> SIGMA: ";
//        print_J(sigma_now);
//        
//        cout<<"\nEupper: "<<this_point.get<2>()<<" Elower: "<<this_point.get<3>();
//        
//    }
    
    
    

}

void initialize_polytope(int n_rand_init, bool opt_init,
                        map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                        set< map<int, double> > &Jqueued){
    map <int, double> J;
    if(opt_init){
        // Optimal vectors
        double all_norm = sqrt(bases.size());
        for (uint b = 0; b < bases.size(); b++){ J[b] = 1.0/ all_norm; }
        Jqueued.insert(J);
        J.clear();
        for (uint b = 0; b < bases.size(); b++){ J[b] = -1.0/ all_norm; }
        Jqueued.insert(J);
        J.clear();
        map<int, double> Jp;
        map<int, double> Jn;
        for (uint bo = 1; bo < bases.size()-1; bo++){
            for (uint bi = 0; bi < bases.size()-bo; bi++){
                for (uint b = 0; b < bases.size(); b++){
                    Jp[b] = 0;
                    Jn[b] = 0;
                }

                Jp[bi] = 1.0/ sqrt(2.0);
                Jn[bi] = -1.0/ sqrt(2.0);
                Jp[bi+bo] = -1.0/ sqrt(2.0);
                Jn[bi+bo] = 1.0/ sqrt(2.0);

                Jqueued.insert(Jp);
                Jqueued.insert(Jn);
                Jp.clear();
                Jn.clear();
            }
        }
    }

    // Random vectors
    for (int b = 0; b < n_rand_init; b++){
        double norm = 0;
        for (uint i = 0; i < bases.size(); i++){
            J[i] = (rand() % 200) - 100;
            norm += J[i] * J[i];
        }
        norm = sqrt(norm);
        for (uint i = 0; i < bases.size(); i++){ J[i] = J[i] / norm; }
        Jqueued.insert(J);
        J.clear();
    }

}

void evaluate_queue(int max_sites, int num_loops, double prec, int mode,
                    map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                    set< map<int, double> > &Jqueued, set< map<int, double> > &Jdone,
                    int &n_iter, int save_freq,
                    vector< tuple< map<int, double>, double > > &lowerbound_points,
                    vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                    int n_procs, boost::mpi::environment &env, boost::mpi::communicator &world,int use_new_pair_terms,
                    int use_new_triplet_terms){

    while(Jqueued.size() > 0){
        
        vector< vector< map<int, double> > > J_distrib;
        for (int p = 0; p < n_procs; p++){
            vector< map<int, double> > proc_J;
            J_distrib.push_back(proc_J);
        }

        int n_distributed = 0;
        set<map<int, double> >::iterator Jit;
        for (Jit = Jqueued.begin(); Jit != Jqueued.end() && n_distributed < save_freq; ++Jit){
            map<int, double> J = *Jit;
            J_distrib[n_distributed%n_procs].push_back(J);
            Jdone.insert(J);
            n_distributed++;
        }
        Jqueued.erase(Jqueued.begin(), Jit);

        for (int p = 1; p < n_procs; p++){
            world.send(p, n_iter, J_distrib[p]);
        }

        // Temporate storage for UB, LB vectors during this execution; later added to general UB/LB lists.
        vector< tuple< map<int, double>, map<int, double>, double, double,
                                    tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > batch_upperbound_points;
        vector< tuple< map<int, double>, double > > batch_lowerbound_points;

        // Evaluate this batch of solutions through the queue
        eval_queue(max_sites, num_loops, prec, mode, bases, J_distrib[0], batch_lowerbound_points, batch_upperbound_points, 0,  use_new_pair_terms,
                   use_new_triplet_terms);

        // Merge temporary batch queues with master queue, but keep around the batch for wring to outfile.
        lowerbound_points.reserve(lowerbound_points.size() + distance(batch_lowerbound_points.begin(), batch_lowerbound_points.end()));
        lowerbound_points.insert(lowerbound_points.end(), batch_lowerbound_points.begin(), batch_lowerbound_points.end());
        upperbound_points.reserve(upperbound_points.size() + distance(batch_upperbound_points.begin(), batch_upperbound_points.end()));
        upperbound_points.insert(upperbound_points.end(), batch_upperbound_points.begin(), batch_upperbound_points.end());
        for (int p = 1; p < n_procs; p++){
            vector< tuple< map<int, double>, double> > p_lowerbound_points;
            vector< tuple< map<int, double>, map<int, double>, double, double,
                                    tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > p_upperbound_points;
            world.recv(p, n_iter * n_procs * 2 + p, p_lowerbound_points);
            world.recv(p, n_iter * n_procs * 2 + p + 1, p_upperbound_points);
            lowerbound_points.reserve(lowerbound_points.size() + distance(p_lowerbound_points.begin(), p_lowerbound_points.end()));
            lowerbound_points.insert(lowerbound_points.end(), p_lowerbound_points.begin(), p_lowerbound_points.end());
            upperbound_points.reserve(upperbound_points.size() + distance(p_upperbound_points.begin(), p_upperbound_points.end()));
            upperbound_points.insert(upperbound_points.end(), p_upperbound_points.begin(), p_upperbound_points.end());
            batch_upperbound_points.reserve(batch_upperbound_points.size() + distance(p_upperbound_points.begin(), p_upperbound_points.end()));
            batch_upperbound_points.insert(batch_upperbound_points.end(), p_upperbound_points.begin(), p_upperbound_points.end());
        }
        save_calculations(bases, batch_upperbound_points, n_iter);
        n_iter++;
    
    }
    
}

void construct_polytope(map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                        vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                        set< int > &polytope_vertices,
                        set< set<int> > &polytope_facets,
                        set< tuple<map<int, double>, double> > &polytope_normals){
    // Run Qhull
    std::string qhull_input = "";
    qhull_input += to_string_is((int)bases.size());
    qhull_input += "\n";
    qhull_input += to_string_is((int)upperbound_points.size());
    qhull_input += "\n";
    for (uint i=0; i < upperbound_points.size(); i++){
        map<int, double> sigma = boost::get<1>(upperbound_points[i]);
        for (uint b=0; b < bases.size(); b++){
            qhull_input += to_string_is((double)sigma[b]);
            qhull_input += " ";
        }
        qhull_input += "\n";
    }
    ofstream qhull_in;
    qhull_in.open(".qhull_tmp.in");
    qhull_in << qhull_input;
    qhull_in.close();

    std::string qhull_cmd = "qconvex Fx Fv n TO \'.qhull_tmp.out\' < .qhull_tmp.in";
    exec(qhull_cmd);

    std::string result;
    std::string line;
    ifstream outfile(".qhull_tmp.out");
    if (outfile.is_open()){
        while ( getline (outfile,line) ){
            result += line + "\n";
        }
        outfile.close();
    }

    remove(".qhull_tmp.in");
    remove(".qhull_tmp.out");

    // Parse Qhull output
    vector<std::string> lines;
    split_is(result, '\n', lines);

    polytope_vertices.clear();
    polytope_facets.clear();
    polytope_normals.clear();
    
    

    int vertices_start = 1;
    int n_vertices = stoi_is(lines[0]);
    for(int i = vertices_start; i < vertices_start + n_vertices; i++){
        std::string line = lines[i];
        boost::trim(line);
        polytope_vertices.insert(stoi_is(line));
    }

    int facets_start = vertices_start + n_vertices + 1;
    int n_facets = stoi_is(lines[facets_start - 1]);
    for(int i = facets_start; i < facets_start + n_facets; i++){
        std::string line = lines[i];
        boost::trim(line);
        vector<std::string> vertices_str;
        split_is(line, ' ', vertices_str);
        int n_vertices = stoi_is(vertices_str[0]);
        set< int > vertices;
        for(int j = 1; j < 1 + n_vertices; j++){
            vertices.insert(stoi_is(vertices_str[j]));
        }
        polytope_facets.insert(vertices);
    }

    
    

    int normals_start = facets_start + n_facets + 2;
    //int n_dim = stoi_is(lines[normals_start - 2]);
    int n_nfacets = stoi_is(lines[normals_start - 1]);
    for(int i = normals_start; i < normals_start + n_nfacets; i++){
        std::string line = lines[i];
        boost::trim(line);
        vector<std::string> normals_str;
        split_is(line, ' ', normals_str);
        map<int, double> normal;
        
        for(uint j = 0; j < bases.size(); j++){
            normal[j] = stod_is(normals_str[j]);
        }
        
//    Wenxuan: I would straight away do a normalization here
        double normalization_constant=0;
        for(uint j = 0; j < bases.size(); j++){
            normalization_constant+=normal[j]*normal[j];
        }
    
        for(uint j = 0; j < bases.size(); j++){
            normal[j] = normal[j]/normalization_constant;
        }
        
        
//    Wenxuan: I think there might be a problem here, I am going to forget about the offset and calculate the upper bound myself....
//        double offset = stod_is(normals_str[bases.size()]);
        
        double WX_upper_bound=1e10;
//        cout<<"\n Wenxuan debug here 4816428123";
//        cout<<"\n what is normal here";
//        printmap(normal);
        
        for (int ind1=0; ind1<upperbound_points.size(); ind1++) {
            
            tuple< map<int, double>, map<int, double>, double, double,tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> up_point_now=upperbound_points[ind1];
            map<int, double> sigma_now=up_point_now.get<1>();
            
            //            cout<<"\n what is normal here";
            //            printmap(normal);
            //            cout<<"\n what is sigmanow here";
            //            printmap(sigma_now);
            
            
            double upperbound_temp=0;
            for (map<int, double>::iterator it1=normal.begin(); it1!=normal.end(); it1++) {
                int index_now=it1->first;
                double value_now=it1->second;
                double sigma_value=sigma_now[index_now];
                upperbound_temp-=value_now*sigma_value;
            }
            
            //            cout<<"\n what is upper bound temp then?"<<upperbound_temp;
            if (WX_upper_bound>upperbound_temp) {
                WX_upper_bound=upperbound_temp;
            }
            
            
            
        }
        

//        cout<<"\nafter all what is WX_upper_bound"<<WX_upper_bound;

        
        
//        End of wenxuan edit
        
        polytope_normals.insert(make_tuple(normal, WX_upper_bound));

//        polytope_normals.insert(make_tuple(normal, offset));
    }
}

void generate_new_J_vectors(map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                            set< tuple<int,int,double> > &constraints,
                            vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                            set< tuple<map<int, double>, double> > &polytope_normals,
                            vector< tuple< map<int, double>, double > > &lowerbound_points,
                            set< map<int, double> > &Jqueued,
                            set< map<int, double> > &Jdone){
    set<tuple<map<int, double>, double> >::iterator Nit;
    for (Nit = polytope_normals.begin(); Nit != polytope_normals.end(); ++Nit){
        tuple<map<int, double>, double> normal = *Nit;
        map<int,double> n = boost::get<0>(normal);
        double Eu = boost::get<1>(normal);
        // Construct J vector as the negative of the normalized normal
        double norm = 0;
        for (uint i=0; i < n.size(); i++){ norm += n[i] * n[i]; }
        norm = sqrt(norm);
        map<int, double> J;
        for (uint i=0; i < n.size(); i++){
            J[i] = -1.0 * n[i] / norm;
        }

        // Check that we haven't run anything too similar to this J vector
        set<map<int, double> >::iterator Jit;
        bool unique = true;
        for (Jit = Jdone.begin(); Jit != Jdone.end(); ++Jit){
            map<int, double> Jold = *Jit;
            double dot = 0.0;
            for (uint b = 0; b < bases.size(); b++){ dot += Jold[b] * J[b]; }
            if (abs(dot - 1.0) < 1e-13) {
                unique = false;
                break;
            }
        }

        if (!unique){
            continue;
        }

        //map<int, double> Jneg;
        //for (uint i=0; i < J.size(); i++){ Jneg[i] = -1.0*J[i]; }
        if (!check_J_ptope(J, Eu, constraints, upperbound_points)){
            continue;
        }

        Jqueued.insert(J);
        //Jqueued.insert(Jneg);

        J.clear();
    }
}

double get_lower_bound(map<int, double> Jin, vector< tuple< map<int, double>, double > > &lowerbound_points){
    GRBEnv env = GRBEnv();
    GRBModel m = GRBModel(env);
    map<int, GRBVar> vars;


    
    for(uint b = 0; b < Jin.size(); b++){
        vars[b] = m.addVar(0, 10000000, 0, GRB_CONTINUOUS);
    }
    m.update();
    for(uint i = 0; i < lowerbound_points.size(); i++){
        map<int, double> J = boost::get<0>(lowerbound_points[i]);
        double offset =  boost::get<1>(lowerbound_points[i]);
        GRBLinExpr expr = 0;
        for(uint b = 0; b < vars.size(); b++){
            expr += vars[b] * J[b];
        }
        m.addConstr(expr, GRB_GREATER_EQUAL, offset);
    }
    m.update();
    GRBLinExpr objective=0;
    for(uint b=0; b < vars.size(); b++){
        objective+=vars[b]*Jin[b];
    }
    m.setObjective(objective,GRB_MINIMIZE);
    m.getEnv().set(GRB_IntParam_OutputFlag, 0);
    m.optimize();
    map<int, double> sigma;
    for(uint b=0; b < vars.size(); b++){
        sigma[b] = vars[b].get(GRB_DoubleAttr_X);
    }
    double energy = 0;
    for(uint b=0; b < vars.size(); b++){
        energy += sigma[b] * Jin[b];
    }
    return energy;
}

void save_calculations(map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                       vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                       int n_iter){

    std::string out_string ="";
    out_string += "BASES\n";
    for(uint b = 0; b < bases.size(); b++){
        out_string += "SYMMETRY_GROUP ";
        out_string += to_string_is(int(b));
        out_string += "\n";
        for(uint bs = 0; bs < bases[b].size(); bs++){
            set< tuple<int,int,int,int,int> >::iterator Bit;
            for (Bit = bases[b][bs].begin(); Bit != bases[b][bs].end(); ++Bit){
                tuple<int,int,int,int,int> base = *Bit;
                out_string += to_string_is(boost::get<0>(base));
                out_string += " ";
                out_string += to_string_is(boost::get<1>(base));
                out_string += " ";
                out_string += to_string_is(boost::get<2>(base));
                out_string += " ";
                out_string += to_string_is(boost::get<3>(base));
                out_string += " ";
                out_string += to_string_is(boost::get<4>(base));
                out_string += "/";
            }
            out_string += "\n";
        }
    }

    out_string += "CONFIGURATIONS\n";
    for(uint p=0; p < upperbound_points.size(); p++){
        map<int, double> Jin = boost::get<0>(upperbound_points[p]);
        map<int, double> sigma = boost::get<1>(upperbound_points[p]);
        double Eu = boost::get<2>(upperbound_points[p]);
        double El = boost::get<3>(upperbound_points[p]);
        tuple<int,int,int,int,int,int> pc = boost::get<4>(upperbound_points[p]);
        map<tuple<int,int,int,int>, int> uc = boost::get<5>(upperbound_points[p]);
        int timing = boost::get<6>(upperbound_points[p]);

        for(uint b=0; b<bases.size(); b++){
            out_string += to_string_is(Jin[b]);
            if (b<bases.size()-1) out_string += " ";
        }
        out_string +="|";
        for(uint b=0; b<bases.size(); b++){
            out_string += to_string_is(sigma[b]);
            if (b<bases.size()-1) out_string += " ";
        }
        out_string +="|";
        out_string +=to_string_is(Eu);
        out_string +="|";
        out_string +=to_string_is(El);
        out_string +="|";
        out_string +=to_string_is(pc);
        out_string +="|";
        out_string +=to_string_is(uc);
        out_string +="|";
        out_string +=to_string_is(timing);
        out_string +="%\n";
    }
    std::string filename = "configurations_batch_";
    filename += to_string_is(n_iter);
    filename += ".out";
    ofstream outfile;
    outfile.open(filename.c_str());
    outfile << out_string;
    outfile.close();
}

void save_polytope(std::string prefix, double ediff,
                   map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                   set< tuple<int,int,double> > &constraints,
                   vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                   vector< tuple< map<int, double>, double > > &lowerbound_points,
                   set< int > &polytope_vertices,
                   set< set<int> > &polytope_facets,
                   set< tuple<map<int, double>, double> > &polytope_normals){

    std::string out_string ="";
    out_string += "BASES\n";
    for(uint b = 0; b < bases.size(); b++){
        out_string += "SYMMETRY_GROUP ";
        out_string += to_string_is(int(b));
        out_string += "\n";
        for(uint bs = 0; bs < bases[b].size(); bs++){
            set< tuple<int,int,int,int,int> >::iterator Bit;
            for (Bit = bases[b][bs].begin(); Bit != bases[b][bs].end(); ++Bit){
                tuple<int,int,int,int,int> base = *Bit;
                out_string += to_string_is(boost::get<0>(base));
                out_string += " ";
                out_string += to_string_is(boost::get<1>(base));
                out_string += " ";
                out_string += to_string_is(boost::get<2>(base));
                out_string += " ";
                out_string += to_string_is(boost::get<3>(base));
                out_string += " ";
                out_string += to_string_is(boost::get<4>(base));
                out_string += "/";
            }
            out_string += "%\n";
        }
    }

    out_string += "CONFIGURATIONS\n";
    for(uint p=0; p < upperbound_points.size(); p++){
        map<int, double> Jin = boost::get<0>(upperbound_points[p]);
        map<int, double> sigma = boost::get<1>(upperbound_points[p]);
        double Eu = boost::get<2>(upperbound_points[p]);
        double El = boost::get<3>(upperbound_points[p]);
        tuple<int,int,int,int,int,int> pc = boost::get<4>(upperbound_points[p]);
        map<tuple<int,int,int,int>, int> uc = boost::get<5>(upperbound_points[p]);
        int timing = boost::get<6>(upperbound_points[p]);

        for(uint b=0; b<bases.size(); b++){
            out_string += to_string_is(Jin[b]);
            if (b<bases.size()-1) out_string += " ";
        }
        out_string +="|";
        for(uint b=0; b<bases.size(); b++){
            out_string += to_string_is(sigma[b]);
            if (b<bases.size()-1) out_string += " ";
        }
        out_string +="|";
        out_string +=to_string_is(Eu);
        out_string +="|";
        out_string +=to_string_is(El);
        out_string +="|";
        out_string +=to_string_is(pc);
        out_string +="|";
        out_string +=to_string_is(uc);
        out_string +="|";
        out_string +=to_string_is(timing);
        out_string +="%\n";
    }

    out_string += "BOUNDS\n";
    for(uint p=0; p < lowerbound_points.size(); p++){
        map<int, double> Jin = boost::get<0>(lowerbound_points[p]);
        double El = boost::get<1>(lowerbound_points[p]);
        for(uint b=0; b<bases.size(); b++){
            out_string += to_string_is(Jin[b]);
            if (b<bases.size()-1) out_string += " ";
        }
        out_string +="|";
        out_string +=to_string_is(El);
        out_string +="%\n";
    }

    out_string += "VERTICES\n";
    set< int >::iterator Vit;
    for(Vit=polytope_vertices.begin(); Vit != polytope_vertices.end(); ++Vit){
        int vertex = *Vit;
        out_string += to_string_is(vertex) + " ";
    }
    out_string += "%\n";

    out_string += "FACETS\n";
    set< set<int> >::iterator Fit;
    for(Fit=polytope_facets.begin(); Fit != polytope_facets.end(); ++Fit){
        set<int> facet = *Fit;
        set<int>::iterator Fvit;
        for(Fvit=facet.begin(); Fvit != facet.end(); ++Fvit){
            int vertex = *Fvit;
            out_string += to_string_is(vertex) + " ";
        }
        out_string += "%\n";
    }

    out_string += "NORMALS\n";
    set< tuple<map<int, double>, double> >::iterator Nit;
    for(Nit=polytope_normals.begin(); Nit != polytope_normals.end(); ++Nit){
        tuple<map<int, double>, double> normal = *Nit;

        map<int, double> n = boost::get<0>(normal);
        // Construct J vector as the negative of the normalized normal
        double normalization = 0;
        for (uint i=0; i < n.size(); i++){ normalization += n[i] * n[i]; }
        normalization = sqrt(normalization);
        map<int, double> J;
        for (uint i=0; i < n.size(); i++){ J[i] = -1.0 * n[i] / normalization; }

        double Eu = boost::get<1>(normal);
        for(uint b=0; b<bases.size(); b++){
            out_string += to_string_is(n[b]);
            if(b < bases.size()-1) out_string += " ";
        }
        out_string += "|";
        out_string += to_string_is(Eu);
        out_string += "|";
        double El = get_lower_bound(J, lowerbound_points);
        bool valid = check_J_ptope(J, Eu, constraints, upperbound_points);
        if (!valid) {
            out_string += "oor";
            continue;
        }
        bool converged = (abs(El - Eu) < ediff);
        if (converged) { out_string += "true"; }
        else { out_string +="false"; }
        out_string += "%\n";
    }

    std::string filename = prefix;
    filename += ".out";
    ofstream outfile;
    outfile.open(filename.c_str());
    outfile << out_string;
    outfile.close();
}

void load_polytope(std::string infile, double ediff,
                   bool include_nonconverged,
                   map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                   vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                   vector< tuple< map<int, double>, double > > &lowerbound_points,
                   set< int > &polytope_vertices,
                   set< set<int> > &polytope_facets,
                   set< tuple<map<int, double>, double> > &polytope_normals,
                   set< map<int, double> > &Jdone){

    std::string in_str = get_file_contents(infile.c_str());
    vector<std::string> in_sections;
    split_is(in_str, "CONFIGURATIONS", in_sections);
    std::string bases_str = in_sections[0];
    in_str = in_sections[1];
    in_sections.clear();

    split_is(in_str, "BOUNDS", in_sections);
    std::string configurations_str = in_sections[0];
    in_str = in_sections[1];
    in_sections.clear();

    split_is(in_str, "VERTICES", in_sections);
    std::string bounds_str = in_sections[0];
    in_str = in_sections[1];
    in_sections.clear();

    split_is(in_str, "FACETS", in_sections);
    std::string vertices_str = in_sections[0];
    in_str = in_sections[1];
    in_sections.clear();

    split_is(in_str, "NORMALS", in_sections);
    std::string facets_str = in_sections[0];
    std::string normals_str = in_sections[1];
    in_sections.clear();

    // Load basis
    bases.clear();
    map<int, set< tuple<int,int,int,int,int> > > basis_group;
    set< tuple<int,int,int,int,int> > basis;
    bool in_group = false;
    vector<std::string> base_lines;
    split_is(bases_str, '\n', base_lines);
    for(uint l = 0; l < base_lines.size(); l++){
        std::string line = base_lines[l];
        boost::erase_all(line, "%");
        if(line.find("BASES") != string::npos){
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

    // Configurations
    upperbound_points.clear();
    vector<std::string> conf_lines;
    split_is(configurations_str, '%', conf_lines);
    for(uint l = 0; l < conf_lines.size(); l++){
        if (conf_lines[l] == "\n"){ continue; }
        vector<std::string> conf_el;
        split_is(conf_lines[l], '|', conf_el);
        map<int, double> J = read_J(conf_el[0]);
        map<int, double> sigma = read_J(conf_el[1]);
        double Eu = stod_is(conf_el[2]);
        double El = stod_is(conf_el[3]);

        if(!include_nonconverged && abs(Eu-El) > ediff){
            continue;
        }

        // Periodicity
        std::string p_string = conf_el[4];
        boost::erase_all(p_string, "(");
        boost::erase_all(p_string, ")");
        vector<std::string> p_el;
        split_is(p_string, ' ', p_el);
        tuple<int,int,int,int,int,int> p = make_tuple(stoi_is(p_el[0]), stoi_is(p_el[1]), stoi_is(p_el[2]),
                                                      stoi_is(p_el[3]), stoi_is(p_el[4]), stoi_is(p_el[5]));

        // Unit cell
        map< tuple<int,int,int,int>,int> uc;
        vector<std::string> z_blocks;
        split_is(conf_el[5], "\n\n", z_blocks);
        //cout << conf_el[5] << endl;
        for(int z = int(z_blocks.size())-1; z >= 0; z--){
            //cout << "------------------------------------" <<endl;
            //cout << z_blocks[z] << endl;
            vector<std::string> y_blocks;
            split_is(z_blocks[z], '\n', y_blocks);
            for(int y = int(y_blocks.size())-1; y >= 0; y--){
                //cout << "yblock: "<< y_blocks[y] <<endl;
                vector<std::string> x_blocks;
                split_is(y_blocks[y], '(', x_blocks);
                for(int x = 0; x < int(x_blocks.size()); x++){
                    //cout << "xblock: " << x_blocks[x] << endl;
                    std::string uc_str = x_blocks[x];
                    boost::erase_all(uc_str, "(");
                    boost::erase_all(uc_str, ")");
                    vector<std::string> locs;
                    split_is(uc_str, ' ', locs);
                    for(uint p = 0; p < locs.size(); p++){
                        uc[make_tuple(x+1, (y_blocks.size()-y), (z_blocks.size()-z), p+1)] = stoi_is(locs[p]);
                    }
                }
            }
        }


        double time = stod_is(conf_el[6]);

        upperbound_points.push_back(make_tuple(J, sigma, Eu, El, p, uc, time));
        Jdone.insert(J);
    }

    // Bounds
    lowerbound_points.clear();
    vector<std::string> bound_lines;
    split_is(bounds_str, '%', bound_lines);
    for(uint l = 0; l < bound_lines.size(); l++){
        if (bound_lines[l] == "\n"){ continue; }
        vector<std::string> bound_el;
        split_is(bound_lines[l], '|', bound_el);
        map<int, double> J = read_J(bound_el[0]);
        double El = stod_is(bound_el[1]);
        lowerbound_points.push_back(make_tuple(J, El));
    }

    // Vertices
    polytope_vertices.clear();
    vector<std::string> vert_el;
    boost::erase_all(vertices_str, "%");
    boost::erase_all(vertices_str, "\n");
    split_is(vertices_str, ' ', vert_el);
    for(uint v = 0; v < vert_el.size(); v++){
        polytope_vertices.insert(stoi_is(vert_el[v]));
    }

    // Facets
    polytope_facets.clear();
    vector<std::string> facets_lines;
    split_is(facets_str, '%', facets_lines);
    for(uint l = 0; l < facets_lines.size(); l++){
        if (facets_lines[l] == "\n"){ continue; }
        vector<std::string> facet_el;
        std::string facet_line = facets_lines[l];
        boost::erase_all(facet_line, "\n");
        split_is(facet_line, ' ', facet_el);
        set<int> facet;
        for(uint v = 0; v < facet_el.size(); v++){
            facet.insert(stoi_is(facet_el[v]));
        }
        polytope_facets.insert(facet);
    }

    // Normals
    polytope_normals.clear();
    vector<std::string> normals_lines;
    split_is(normals_str, '%', normals_lines);
    for(uint l = 0; l < normals_lines.size(); l++){
        if (normals_lines[l] == "\n"){ continue; }
        vector<std::string> normal_el;
        split_is(normals_lines[l], '|', normal_el);
        map<int,double> J = read_J(normal_el[0]);
        double offset = stod_is(normal_el[1]);
        polytope_normals.insert(make_tuple(J, offset));
    }

}

map<int, double> read_J(std::string J_str){
    boost::erase_all(J_str, "\n");
    vector<std::string> J_el;
    split_is(J_str, ' ', J_el);
    map<int, double> J;
    for(uint i = 0; i < J_el.size(); i++){
        J[int(i)] = stod_is(J_el[i]);
    }
    return J;
}

bool check_J_ptope(map<int, double> J, double Eu, set< tuple<int,int,double> > &constraints,
                   vector< tuple< map<int, double>, map<int, double>, double, double,
                                tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points){
    // Check if J satisfies constraints
    bool valid = check_J(J, constraints); // If the J itself satisfies constraints, great
    // Otherwise, find facet it came from and check if any of its vertices satisfy constraints
    for(uint p = 0; p < upperbound_points.size(); p++){
        map<int, double> sigma = boost::get<1>(upperbound_points[p]);
        double Ep = 0;
        for(uint pp = 0; pp < sigma.size(); pp++){ Ep += sigma[pp] * J[pp]; }
        //cout << "Ep: " << Ep << "; Eu: " << Eu << endl;
        if(abs(Ep - Eu) < 1e-8){ valid = (valid || check_J(boost::get<0>(upperbound_points[p]),constraints)); }
        if(valid) { break; }
    }
    
    return valid;
}

bool check_J(map<int, double> J, set< tuple<int,int,double> > &constraints){
    set< tuple<int,int,double> >::iterator Cit;
    bool valid = true;
    for(Cit = constraints.begin(); Cit != constraints.end(); ++Cit){
        tuple<int,int,double> constraint = *Cit;
        int J1 = boost::get<0>(constraint);
        int J2 = boost::get<1>(constraint);
        double k = boost::get<2>(constraint);
        if (J2 >= 0) {
            valid = valid && (abs(J[J1]) <= k * abs(J[J2]));
        }else{
            valid = valid && J[J1] <= k;
        }
        if(!valid){
            break;
        }
    }
    return valid;
}

std::string to_string_is(map<tuple<int,int,int,int>,int> &spin){
    set<tuple<int,int,int,int> > keylist;
    int imax = getimax<0>(spin);
    int jmax = getimax<1>(spin);
    int zmax = getimax<2>(spin);
    int location_max = getimax<3>(spin);
    int location_min = getimin<3>(spin);
    std::string spinstring="";
    for (int z = zmax; z > 0; z--) {
        spinstring += "\n";
        for (int j = jmax; j > 0; j--) {
            spinstring = spinstring + "\n";
            for (int i = 1; i < imax + 1; i++) {
                spinstring += "(";
                for (int location = location_min; location < location_max + 1; location++) {
                    if (spin.count(make_tuple(i, j, z, location)) == 1){
                        spinstring = spinstring + " " + lexical_cast<string>(spin[(make_tuple(i, j, z, location))]);
                    }else{
                        spinstring = spinstring + " " + "x";
                    }
                }
                spinstring += ")";
            }
        }
    }
    return spinstring;
}

std::string to_string_is(tuple<int,int,int,int,int,int> &pc){
    return boost::lexical_cast<std::string>(pc);
}

void print_basis_set(map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases){
    for (uint b = 0; b < bases.size(); b++){
        map<int, set<tuple<int,int,int,int,int> > > symm_basis = bases[b];
        cout << "Symmetric Basis Group #" << b << "\n";
        for (uint bs = 0; bs < symm_basis.size(); bs++){
            set<tuple<int,int,int,int,int> > basis = symm_basis[bs];
            set<tuple<int,int,int,int,int> >::iterator it;
            cout << "| ";
            for (it = basis.begin(); it != basis.end(); ++it){
                tuple<int,int,int,int,int> cluster_el = *it;
                cout << "(" << boost::get<0>(cluster_el) << ", " << boost::get<1>(cluster_el)
                     << ", "<< boost::get<2>(cluster_el) << ", "<< boost::get<3>(cluster_el)
                     << "): " << boost::get<4>(cluster_el) << " | ";
            }
            cout<<"\n";
        }
        cout<<"\n";
    }
}

void print_constraints(set< tuple<int,int,double> > &constraints){
    set< tuple<int,int,double> >::iterator Cit;
    for(Cit = constraints.begin(); Cit != constraints.end(); ++Cit){
        tuple<int,int,double> constraint = *Cit;
        int J1 = boost::get<0>(constraint);
        int J2 = boost::get<1>(constraint);
        double k = boost::get<2>(constraint);
        if(J2 >= 0){
            cout << "#" << J1 << " < " << k << " * #" << J2 << endl;
        }else{
            cout << "#" << J1 << " < " << k << endl;
        }
    }
}

void print_vertices(vector< tuple< map<int, double>, map<int, double>, double, double,
                                   tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                    set< int > &polytope_vertices){
    set< int >::iterator Vit;
    int vertex_num = 1;
    for (Vit = polytope_vertices.begin(); Vit != polytope_vertices.end(); ++Vit)
    {
        int V = *Vit;
        cout << "#" << vertex_num++ << "(point "<< V << "):\n";
        cout << "J: ";
        print_J(boost::get<0>(upperbound_points[V]));
        cout <<"; SIGMA: ";
        print_J(boost::get<1>(upperbound_points[V]));
        cout <<"\n Energy: " << boost::get<2>(upperbound_points[V]);
        cout <<"\n Periodicity: "<<upperbound_points[V].get<4>();
        cout <<"\n Unit Cell:";
        printblock(boost::get<5>(upperbound_points[V]));
        cout <<"\n\n";
    }
}

bool is_converged(double ediff, map<int, map<int, set< tuple<int, int, int, int> > > > &bases,
                                set< tuple<int,int,double> > &constraints,
                               set< tuple<map<int, double>, double> > &polytope_normals,
                               vector< tuple< map<int, double>, map<int, double>, double, double,
                                    tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                               vector< tuple< map<int, double>, double> > &lowerbound_points){
    set<tuple<map<int, double>, double> >::iterator Nit;
    uint n_converged = 0;
    uint n_relevant = 0;
    for (Nit = polytope_normals.begin(); Nit != polytope_normals.end(); ++Nit){
        tuple<map<int, double>, double> normal = *Nit;
        map<int,double> n = boost::get<0>(normal);
        double Eupper = boost::get<1>(normal);
        // Construct J vector as the negative of the normalized normal
        double normalization = 0;
        for (uint i=0; i < n.size(); i++){ normalization += n[i] * n[i]; }
        normalization = sqrt(normalization);
        map<int, double> J;
        for (uint i=0; i < n.size(); i++){ J[i] = -1.0 * n[i] / normalization; }
        bool valid = check_J_ptope(J, Eupper, constraints, upperbound_points);
        if (!valid) { continue; }
        double Elower = get_lower_bound(J, lowerbound_points);

        bool converged = (abs(Elower - Eupper) < ediff);
        if (converged){
            n_converged += 1;
        }
        n_relevant += 1;
    }
    return (n_converged == polytope_normals.size());
}

void print_normals_convergence(double ediff,
                               map<int, map<int, set< tuple<int,int,int,int,int> > > > &bases,
                               set< tuple<int,int,double> > &constraints,
                               set< tuple<map<int, double>, double> > &polytope_normals,
                               vector< tuple< map<int, double>, map<int, double>, double, double,
                                    tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points,
                               vector< tuple< map<int, double>, double> > &lowerbound_points,
                               bool print_all,
                               bool print_non_conv,
                               bool print_conv){
    set<tuple<map<int, double>, double> >::iterator Nit;
    int n_converged = 0;
    int n_relevant = 0;
    for (Nit = polytope_normals.begin(); Nit != polytope_normals.end(); ++Nit){
        tuple<map<int, double>, double> normal = *Nit;
        map<int,double> n = boost::get<0>(normal);
        double Eupper = boost::get<1>(normal);
        // Construct J vector as the negative of the normalized normal
        double normalization = 0;
        for (uint i=0; i < n.size(); i++){ normalization += n[i] * n[i]; }
        normalization = sqrt(normalization);
        map<int, double> J;
        for (uint i=0; i < n.size(); i++){ J[i] = -1.0 * n[i] / normalization; }
        
        
        //debug from Wenxuan, a key here forget to normalize Eupper!
        Eupper=Eupper/normalization;
        
        // Check if this facet is relevant
        bool valid = check_J_ptope(J, Eupper, constraints, upperbound_points);
        if (!valid) { continue; }

        
        // Check if this facet is converged
        double Elower = get_lower_bound(J, lowerbound_points);
        bool converged = (abs(Elower - Eupper) < ediff);
        if ((converged && print_conv) || (!converged && print_non_conv) || print_all){
            cout << "J: ";
            print_J(J);
            cout << "\n\t\tEu: " << Eupper<< "\t\tEl: "<< Elower << "; dE: " << Eupper-Elower  << endl;
        }
        if (converged){
            if(print_all) cout <<"\t Converged: TRUE\n" << endl;
            n_converged += 1;
            n_relevant += 1;
        }else{
            if(print_all) cout <<"\t Converged: FALSE\n" <<endl;
            n_relevant += 1;
            
            
            cout<<"\n debug from WX ";
            cout<<"\n let's extract this J from upper bound points";
            int index=0;
            map<int, double> supposed_sigma,supposed_J;
            double supposed_UB=0;
            double supposed_LB=0;
            bool extract_successful=false;
            for (int i=0; i<upperbound_points.size(); i++) {
                map<int, double> J_temp=upperbound_points[i].get<0>();
                double diff=0;
                for (int j=0; j<J_temp.size(); j++) {
                    diff+=fabs(J_temp[j]-J[j]);
                }
                if (diff<1e-7) {
                    extract_successful=true;
                    supposed_J=upperbound_points[i].get<0>();
                    supposed_sigma=upperbound_points[i].get<1>();
                    supposed_UB=upperbound_points[i].get<2>();
                    supposed_LB=upperbound_points[i].get<3>();
                    index=i;
                    break;
                }
            }
            
            if (extract_successful) {
                cout<<"\nextraction sucessful, supposed_UB is"<<supposed_UB<<"supposed_LB"<<supposed_LB<<" supposed_sigma is\n";
                print_J(supposed_sigma);
                cout<<"\n then let's check which lower bound calculation is wrong before";
                bool violated_found=false;
                for (int i=0; i<lowerbound_points.size(); i++) {
                    //check if this sigma has lower energy than the lower bound
                    double temp_energy=0;
                    map<int,double> LB_J=lowerbound_points[i].get<0>();
                    double LB=lowerbound_points[i].get<1>();
                    for (int j=0; j<LB_J.size(); j++) {
                        temp_energy+=supposed_sigma[j]*LB_J[j];
                    }
                    
                    if (temp_energy<LB-1e-5) {
                        cout<<"\n violated lower bound is found at index:"<<i;
                        cout<<"\n The corresponding J is\n";
                        print_J(LB_J);
                        cout<<"\n the supposed LB is: "<<LB<<" while we arrive at energy"<<temp_energy<<endl;
                        violated_found=true;
                        break;
                    }
                }
                if (violated_found==false) {
                    cout<<"\n no violation is found";
                }
                
            }
            else{
                cout<<"\nextraction unsuccessful no previous calculation is done before";
            }
            
            
            
            
        }
    }
    if(print_all)
        cout << "Total converged facets: " << n_converged << "/" << n_relevant << endl;
}

void print_Ubound(vector< tuple< map<int, double>, map<int, double>, double, double,
                        tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > &upperbound_points){
    for(uint i = 0; i < upperbound_points.size(); i++){
        map<int, double> J = boost::get<0>(upperbound_points[i]);
        map<int, double> sigma = boost::get<1>(upperbound_points[i]);
        double energy = boost::get<2>(upperbound_points[i]);
        map<tuple<int,int,int,int>, int> unitcell =  boost::get<5>(upperbound_points[i]);
        cout << "Jin: ";
        print_J(J);
        cout <<"\nSIGMA: ";
        print_J(sigma);
        cout <<"\nEnergy: " << energy;
        cout <<"\nUnit Cell:";
        printblock(unitcell);
        cout <<"\n\n";
    }
}

void print_Lbound(vector< tuple< map<int, double>, double> > &lowerbound_points){
    for(uint i = 0; i < lowerbound_points.size(); i++){
        map<int, double> J = boost::get<0>(lowerbound_points[i]);
        double offset =  boost::get<1>(lowerbound_points[i]);
        cout << "Jin: ";
        print_J(J);
        cout <<"\nElower: " << offset <<"\n\n";
    }
}

void print_J_set(set< map<int, double> > &Jset){
    set< map<int, double> >::iterator Jit;
    for (Jit = Jset.begin(); Jit != Jset.end(); ++Jit){
        map<int, double> J = *Jit;
        print_J(J);
        cout << "\n";
    }
}

void print_J(map<int, double> &J){
    cout << "[";
    for (uint i=0; i<J.size(); i++){
        if (J[i] >= 0.0) cout << " ";
        cout << J[i];
        if (i < J.size()-1){ cout << ", "; }
        else{ cout <<"]"; }
    }
}


void print_J_in_file_non_converged(vector< tuple< map<int, double>, map<int, double>, double, double,
                                   tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > upperbound_points,double ediff,map<int, map<int, set< tuple<int,int,int,int,int> > > > bases)
{
    system("rm -r nonconverged_J");
    system("mkdir nonconverged_J");
    
    cout<<"\nediff is "<<ediff;
    
    for (int i=0; i<upperbound_points.size(); i++) {
        double UB=upperbound_points[i].get<2>();
        double LB=upperbound_points[i].get<3>();
//        cout<<"\n UB points number i="<<i<<" UB is:"<<UB<<" LB is:"<<LB<<" fabs(UB-LB)="<<fabs(UB-LB)<<endl;
        
        const char *path="./nonconverged_J/test_file";
        std::ofstream test_file;
        test_file.open(path);
        test_file<<"123";
        
        
        if (fabs(UB-LB)<ediff*0.01) {
            continue;
        }
        else{
            std::ostringstream ss1;
            ss1<<std::setfill('0')<< std::setw(5)<<i;
            
            string J_in_file_name="J_in.in"+ss1.str();
            string J_in_file_full_path="./nonconverged_J/"+J_in_file_name;
            const char *path=J_in_file_full_path.c_str();
            std::ofstream J_in_file;
            J_in_file.open(path);
            
            
            
            stringstream J_in_contents;
            J_in_contents<<setprecision(15);
            J_in_contents<<"Constant 0";
            
            int count_cluster=0;
            
            map<int, double> J_in_temp=upperbound_points[i].get<0>();
            
            for (map<int, double>::iterator it1=J_in_temp.begin(); it1!=J_in_temp.end(); it1++) {
                int outer_index=it1->first;
                double outer_value=it1->second;
                
                map<int, set< tuple<int,int,int,int,int> > >  sub_basis=bases[outer_index];
                for (map<int, set< tuple<int,int,int,int,int > > >::iterator it2=sub_basis.begin(); it2!=sub_basis.end(); it2++) {
                    int inner_index=it2->first;
                    set< tuple<int,int,int,int,int> >  corresponding_set=it2->second;
                    
                    J_in_contents<<"\n\nCluster "<<count_cluster<<"\n";
                    count_cluster++;
                    for (set< tuple<int,int,int,int,int> >::iterator it3=corresponding_set.begin(); it3!=corresponding_set.end(); it3++) {
                        tuple<int,int,int,int,int> tuple_now=*it3;
                        J_in_contents<<tuple_now.get<0>()<<","<<tuple_now.get<1>()<<","<<tuple_now.get<2>()<<","<<tuple_now.get<3>()<<","<<tuple_now.get<4>()<<"    ";
                    }
                    J_in_contents<<"\nJ="<<outer_value;
                    
                    
                }
                
                
            }
 
            
            J_in_file<<J_in_contents.str();
            
            
        }
        
    }
    
    
    
}


void update_lowerbound_points(vector< tuple< map<int, double>, map<int, double>, double, double,
                              tuple<int,int,int,int,int,int>, map< tuple<int,int,int,int>,int>,double> > upperbound_points,vector< tuple< map<int, double>, double> > &lowerbound_points)
{
    lowerbound_points.clear();
    for (int i=0; i<upperbound_points.size(); i++) {
        map<int, double> temp_J=upperbound_points[i].get<0>();
        double LB=upperbound_points[i].get<3>();
        lowerbound_points.push_back(make_tuple(temp_J,LB));
    }
    
}


