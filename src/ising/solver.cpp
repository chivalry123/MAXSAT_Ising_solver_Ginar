#include "solver.h"

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
                double &formation_energy_LB,bool new_cluster_algorithm,int use_new_pair_terms,int use_new_triplet_terms, map<tuple<int,int,int,int,int,int>, map<tuple<int,int,int,int>,int> > &map_periodicity_to_spin, bool use_level_method,bool use_weighted_dual_average,solver_variable &global_parameters,
                double prec,
                int num_loops,
                bool basic_exact_mode,
                bool pseudo_mode,
                bool pseudo_mode_with_proof,
                bool verbose,
                bool very_verbose,
                bool obscenely_verbose,
                bool input_PRIM_output_PRIMOUT_mode,
                double limit_dimension,
                double constant,map< set< tuple<int,int,int,int,int> >, double> J_fixed_part){
// J_fixed_part here is only for binary mu
//    char buffer[8192];
//    
//    setvbuf(stdout, buffer, _IOFBF, sizeof(buffer));

    if (!((basic_exact_mode && !pseudo_mode && !pseudo_mode_with_proof) ||
           (!basic_exact_mode && pseudo_mode && pseudo_mode_with_proof) ||
           (!basic_exact_mode && pseudo_mode && !pseudo_mode_with_proof))){
        cout << "Invalid choice of solver mode! Exiting." << endl;
        exit(1);
    }
    
    if (obscenely_verbose) {
        cout<<"\n J_in is:";
        printmapfromsets(J);
    }

    // If verbose, print back basic input
    if (verbose){
        cout << "(" << id <<") %%%%%%%%%%%%%%%%%%%% Solver %%%%%%%%%%%%%%%%%%%%%" << endl;
        cout << "(" << id <<") -------------------- Input ----------------------" << endl; 
        if(basic_exact_mode){
            cout << "Basic exact mode" << endl;
        }else if (pseudo_mode && !pseudo_mode_with_proof){
            cout << "Pseudo mode (no proof)" << endl;
        }else if (pseudo_mode && pseudo_mode_with_proof){
            cout << "Pseudo mode (with proof)" << endl;
        }
        cout << "Precision: " << prec << endl;
        cout << "(" << id <<") J: " << endl;
//        printmapfromsets(J);
        cout << "(" << id <<") Solver: " << endl;
        cout << "(" << id <<") -------------------------------------------------" << endl;
    }

    // Convert J to integers
    map< set< tuple<int,int,int,int,int> >, double> Ji;
    map< set< tuple<int,int,int,int,int> >, double>::iterator Jit;
    for (Jit = J.begin(); Jit != J.end(); ++Jit){
//        Ji[Jit->first] = (double)((long long)(Jit->second / prec));
        Ji[Jit->first] = Jit->second ;

    }

    // Calculate basic block size and dimension (in species)
    map<int,int> components;
    int x_range, y_range, z_range;
    calculate_range_from_J(Ji, x_range, y_range, z_range, components);

    // Solve optimization problem
    if(verbose){
        cout << "(" << id <<") Solving optimization problem ..." << endl;
    }
//    cout << "In corecode x_range, y_range, z_range "<<make_tuple(x_range,y_range,z_range);

    corecode(components, x_range, y_range, z_range, max_sites, num_loops,
             Ji,J_fixed_part,
             lowerboundclustertype, upperboundclustertype, cellrepresentation,
             lower_bound, upper_bound, unitcell, periodicity, J_for_proof,
             id, pseudo_mode, pseudo_mode_with_proof, basic_exact_mode,
             very_verbose, obscenely_verbose,limit_dimension,new_cluster_algorithm,use_new_pair_terms,use_new_triplet_terms,map_periodicity_to_spin,  use_level_method, use_weighted_dual_average,global_parameters,constant, mu_constant);

    // If necessary, prove lower bound convergence
    if (pseudo_mode_with_proof){
        if(verbose){
            cout << "(" << id <<") Proving correctness ..." << endl;
            cout << "(" << id <<") J for proof: " << endl;
//            printmapfromsets(J_for_proof);
        }
        calculate_range_from_J(J_for_proof, x_range, y_range, z_range, components);
        map< tuple<int, int, int, int>, int> tempspin;
        map< tuple<int, int, int, int, int, int>,
             map<set<tuple<int, int, int, int, int> >, double> > temp_clustertype_periodic;
        double temp_LB_so_far = lower_bound;
        get_lowerbound_sub(x_range, 0, y_range, 0, 0, z_range,
                           J_for_proof,
                           x_range, y_range, z_range, components,
                           tempspin,
                           exact_lower_bound,
                           temp_clustertype_periodic,
                           temp_LB_so_far, id,
                           false, false, true, obscenely_verbose,global_parameters);
    }

//    upper_bound = upper_bound * prec;
//    lower_bound = lower_bound * prec;
//    exact_lower_bound = exact_lower_bound * prec;
    upper_bound = upper_bound ;
    lower_bound = lower_bound ;
    exact_lower_bound = exact_lower_bound ;

    
    
    
    // Output results
    if (verbose){
        cout << "(" << id <<") ------------------- Results --------------------" << endl;
        if(basic_exact_mode){
            cout << "(" << id <<") UB: " << upper_bound << "; LB: " << lower_bound << endl;
        }else if(!pseudo_mode_with_proof){
            cout << "(" << id <<") UB: " << upper_bound << "; Pseudo-LB: " << lower_bound << endl;
        }else{
            cout << "(" << id <<") UB: " << upper_bound << "; Pseudo-LB: " << lower_bound << "; Exact-LB: " << exact_lower_bound << endl;
        }
        
        if(basic_exact_mode){
            cout << "(" << id <<") constant corrected UB: " << upper_bound+constant << "; LB: " << lower_bound +constant<< endl;
        }else if(!pseudo_mode_with_proof){
            cout << "(" << id <<") constant corrected UB: " << upper_bound+constant << "; Pseudo-LB: " << lower_bound+constant << endl;
        }else{
            cout << "(" << id <<") constant corrected UB: " << upper_bound+constant << "; Pseudo-LB: " << lower_bound+constant << "; Exact-LB: " << exact_lower_bound+constant << endl;
        }
        
        
        if (work_with_mu||global_parameters.ternary_alg) {
            double shifted_energy=mu_constant;
            for (map< set< tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
                shifted_energy+=upperboundclustertype[it1->first]*it1->second;
            }
            
            if(basic_exact_mode){
                cout << "(" << id <<") formation energy constant corrected UB: " << upper_bound+constant-shifted_energy << "; LB: " << lower_bound +constant-shifted_energy<< endl;
            }else if(!pseudo_mode_with_proof){
                cout << "(" << id <<") formation energy constant corrected UB: " << upper_bound+constant-shifted_energy << "; Pseudo-LB: " << lower_bound+constant-shifted_energy << endl;
            }else{
                cout << "(" << id <<") formation energy constant corrected UB: " << upper_bound+constant-shifted_energy << "; Pseudo-LB: " << lower_bound+constant-shifted_energy << "; Exact-LB: " << exact_lower_bound+constant-shifted_energy << endl;
            }
            
            
            if(basic_exact_mode){
                formation_energy_UB=upper_bound+constant-shifted_energy;
                formation_energy_LB=lower_bound+constant-shifted_energy;
            }else if(!pseudo_mode_with_proof){
                formation_energy_UB=upper_bound+constant-shifted_energy;
                formation_energy_LB=lower_bound+constant-shifted_energy;
            }else{
                formation_energy_UB=upper_bound+constant-shifted_energy;
                formation_energy_LB=exact_lower_bound+constant-shifted_energy;
            }


        }
        

        
        
        cout << "(" << id <<") Correlations: ";
        printmapfromsets(upperboundclustertype);
        cout << "\n(" << id <<") Unit cell: ";
        printblock(unitcell);
        cout << "\n(" << id <<") Periodicity: " << periodicity << endl;
        cout << "(" << id <<") --------------------- Done --------------------\n" << endl;
    }
    
    
        
        
    if (work_with_mu||global_parameters.ternary_alg) {
        double shifted_energy=mu_constant;
        for (map< set< tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
            shifted_energy+=upperboundclustertype[it1->first]*it1->second;
        }
                
        if(basic_exact_mode){
            formation_energy_UB=upper_bound+constant-shifted_energy;
            formation_energy_LB=lower_bound+constant-shifted_energy;
        }else if(!pseudo_mode_with_proof){
            formation_energy_UB=upper_bound+constant-shifted_energy;
            formation_energy_LB=lower_bound+constant-shifted_energy;
        }else{
            formation_energy_UB=upper_bound+constant-shifted_energy;
            formation_energy_LB=exact_lower_bound+constant-shifted_energy;
        }
        

    }
    
    if (input_PRIM_output_PRIMOUT_mode==1) {
        print_poscar_out(periodicity, unitcell, components);
    }
    

}

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
              double constant,double mu_constant)
{
    
    mpi::communicator mpi_world;
    
    global_parameters.formation_periodic_to_concentration_transfer.clear();
    global_parameters.formation_periodic_to_energy_transfer.clear();
    global_parameters.formation_mu_to_concentration_transfer.clear();
    global_parameters.formation_mu_to_energy_transfer.clear();
    global_parameters.map_mu_to_spin_transfer.clear();
    

    double max_value_of_J=0;
    for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1=J.begin(); it1!=J.end(); it1++) {
        if (fabs(it1->second)>max_value_of_J) {
            max_value_of_J=fabs(it1->second);
        }
    }
    
    map<tuple<int,int,int,int,int,int>, map<tuple<int,int,int,int>,int> > spin_periodic, spin_lowerbound;
    map<double, map<tuple<int,int,int,int>,int> > spin_mu;
    map<tuple<int,int,int,int,int,int>, double > energy_periodic,energy_lowerbound;
    map<double, double > energy_mu;
    

    int dimension;
    double lowerbound_from_compat = -1e18;
    double min_bound = -1e18, max_bound = 0;
    set<tuple<int,int,int,int,int,int> >  min_list, max_list;

    int a0, a1, a2, a3, a4, a5;
    set<tuple<int,int,int,int,int,int> > setoftupletemp;
    map<tuple<int,int,int,int,int,int>, map<set<tuple<int,int,int,int,int> >, double> > clustertype_periodic;
    map<tuple<int,int,int,int,int,int>, map<set<tuple<int,int,int,int,int> >, double> > clustertype_lowerbound;

    if (!((basic_exact_mode && !pseudo_mode && !pseudo_mode_with_proof) ||
           (!basic_exact_mode && pseudo_mode && pseudo_mode_with_proof) ||
           (!basic_exact_mode && pseudo_mode && !pseudo_mode_with_proof))){
        cout << "Invalid choice of solver mode! Exiting." << endl;
        exit(1);
    }

    // Identify dimensionality of the problem
    if (z_range == 1 && y_range == 1) {
        dimension = 1;
        if (very_verbose) { cout << "1D algorithm activated" << endl; }
    }else if (z_range == 1){
        dimension = 2;
        if (very_verbose) { cout << "2D algorithm activated" << endl; }
    }else{
        dimension = 3;
        if (very_verbose) { cout << "3D algorithm activated" << endl;
//            cout<<"\n what is global_parameters.dedicated1D: "<<global_parameters.dedicated1D<<" what is  dimension==1:"<<(dimension==1)<<endl;
 }
    }

    
//    cout<<"\n what is global_parameters.dedicated1D: "<<global_parameters.dedicated1D<<" what is  dimension==1:"<<(dimension==1)<<endl;
    
    //check if you are going to use dedicated 1D algorithm
    if (!(global_parameters.dedicated1D==true&&dimension==1)) {
        // Obtain lower bound estimate on the energy
        
//        cout<<"\n I am still alive at point: 4794sf216hg42136"<<endl;

        if (loopnumber > 0){
            if (!new_cluster_algorithm) {
                if (dimension == 1) {
                    a0=loopnumber;
                    a1=0;
                    a2=1;
                    a3=0;
                    a4=0;
                    a5=1;
                }else if (dimension == 2) {
                    a0=loopnumber;
                    a1=0;
                    a2=loopnumber;
                    a3=0;
                    a4=0;
                    a5=1;
                }else if (dimension == 3) {
                    a0=loopnumber;
                    a1=0;
                    a2=loopnumber;
                    a3=0;
                    a4=0;
                    a5=loopnumber;
                }
                
                get_lowerbound(a0, a1, a2, a3, a4, a5, J, x_range, y_range, z_range,
                               component, spin_lowerbound[make_tuple(a0,a1,a2,a3,a4,a5)],
                               energy_lowerbound[make_tuple(a0,a1,a2,a3,a4,a5)],
                               clustertype_lowerbound, J_for_proof, id,
                               pseudo_mode, pseudo_mode_with_proof, basic_exact_mode,
                               obscenely_verbose,limit_dimension, use_level_method, use_weighted_dual_average,global_parameters);
                
                
                
            }
            else {
                
                //for this one you could choose whether you incorporate range here
                set<set<tuple<int,int,int,int,int> > > J_defined_sets;
                set<tuple<int,int,int,int,int> > set_of_sites;
                int loopnumber_x=loopnumber;
                int loopnumber_y=loopnumber;
                int loopnumber_z=loopnumber;
                if (dimension==1) {
                    loopnumber_y=1;
                    loopnumber_z=1;
                    
                }
                else if (dimension==2){
                    loopnumber_z=1;
                }
                
                
                for (int i = 1; i <  x_range+loopnumber_x; i++) {
                    for (int j = 1; j <  y_range+loopnumber_y; j++) {
                        for (int k = 1; k <  z_range+loopnumber_z; k++) {
                            for (map<int, int>::iterator cit = component.begin();
                                 cit != component.end(); cit++) {
                                int position = cit->first;
                                int element_count = cit->second;
                                for (int mu = 1; mu < element_count; mu++) {
                                    
                                    set_of_sites.insert(make_tuple(i,j,k,position,mu));
                                    
                                }
                            }
                        }
                    }
                }
                
                
                for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1 = J.begin(); it1 != J.end(); it1++) {
                    
                    set<tuple<int,int,int,int,int> > prototype_set=it1->first;
                    set<tuple<int,int,int,int,int> > temp_equivalent_set;
                    convert_to_starting_point_cluster(prototype_set, temp_equivalent_set);
                    
                    
                    J_defined_sets.insert(temp_equivalent_set);
                    
                }
                //
                //            cout<<"\ncheck what is J_defined sets";
                //            printvectorofvector(J_defined_sets);
                
                map<set<tuple<int,int,int,int,int> >, double> J_new_with_zeros=J;
                set<set<tuple<int,int,int,int,int> > > J_new_with_zeros_defined_sets=J_defined_sets;
                
                
                if (use_new_pair_terms>0) {
                    int new_pair_terms_introduced=0;
                    for (set<tuple<int,int,int,int,int> >::iterator it1=set_of_sites.begin(); it1!=set_of_sites.end(); it1++) {
                        for (set<tuple<int,int,int,int,int> >::iterator it2=set_of_sites.begin(); it2!=set_of_sites.end(); it2++) {
                            set<tuple<int,int,int,int,int> > temp_set_of_sites;
                            temp_set_of_sites.insert(*it1);
                            temp_set_of_sites.insert(*it2);
                            set<tuple<int,int,int,int,int> > equivalent_set_of_sites;
                            convert_to_starting_point_cluster(temp_set_of_sites, equivalent_set_of_sites);
                            if (J_new_with_zeros_defined_sets.count(equivalent_set_of_sites)==0&&equivalent_set_of_sites.size()==2) {
                                J_new_with_zeros[equivalent_set_of_sites]+=0;
                                J_new_with_zeros_defined_sets.insert(equivalent_set_of_sites);
                                new_pair_terms_introduced++;
                                if (new_pair_terms_introduced>=use_new_pair_terms) {
                                    goto OUTSIDE_TWO_FOR_LOOPS_739159;
                                }
                            }
                        }
                    }
                    
                OUTSIDE_TWO_FOR_LOOPS_739159:
                    
                    1;
                }
                
                
                if (use_new_triplet_terms>0)
                {
                    int new_triplet_terms_introduced=0;
                    for (set<tuple<int,int,int,int,int> >::iterator it1=set_of_sites.begin(); it1!=set_of_sites.end(); it1++) {
                        for (set<tuple<int,int,int,int,int> >::iterator it2=set_of_sites.begin(); it2!=set_of_sites.end(); it2++) {
                            for (set<tuple<int,int,int,int,int> >::iterator it3=set_of_sites.begin(); it3!=set_of_sites.end(); it3++){
                                set<tuple<int,int,int,int,int> > temp_set_of_sites;
                                temp_set_of_sites.insert(*it1);
                                temp_set_of_sites.insert(*it2);
                                temp_set_of_sites.insert(*it3);
                                set<tuple<int,int,int,int,int> > equivalent_set_of_sites;
                                convert_to_starting_point_cluster(temp_set_of_sites, equivalent_set_of_sites);
                                if (J_new_with_zeros_defined_sets.count(equivalent_set_of_sites)==0&&equivalent_set_of_sites.size()==3) {
                                    J_new_with_zeros[equivalent_set_of_sites]+=0;
                                    J_new_with_zeros_defined_sets.insert(equivalent_set_of_sites);
                                    new_triplet_terms_introduced++;
                                    if (new_triplet_terms_introduced>=use_new_triplet_terms) {
                                        goto OUTSIDE_THREE_FOR_LOOPS_374932;
                                    }
                                }
                            }
                            
                        }
                    }
                    
                OUTSIDE_THREE_FOR_LOOPS_374932:
                    1;
                    
                }
                
                //            cout<<"\let's check the new J's";
                //            printvectorofvector(J_new_with_zeros_defined_sets);
                
                int x_range_new=x_range+loopnumber_x-1;
                int y_range_new=y_range+loopnumber_y-1;
                int z_range_new=z_range+loopnumber_z-1;
                
                a0=1;
                a1=0;
                a2=1;
                a3=0;
                a4=0;
                a5=1;
                
                get_lowerbound(a0, a1, a2, a3, a4, a5, J_new_with_zeros, x_range_new, y_range_new, z_range_new,
                               component, spin_lowerbound[make_tuple(a0,a1,a2,a3,a4,a5)],
                               energy_lowerbound[make_tuple(a0,a1,a2,a3,a4,a5)],
                               clustertype_lowerbound, J_for_proof, id,
                               pseudo_mode, pseudo_mode_with_proof, basic_exact_mode,
                               obscenely_verbose,limit_dimension, use_level_method, use_weighted_dual_average,global_parameters );
                
                
            }
            
            
            
            
            
            findminmax(energy_lowerbound, min_bound, max_bound);
            matchnumber(energy_lowerbound, max_bound, max_list);
            lowerbound_from_compat=max_bound;
            
            if(very_verbose){
                cout << "******************************************************************************" << endl;
                cout << "LB: " << endl;
                printmap(energy_lowerbound);
                
                cout << "LB solution: " << endl;
                for (set<tuple<int,int,int,int,int,int> >::iterator ait=max_list.begin(); ait!=max_list.end(); ait++){
                    tuple<int,int,int,int,int,int> a;
                    a=*ait;
                    cout << "\n" << a;
                    printblock(spin_lowerbound[a]);
                }
                
                cout << "\nLB maximum is: " << max_bound << endl;
                cout << "LB maximum list: ";
                printvector(max_list);
                cout << "\nAbsolute LB: " << lowerbound_from_compat << endl;
            }
        }
        
        // Obtain upper bound estimate on the energy
        if (max_sites > 0){
            
            vector<tuple<int,int,int,int,int,int> > periodicity_vector;
            periodicity_vector.clear();
            
//            cout<<"\n I am still alive at point: vn09n23biwa"<<endl;
            vector<double> mu_vector_numeric;
            if (global_parameters.activate_large_cell_algo_binary)
            {
                periodicity_vector.push_back(make_tuple(global_parameters.large_cell_binary_cubic_periodicity,0,global_parameters.large_cell_binary_cubic_periodicity,0,0,global_parameters.large_cell_binary_cubic_periodicity));
                
                
                for (int N=0 ; N< (global_parameters.large_cell_binary_mu_final_todo-global_parameters.large_cell_binary_mu_init_todo)/global_parameters.large_cell_binary_mu_step_todo+0.1; N++) {
                    double mu_now=global_parameters.large_cell_binary_mu_init_todo+global_parameters.large_cell_binary_mu_step_todo*N;
                    mu_vector_numeric.push_back(mu_now);
                }
                
                
            }
            else if (!global_parameters.activate_large_cell_algo_binary)
            {
                for (int alpha=1; alpha < max_sites+1; alpha++){
                    if (dimension == 1) {
                        a0=alpha;
                        a1=0;
                        a2=1;
                        a3=0;
                        a4=0;
                        a5=1;
                        
                        periodicity_vector.push_back(make_tuple(a0,a1,a2,a3,a4,a5));
                        
                    }else if (dimension == 2) {
                        a3=0;
                        a4=0;
                        a5=1;
                        for (int a2 = 1; a2 < alpha+1; a2++){
                            if (alpha % a2 == 0) {
                                a0 = alpha/a2;
                                for (int a1=0; a1<a0; a1++) {
                                    
                                    periodicity_vector.push_back(make_tuple(a0,a1,a2,a3,a4,a5));
                                    
                                    
                                }
                            }
                        }
                        
                        
                    }else if (dimension == 3){
                        
                        for (int beta = 1; beta < alpha + 1; beta++) {
                            if (alpha % beta == 0) {
                                a0 = alpha/beta;
                                for (int gamma = 1; gamma < beta+1; gamma++){
                                    if (beta % gamma == 0) {
                                        a2 = beta/gamma;
                                        a5 = gamma;
                                        
                                        if (global_parameters.Major_periodicity)
                                        {
                                            periodicity_vector.push_back(make_tuple(a0,0,a2,0,0,a5));
                                        }
                                        else{
                                            for (a4 = 0; a4 < a2; a4++) {
                                                for (a3 = 0; a3 < a0; a3++) {
                                                    
                                                    for (a1 = 0; a1 < a0; a1++) {
                                                        
                                                        periodicity_vector.push_back(make_tuple(a0,a1,a2,a3,a4,a5));
                                                        
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                    }
                                }
                            }
                            
                        }
                    }
                }
            }

            
            bool done=false;
            
//            cout<<"\n debug 02819374 periodicity_vector.size(): "<<periodicity_vector.size()<<endl;
            
            
            bool kill_all_slave=false;
//            usleep(100000);
//            cout<<"[MASTER] I am broadcasting kill_all_slave"<<endl;

            broadcast(mpi_world,kill_all_slave,0);

//            usleep(100000);
//            cout<<"[MASTER] I am broadcasting periodicity_vector"<<endl;

            broadcast(mpi_world,periodicity_vector,0);
            
            broadcast(mpi_world,mu_vector_numeric,0);


//            usleep(1000000);
            
            
//            cout<<"[MASTER] I am broadcasting J"<<endl;

            broadcast(mpi_world,J,0);
            broadcast(mpi_world,J_fixed_part,0);
            broadcast(mpi_world,constant,0);
            broadcast(mpi_world,mu_constant,0);

            

            usleep(1000000);
//            cout<<"[MASTER] I am broadcasting x_range"<<endl;

            broadcast(mpi_world,x_range,0);
            usleep(100000);
//            cout<<"[MASTER] I am broadcasting y_range"<<endl;

            broadcast(mpi_world,y_range,0);
            usleep(100000);
//            cout<<"[MASTER] I am broadcasting z_range"<<endl;

            broadcast(mpi_world,z_range,0);
            usleep(100000);
//            cout<<"[MASTER] I am broadcasting component"<<endl;
//            cout << "In run_solver x_range, y_range, z_range "<<make_tuple(x_range,y_range,z_range);

            broadcast(mpi_world,component,0);
            usleep(100000);
//            cout<<"[MASTER] I am broadcasting minbound"<<endl;

            broadcast(mpi_world,min_bound,0);
            usleep(100000);
//            lets write something here, for loop starts here:
            set<int> killed_slave;
            map<int, string > status_of_slave;
            
            // Initialize requests
            unsigned int job_id = 0;
            std::vector<mpi::request> reqs(mpi_world.size());
            
//            cout<<"\n hello I am root: I am going to send to receieve status status"<<endl;
            
            // Send initial jobs
            for (unsigned int dst_rank = 1; dst_rank < mpi_world.size(); ++dst_rank) {

                // Post receive request for new jobs requests by slave [nonblocking]
//                cout<<"[MASTER] I am i receving status_of_slave: "<<dst_rank<<endl;

                reqs[dst_rank] = mpi_world.irecv(dst_rank, 0,status_of_slave[dst_rank]);
//                ++job_id;
            }
            
//            cout<<"\n hello I am root: I finish status request"<<endl;

            
            // Send jobs as long as there is job left
            
            long long to_loop_over_object_size;
            if (global_parameters.activate_large_cell_algo_binary) {
                to_loop_over_object_size=mu_vector_numeric.size();
            }
            else if (!global_parameters.activate_large_cell_algo_binary){
                to_loop_over_object_size=periodicity_vector.size();
            }
                
                ;
            
            
            while(job_id < to_loop_over_object_size) {
                bool stop;
                for (unsigned int dst_rank = 1; dst_rank < mpi_world.size(); ++dst_rank) {
                    // Check if dst_rank is done
//                    cout<<"status_of_slave["<<dst_rank<<"] is "<<status_of_slave[dst_rank] <<endl;
                    if (reqs[dst_rank].test()) {
//                        cout<<"[MASTER] I receive reqs[dst_rank].test() is: "<<1<<endl;
//                        cout<<"status_of_slave["<<dst_rank<<"] is "<<status_of_slave[dst_rank] <<endl;

                        if(status_of_slave[dst_rank]=="finish") {
//note here, I deliberated exclude cluster_type_periodic which may be the culprit of a lot of bad computing performance for realistic system

//                            cout<<"[MASTER] I receive finish status "<<endl;

//                            remember to receive more detail output;
//                            and update the things you want to update;
                            map<tuple<int,int,int,int>,int>  spin_tmp;
                            double energy_tmp=0;
                            mpi_world.recv(dst_rank, 1,spin_tmp);
//                            cout<<"[MASTER] I receive spin_tmp status "<<endl;

                            mpi_world.recv(dst_rank, 2,energy_tmp);
//                            cout<<"[MASTER] I receive energy_tmp status "<<endl;

                            tuple<int,int,int,int,int,int>  periodicity_now;
                            
                            mpi_world.recv(dst_rank, 3,periodicity_now);
//                            cout<<"[MASTER] I receive periodicity_now status "<<endl;
                            
                            double mu_now;
                            if (global_parameters.activate_large_cell_algo_binary) {
                                mpi_world.recv(dst_rank, 6,mu_now);
                            }
                            
                            if (!global_parameters.ternary_alg) {
                                double formation_energy_out, concentration_out;
                                mpi_world.recv(dst_rank, 4,formation_energy_out);
//                                cout<<"[MASTER] I receive formation_energy_out status "<<endl;
                                
                                mpi_world.recv(dst_rank, 5,concentration_out);
//                                cout<<"[MASTER] I receive concentration_out status "<<endl;
                                
                                global_parameters.formation_periodic_to_concentration_transfer[periodicity_now]=concentration_out;
                                
                                global_parameters.formation_periodic_to_energy_transfer[periodicity_now]=formation_energy_out;
                                
                                if (global_parameters.activate_large_cell_algo_binary) {
                                    global_parameters.formation_mu_to_concentration_transfer[mu_now]=concentration_out;
                                    global_parameters.formation_mu_to_energy_transfer[mu_now]=formation_energy_out;
                                    
                                }
                            }
                            
                            
                            

                            
                            
                            
                            
                            
                            spin_periodic[periodicity_now]=spin_tmp;
                            energy_periodic[periodicity_now]=energy_tmp;
                            
                            if (global_parameters.activate_large_cell_algo_binary) {
                                spin_mu[mu_now]=spin_tmp;
                                energy_mu[mu_now]=energy_tmp;
                            }
                            
                            
                            bool activate_pruning=false; //I disable pruning since it is unstable in ising #109
                            if (activate_pruning)
                            {
                                double dump1;
                                findminmax(energy_periodic, min_bound, dump1);
                                matchnumber(energy_periodic, min_bound, min_list);
                                
                                if (very_verbose){
                                    cout << "\nUB: " << min_bound << ", LB: " << lowerbound_from_compat << "; Periodicity: ";
                                    printvector(min_list);
                                    cout << endl;
                                }
                                
                                if (min_bound <= lowerbound_from_compat+1e-5*max_value_of_J) {
                                    done=true;
                                }
                                
                            }
                            
                        }
                        
//                        std::cout << "[MASTER] Rank " << dst_rank << " is done.\n"<<endl;
                        // Check if there is remaining jobs
                        if (job_id  < to_loop_over_object_size) {
                            // Tell the slave that a new job is coming.
                            stop = false;
                            mpi_world.isend(dst_rank, 0, stop);
                            // Send the new job.
//                            std::cout << "[MASTER] Sending new job (" << job_id
//                            << ") to SLAVE " << dst_rank << ".\n"<<endl;
                            mpi_world.isend(dst_rank, 0, job_id);
                            mpi_world.isend(dst_rank, 1, min_bound);

                            reqs[dst_rank] = mpi_world.irecv(dst_rank, 0,status_of_slave[dst_rank]);
                            ++job_id;
                        }
                        else {
                            // Send stop message to slave.
                            stop = true;
                            mpi_world.send(dst_rank, 0, stop);
                            killed_slave.insert(dst_rank);
                        }
                    }
                }
                usleep(1000);
            }
            std::cout << "[MASTER] Sent all jobs.\n";
            
            // Listen for the remaining jobs, and send stop messages on completion.
            bool all_done = false;
            while (!all_done) {
                all_done = true;
                for (unsigned int dst_rank = 1; dst_rank < mpi_world.size(); ++dst_rank) {
                    
//                    cout<<"[MASTER] check if Slave "<<dst_rank<<" is killed?"<<endl;
                    
                    if(killed_slave.count(dst_rank)==0)
                    {
//                        cout<<"[MASTER] Slave "<<dst_rank<<" is not yet killed"<<endl;
//
//                        cout<<"[MASTER] last round; checking status message from slave "<<dst_rank<<endl;
                        
                        if (reqs[dst_rank].test()) {
                            
//                            cout<<"[MASTER] last round; checking status message from slave "<<dst_rank<<endl;
//                            cout<<"status_of_slave["<<dst_rank<<"] is "<<status_of_slave[dst_rank] <<endl;
                            
                            
                            if(status_of_slave[dst_rank]=="finish") {
                                //note here, I deliberated exclude cluster_type_periodic which may be the culprit of a lot of bad computing performance for realistic system
                                
//                                cout<<"[MASTER] I receive finish status from slave "<<dst_rank<<endl;
                                
                                //                            remember to receive more detail output;
                                //                            and update the things you want to update;
                                map<tuple<int,int,int,int>,int>  spin_tmp;
                                double energy_tmp=0;
                                mpi_world.recv(dst_rank, 1,spin_tmp);
//                                cout<<"[MASTER] I receive spin_tmp status "<<endl;
                                
                                mpi_world.recv(dst_rank, 2,energy_tmp);
//                                cout<<"[MASTER] I receive energy_tmp status "<<endl;
                                
                                tuple<int,int,int,int,int,int>  periodicity_now;
                                
                                mpi_world.recv(dst_rank, 3,periodicity_now);
//                                cout<<"[MASTER] I receive periodicity_now status "<<endl;
                                
                                double mu_now;
                                if (global_parameters.activate_large_cell_algo_binary) {
                                    mpi_world.recv(dst_rank, 6,mu_now);
                                }
                                
                                if (!global_parameters.ternary_alg) {
                                    double formation_energy_out, concentration_out;
                                    mpi_world.recv(dst_rank, 4,formation_energy_out);
//                                    cout<<"[MASTER] I receive formation_energy_out status "<<endl;
                                    
                                    mpi_world.recv(dst_rank, 5,concentration_out);
//                                    cout<<"[MASTER] I receive concentration_out status "<<endl;
                                    
                                    
                                    
                                    
                                    global_parameters.formation_periodic_to_concentration_transfer[periodicity_now]=concentration_out;
                                    
                                    global_parameters.formation_periodic_to_energy_transfer[periodicity_now]=formation_energy_out;
                                    
                                    if (global_parameters.activate_large_cell_algo_binary) {
                                        global_parameters.formation_mu_to_concentration_transfer[mu_now]=concentration_out;
                                        global_parameters.formation_mu_to_energy_transfer[mu_now]=formation_energy_out;
                                        
                                    }
                                }
                                
                                
                                
                                
                                
                                spin_periodic[periodicity_now]=spin_tmp;
                                energy_periodic[periodicity_now]=energy_tmp;
                                
                                
                                if (global_parameters.activate_large_cell_algo_binary) {
                                    spin_mu[mu_now]=spin_tmp;
                                    energy_mu[mu_now]=energy_tmp;
                                }
                                
                                
                                
//                                cout<<"[MASTER] I updated spin_periodic and energy_periodic "<<endl;

                                
                                bool activate_pruning=false; //I disable pruning since it is unstable in ising #109
                                if (activate_pruning) {
                                    double dump1;
                                    findminmax(energy_periodic, min_bound, dump1);
                                    matchnumber(energy_periodic, min_bound, min_list);
                                    
                                    if (very_verbose){
                                        cout << "\nUB: " << min_bound << ", LB: " << lowerbound_from_compat << "; Periodicity: ";
                                        printvector(min_list);
                                        cout << endl;
                                    }
                                    
                                    if (min_bound <= lowerbound_from_compat+1e-5*max_value_of_J) {
                                        done=true;
                                    }
                                    
                                }
                                
                                
//                                cout<<"[MASTER] I updated minlist and minbound "<<endl;

                                
                            }
                            // Tell the slave that it can exit.
                            bool stop = true;
                            mpi_world.send(dst_rank, 0, stop);
                            killed_slave.insert(dst_rank);
//                            cout<<"[MASTER] I stopped slave "<<dst_rank<<endl;

                        }
                        else {
                            all_done = false;
                        }
                    }
                    

                }
                usleep(1000);
            }
            std::cout << "[MASTER] Handled all jobs, killed every process.\n";
        }

        
        // Output final results
        
        cout<<"\n I am still alive at point: sjoqno2hcb02 "<<endl;

        


//        cout<<"\n I am still alive at point: jpxa0vyqbutdjvckqo"<<endl;

        
        int large_dump=1000000;
        
//        printmap(energy_periodic);
        
        double dump1;
        
        
        
        findminmax(energy_periodic, min_bound, dump1);
        matchnumber(energy_periodic, min_bound, min_list);
//         energy_periodic,energy_lowerbound;
//
        lower_bound = lowerbound_from_compat;
        upper_bound = min_bound;
        
        lowerboundclustertype = clustertype_lowerbound[*(max_list.begin())];

        
        
        for (set<tuple<int,int,int,int,int,int> >::iterator it1=min_list.begin(); it1!=min_list.end(); it1++) {
            tuple<int,int,int,int,int,int> periodicity_now=*it1;
            int a0=periodicity_now.get<0>();
            int a2=periodicity_now.get<2>();
            int a5=periodicity_now.get<5>();
            if (a0*a2*a5<large_dump) {
                periodicity=periodicity_now;
                large_dump=a0*a2*a5;
            }
        }
        
        
        
        
        
//        cout<<"\n I am still alive at point: bpahvfgwwweof"<<endl;

        
//        upperboundclustertype = clustertype_periodic[periodicity];
        unitcell = spin_periodic[periodicity];
        
        
        map<set<tuple<int, int, int, int, int> >, double> cluster_type_here;
        
        double upper_bound_temp=0;
        
//        cout<<"\n I am still alive at point: qwrwqsacbbi218"<<endl;

//        unsigned int microseconds=100000000;
//
//        usleep(microseconds);
        
//        cout<<"\n I am still alive at point: dasoihob2r281apdb25"<<endl;
        
//        cout<<"\nperiodicity is "<<periodicity;
//
//        cout<<endl;
        
        
        if (!global_parameters.activate_large_cell_algo_binary)
        {
            calculate_cluster_type_and_energy_periodic(J, periodicity, upper_bound_temp, cluster_type_here, unitcell);
        }
            
        
        
//        cout<<"\n I am still alive at point: 12iiv9310hidif"<<endl;

        upperboundclustertype = cluster_type_here;

        map_periodicity_to_spin=spin_periodic;
        
        global_parameters.map_mu_to_spin_transfer=spin_mu;

        
//        cout<<"\n I am still alive at point: vo009uidg98613ihk"<<endl;
        
        
    }
    else if(global_parameters.dedicated1D==true&&dimension==1){
        
        spin_periodic_struct spin_struct_min;
        dedicated1D_alg(J,x_range,component,spin_struct_min,id,global_parameters,lowerbound_from_compat,obscenely_verbose);
        
        //output
        lower_bound = lowerbound_from_compat;
        map<set<tuple<int, int, int, int, int> >, double> cluster_type_here;
        calculate_cluster_type_and_energy_periodic(J, spin_struct_min.periodicity, upper_bound, cluster_type_here, spin_struct_min.spin);
//        upper_bound = min_bound;
        
//        lowerboundclustertype = clustertype_lowerbound[*(max_list.begin())];
        
        upperboundclustertype = cluster_type_here;
        unitcell = spin_struct_min.spin;
        periodicity=spin_struct_min.periodicity;
//        map_periodicity_to_spin=spin_periodic;
        
    }
        
}





void simulated_annealing(map< set<tuple<int,int,int,int,int> >, double> J,solver_variable global_parameters,double &monte_energy, spin_periodic_struct &monte_spin_struct, double constant)
{
    cout<<"\n what is constant "<<constant;
    bool faster_alg=true,slow_alg=false;
    srand (time(NULL));
    map<int,int> components;
    int x_range, y_range, z_range;
    calculate_range_from_J(J, x_range, y_range, z_range, components);

    tuple<int,int,int,int,int,int> periodicity=global_parameters.monte_carlo_periodicity;
    int a0=global_parameters.monte_carlo_periodicity.get<0>();
    int a1=global_parameters.monte_carlo_periodicity.get<1>();
    int a2=global_parameters.monte_carlo_periodicity.get<2>();
    int a3=global_parameters.monte_carlo_periodicity.get<3>();
    int a4=global_parameters.monte_carlo_periodicity.get<4>();
    int a5=global_parameters.monte_carlo_periodicity.get<5>();
    
    
    
    map<tuple<int,int,int,int>,int> periodic_spin_now,periodic_spin_next,periodic_spin_best ;
    map<tuple<int,int,int,int>,int> periodic_spin_now_fast,periodic_spin_next_fast,periodic_spin_best_fast ;
//    random_initialization
    
    for (int i = 1; i < a0 + 1; i++) {
        for (int j = 1; j < a2 + 1; j++) {
            for (int k = 1; k < a5 + 1; k++) {
                for (map<int, int>::iterator cit = components.begin();
                     cit != components.end(); cit++) {
                    int position = cit->first;
                    int element_count = cit->second;
                    periodic_spin_now[make_tuple(i,j,k,position)]=rand()%element_count;
                    
                    
                }
            }
        }
    }
    periodic_spin_now_fast=periodic_spin_now;
    
//    end of random initialization

    double initial_E,next_E,best_E=1e30;
    double initial_E_fast,next_E_fast,best_E_fast=1e30;
    map<set<tuple<int,int,int,int,int> >, double> cluster_type_dump;
    calculate_cluster_type_and_energy_periodic(J, periodicity, initial_E, cluster_type_dump, periodic_spin_now);
    
    initial_E_fast=initial_E;
    
    periodic_spin_next=periodic_spin_now;
    periodic_spin_next_fast=periodic_spin_now;
    
    
    
    
    
    
    
// calculate
    
    long long grand_count=0;
    for (double T_now=global_parameters.initial_T; T_now>=global_parameters.final_T; T_now=T_now-global_parameters.T_decrement) {
        
//        cout<<"\nT_now is"<<T_now;
        
        for (long long loop_fixed_T=0; loop_fixed_T<global_parameters.each_T_how_many_runs; loop_fixed_T++) {
            grand_count++;
            int random_i=rand()%a0+1;
            int random_j=rand()%a2+1;
            int random_k=rand()%a5+1;
            int random_position=rand()%components.size()+1;
            int random_element=((rand()%(components[random_position]-1)+1)+periodic_spin_now[make_tuple(random_i,random_j,random_k,random_position)])%components[random_position];
            
            

            long long random_number_now=rand();
            
            if (slow_alg) {
                
                periodic_spin_next[make_tuple(random_i,random_j,random_k,random_position)]=random_element;
                
                cluster_type_dump.clear();
                calculate_cluster_type_and_energy_periodic(J, periodicity, next_E, cluster_type_dump, periodic_spin_next);
                double dE=next_E-initial_E;
                
//                cout<<"\n T_now: "<<T_now<<" loop_fixed_T: "<<loop_fixed_T;
//                cout<<"\n dE is :"<<dE ;
//                cout<<"\n T_now*8.6e-5 is :"<<T_now*8.6e-5<<endl;
//                cout<<"exp(-dE/(T_now*8.6e-5)) is :"<<exp(-dE/(T_now*8.6e-5))<<endl;
//                cout<<"(random_number_now%100)/100.0 is :"<<(random_number_now%100)/100.0<<endl;

//                cout<<"\n initial spin";
//                printblock(periodic_spin_now);
//                cout<<"\n next spin";
//                printblock(periodic_spin_next);
                
                if (exp(-dE/(T_now*8.6e-5))>(random_number_now%100)/100.0) {
                    periodic_spin_now=periodic_spin_next;
                    initial_E=next_E;
//                    cout<<"change accpeted"<<endl;
                }
                else{
                    periodic_spin_next=periodic_spin_now;
//                    cout<<"change rejected"<<endl;

                }
                if (next_E<best_E) {
                    best_E=next_E;
                    periodic_spin_best=periodic_spin_now;
                }
                
//                cout<<"\n best spin energy is"<<best_E;
//                cout<<"\nbest configuration is";
//                printblock(periodic_spin_best);

            }
            
            if (faster_alg) {
                periodic_spin_next_fast[make_tuple(random_i,random_j,random_k,random_position)]=random_element;
                double dE_faster=0;
                
                set<tuple<int,int,int,int> >changed_sites;
                
                for (int i = 1; i < a0 +x_range; i++) {
                    for (int j = 1; j < a2 + y_range;j++) {
                        for (int k = 1; k < a5 + z_range; k++){
                            
                            int i0 = positive_modulo(((i-1)-floor_int_division(j-1-(k-1)/a5*a4,a2)*a1-(k-1)/a5*a3), a0) + 1;
                            int j0 = positive_modulo(((j-1)-(k-1)/a5*a4), a2) + 1;
                            int k0 = ((k-1) % a5) + 1;
                            
                            if (i0==random_i&&j0==random_j&&k0==random_k) {
                                changed_sites.insert(make_tuple(i,j,k,random_position));
                            }
                        }
                    }
                    
                }
                
                
                for (int i1 = 1; i1 < a0 +1; i1++) {
                    for (int j1 = 1; j1 < a2 + 1;j1++) {
                        for (int k1 = 1; k1 < a5 + 1; k1++){
                            
                            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J.begin(); it1!=J.end(); it1++) {
                                set<tuple<int,int,int,int,int> > set_here= it1->first;
                                double J_value_now=it1->second;
                                bool involved_changed_sites=false;
                                for (set<tuple<int,int,int,int,int> >::iterator it2=set_here.begin(); it2!=set_here.end(); it2++){
                                    tuple<int,int,int,int,int> this_tuple=*it2;
                                    int i=this_tuple.get<0>()+i1-1;
                                    int j=this_tuple.get<1>()+j1-1;
                                    int k=this_tuple.get<2>()+k1-1;
                                    int p=this_tuple.get<3>();
                                    int i0 = positive_modulo(((i-1)-floor_int_division(j-1-(k-1)/a5*a4,a2)*a1-(k-1)/a5*a3), a0) + 1;
                                    int j0 = positive_modulo(((j-1)-(k-1)/a5*a4), a2) + 1;
                                    int k0 = ((k-1) % a5) + 1;
                                    
                                    if (changed_sites.count(make_tuple(i0,j0,k0,p))) {
                                        involved_changed_sites=true;
                                    }
                                    

                                }
                                
                                
                                if (involved_changed_sites==true) {
                                    int initial_satisfied=1,next_satisfied=1;
                                    for (set<tuple<int,int,int,int,int> >::iterator it2=set_here.begin(); it2!=set_here.end(); it2++){
                                        tuple<int,int,int,int,int> this_tuple=*it2;
                                        int i=this_tuple.get<0>()+i1-1;
                                        int j=this_tuple.get<1>()+j1-1;
                                        int k=this_tuple.get<2>()+k1-1;
                                        int p=this_tuple.get<3>();
                                        int c=this_tuple.get<4>();
                                        int i0 = positive_modulo(((i-1)-floor_int_division(j-1-(k-1)/a5*a4,a2)*a1-(k-1)/a5*a3), a0) + 1;
                                        int j0 = positive_modulo(((j-1)-(k-1)/a5*a4), a2) + 1;
                                        int k0 = ((k-1) % a5) + 1;

                                        int initial_match,next_match;
                                        if (periodic_spin_now_fast[make_tuple(i0,j0,k0,p)]==c) {
                                            initial_match=1;
                                        }
                                        else{
                                            initial_match=0;
                                        }
                                        
                                        if (periodic_spin_next_fast[make_tuple(i0,j0,k0,p)]==c) {
                                            next_match=1;
                                        }
                                        else{
                                            next_match=0;
                                        }
                                        initial_satisfied*=initial_match;
                                        next_satisfied*=next_match;
                                        
                                    }
                                    

                                    
                                    dE_faster+=(next_satisfied-initial_satisfied)*J_value_now/a0/a2/a5;

                                    
                                    
                                }
                            }
                        }
                    }
                    
                }
                
                //remember to cancel this for speed
//                bool I_am_just_debugging=true;
//                if (I_am_just_debugging) {
//
//
//                    cluster_type_dump.clear();
//                    calculate_cluster_type_and_energy_periodic(J, periodicity, next_E, cluster_type_dump, periodic_spin_next_fast);
//                    double dE=next_E-initial_E_fast;
//                    initial_E=next_E;
//                    cout<<"\n dE is "<<dE<<" and dE_faster is "<<dE_faster;
//
//                }
                
//                cout<<"\n fast alg\n dE_faster:"<<dE_faster ;
//                cout<<"\n initial spin";
//                printblock(periodic_spin_now_fast);
//                cout<<"\n next spin";
//                printblock(periodic_spin_next_fast);
                
                next_E_fast=dE_faster+initial_E_fast;
                
                if (exp(-dE_faster/(T_now*8.6e-5))>(random_number_now%100)/100.) {
//                    periodic_spin_now_fast=periodic_spin_next_fast;
                    periodic_spin_now_fast[make_tuple(random_i,random_j,random_k,random_position)]=random_element;
                    initial_E_fast=next_E_fast;
                }
                else{
                    periodic_spin_next_fast=periodic_spin_now_fast;
                }
                
                if (next_E_fast<best_E_fast-1e-8) {
                    best_E_fast=next_E_fast;
                    periodic_spin_best_fast=periodic_spin_now_fast;
                }
                
            }
        }
    }
    
    
    if (slow_alg) {
        best_E+=constant;
        monte_spin_struct.spin=periodic_spin_best;
        monte_spin_struct.periodicity=periodicity;
//        cout<<"\n use slow alg";
//        cout<<"\n best spin energy is"<<best_E;
//        cout<<"\nbest configuration is";
//        printblock(periodic_spin_best);
        
    }
    if (faster_alg) {
        best_E+=constant;
        monte_spin_struct.spin=periodic_spin_best;
        monte_spin_struct.periodicity=periodicity;

//        cout<<"\nuse fast alg";
//        cout<<"\n best spin energy is"<<best_E_fast;
//        cout<<"\nbest configuration is";
//        printblock(periodic_spin_best_fast);
    }

//



}






