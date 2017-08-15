#include "treeofdevil.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


using namespace std;
using namespace ::boost::tuples;
using namespace ::boost;

void get_lowerbound(int a0, int a1, int a2, int a3, int a4, int a5,
                    map<set<tuple<int,int,int,int,int> >, double> &J,
                    int x_range, int y_range, int z_range,
                    map<int,int> componentnumber,
                    map<tuple<int,int,int,int>,int> &spin,
                    double &lower_bound,
                    map<tuple<int,int,int,int,int,int>, map<set<tuple<int,int,int,int,int> >, double> > &clustertypeperiodic,
                    map<set<tuple<int,int,int,int,int> >, double> &J_for_proof,
                    std::string id,
                    bool pseudo_mode,
                    bool pseudo_mode_with_proof,
                    bool basic_exact_mode,
                    bool obscenely_verbose,
                    double limit_dimension, bool use_level_method,bool use_weighted_dual_average ,solver_variable &global_parameters){
    
    if (!((basic_exact_mode && !pseudo_mode && !pseudo_mode_with_proof) ||
          (!basic_exact_mode && pseudo_mode && pseudo_mode_with_proof) ||
          (!basic_exact_mode && pseudo_mode && !pseudo_mode_with_proof))){
        cout << "Invalid choice of solver mode! Exiting." << endl;
        exit(1);
    }
    
    if ((int)use_weighted_dual_average+(int)use_level_method>=2) {
        cout << "Invalid choice of lower bound method! Exiting." << endl;
        exit(1);
    }
    
    
    
    if (obscenely_verbose){
        cout << "---------------------- Computing lower bound ---------------------" << endl;
    }
    
    if (use_level_method||global_parameters.use_gradual_introduction_new_variables) {
        // Define optimization range (based on loopnumber and dimensionality)
        x_range += a0-1;
        y_range += a2-1;
        z_range += a5-1;
        
        // Set up Gurobi optimizer
        
        long long equivalent_sets_count = 0;
        map<long long, set<set<tuple<int,int,int,int,int> > > > index_to_equivalent_sets;
        map<set<set<tuple<int,int,int,int,int> > >,long long> revert_index_to_equivalent_sets;
        map<long long, double> value_of_equivalent_sets;
        double max_value_of_J = 0;
        
        for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1 = J.begin(); it1 != J.end(); it1++) {
            set<tuple<int,int,int,int,int> > prototype_set=it1->first;
            double value = it1->second;
            
            // Set some limits on x, y, z range
            int x_max = -1e5;
            int x_min = 1e5;
            int y_max = -1e5;
            int y_min = 1e5;
            int z_max = -1e5;
            int z_min = 1e5;
            for (set<tuple<int,int,int,int,int> >::iterator it2=prototype_set.begin(); it2!=prototype_set.end(); it2++) {
                tuple<int,int,int,int,int> tuple_element=*it2;
                int x_position=tuple_element.get<0>();
                int y_position=tuple_element.get<1>();
                int z_position=tuple_element.get<2>();
                if (x_min > x_position) { x_min = x_position; }
                if (x_max < x_position) { x_max = x_position; }
                if (y_min > y_position) { y_min = y_position; }
                if (y_max < y_position) { y_max = y_position; }
                if (z_min > z_position) { z_min = z_position; }
                if (z_max < z_position) { z_max = z_position; }
            }
            
            // Find equivalent sets of ECIs
            set<set<tuple<int,int,int,int,int> > > set_of_equivalent;
            for (int translation_z = 1-z_min; translation_z <= z_range-z_max; translation_z++) {
                for (int translation_y = 1-y_min; translation_y <= y_range-y_max; translation_y++) {
                    for (int translation_x = 1-x_min; translation_x <= x_range-x_max; translation_x++) {
                        set<tuple<int,int,int,int,int> > temp_equivalent_set;
                        
                        for (set<tuple<int,int,int,int,int> >::iterator it2=prototype_set.begin(); it2!=prototype_set.end(); it2++) {
                            tuple<int,int,int,int,int> tuple_element=*it2;
                            tuple<int,int,int,int,int> new_tuple=make_tuple(tuple_element.get<0>() +
                                                                            translation_x,tuple_element.get<1>() +
                                                                            translation_y,tuple_element.get<2>() +
                                                                            translation_z,tuple_element.get<3>(),
                                                                            tuple_element.get<4>());
                            temp_equivalent_set.insert(new_tuple);
                        }
                        set_of_equivalent.insert(temp_equivalent_set);
                    }
                }
            }
            
            if (revert_index_to_equivalent_sets.count(set_of_equivalent) == 0) {
                index_to_equivalent_sets[equivalent_sets_count] = set_of_equivalent;
                revert_index_to_equivalent_sets[set_of_equivalent] = equivalent_sets_count;
                value_of_equivalent_sets[equivalent_sets_count] = value;
                equivalent_sets_count++;
            } else {
                long long number = revert_index_to_equivalent_sets[set_of_equivalent];
                value_of_equivalent_sets[number] += value;
            }
        }
        
        for (map<long long, double>::iterator it1 = value_of_equivalent_sets.begin(); it1 != value_of_equivalent_sets.end(); it1++) {
            if (max_value_of_J < abs(it1->second)) { max_value_of_J = abs(it1->second); }
        }
        
        map<tuple<long long,long long>, set<tuple<int,int,int,int,int> > > index_to_subcluster;
        map<set<tuple<int,int,int,int,int> >, tuple<long long,long long> > revert_index_to_subcluster;
        map<tuple<long long,long long>, double> subcluster_value;
        map<long long,long long> prenum_to_subnum; //Last one does not count
        
        //WX 2015-03-17 begin
        
        long long count_dimension=0;
        vector<tuple<double,long long> > vector_pair_for_value_of_equivalent_sets;
        set<long long> significant_index;
        map<set<tuple<int,int,int,int,int> >, double> J_other_fixed;
        
        for (map<long long, double>::iterator it1=value_of_equivalent_sets.begin(); it1!=value_of_equivalent_sets.end();it1++ ) {
            vector_pair_for_value_of_equivalent_sets.push_back(make_tuple(fabs(it1->second),it1->first));
        }
        
        sort(vector_pair_for_value_of_equivalent_sets.rbegin(), vector_pair_for_value_of_equivalent_sets.rend());
        
        for (vector<tuple<double,long long> >::iterator it1=vector_pair_for_value_of_equivalent_sets.begin(); it1!=vector_pair_for_value_of_equivalent_sets.end(); it1++) {
            tuple<double,long long> tuple_now=*it1;
            long long index=tuple_now.get<1>();
            long long size=index_to_equivalent_sets[index].size();
            if (count_dimension<=limit_dimension) {
                significant_index.insert(index);
                count_dimension+=size;
            }
            else if (count_dimension>limit_dimension){
                break;
            }
        }
        
        //        a key difference here
        //        significant_index.clear();
        //note use significant_index and index_to_equivalent_sets to activate the progra
        
        
        //        initial set up done;
        //        loop over something;
        //            solves the predefined convex optimization problem;
        //            define the new convex model;
        //            check if you would like to stop;
        
        vector<map<tuple<int,int,int,int>,int> > spin_collections_vector;
        set<spin_no_periodic_struct>  warm_start_spin_structs ;
        bool warm_restart=false;
        map<set<tuple<int,int,int,int,int> >, double> warm_start_J;
        
        
        int model_number_max=global_parameters.model_numbers;
        if (use_level_method) {
            model_number_max=1;
        };
        
        
        for (long long model_update_loops=0;model_update_loops<model_number_max ; model_update_loops++) {
            
            //solve the convex problem
            //note use significant_index and index_to_equivalent_sets to activate the progra
            
            spin_collections_vector.clear();
            map<set<tuple<int,int,int,int,int> >, double> J_input_best_memory;
            map<spin_no_periodic_struct, double> PI_value;
            
            
            
            {
                index_to_subcluster.clear();
                revert_index_to_subcluster.clear();
                prenum_to_subnum.clear();
                J_other_fixed.clear();
                
                map<spin_no_periodic_struct, GRBConstr> spin_struct_to_Constr;
                
                GRBEnv env = GRBEnv();
                env.set(GRB_IntParam_Method, 1);
                env.set(GRB_IntParam_ScaleFlag, 0);
                env.set(GRB_DoubleParam_ObjScale, 1);
                env.set(GRB_IntParam_DualReductions, 0);
                env.set(GRB_IntParam_Presolve, 0);
                
                GRBEnv env_quadratic = GRBEnv();
                env_quadratic.set(GRB_IntParam_DualReductions, 0);
                env_quadratic.set(GRB_IntParam_Method, 0);
                env_quadratic.set(GRB_IntParam_Quad, 0);
                env_quadratic.set(GRB_DoubleParam_TimeLimit, 10);
                
                GRBModel m = GRBModel(env);
                GRBModel m_quadratic = GRBModel(env_quadratic);
                if (!obscenely_verbose){
                    m.getEnv().set(GRB_IntParam_OutputFlag, 0);
                    m_quadratic.getEnv().set(GRB_IntParam_OutputFlag, 0);
                }
                GRBVar mu = m.addVar(-1e25, 1e25, 1, GRB_CONTINUOUS);
                GRBVar mu_quadratic = m_quadratic.addVar(-1e17, 1e17, 0, GRB_CONTINUOUS);
                map<tuple<long long,long long>, GRBVar> subcluster_variable,subcluster_variable_quadratic;
                
                m.update();
                m_quadratic.update();
                m.set(GRB_IntAttr_ModelSense,-1);
                m_quadratic.set(GRB_IntAttr_ModelSense,-1);
                
                for (map<long long, set<set<tuple<int,int,int,int,int> > > >::iterator it1 = index_to_equivalent_sets.begin(); it1 != index_to_equivalent_sets.end(); it1++) {
                    if (significant_index.count(it1->first)) {
                        long long sub_number = 0;
                        long long pre_number = it1->first;
                        set<set<tuple<int,int,int,int,int> > > equivalent_set_here = it1->second;
                        for (set<set<tuple<int,int,int,int,int> > >::iterator it2=equivalent_set_here.begin(); it2 != equivalent_set_here.end(); it2++) {
                            set<tuple<int,int,int,int,int> > subset_here = *it2;
                            index_to_subcluster[make_tuple(pre_number,sub_number)] = subset_here;
                            revert_index_to_subcluster[subset_here]=make_tuple(pre_number, sub_number);
                            subcluster_variable[make_tuple(pre_number,sub_number)] = m.addVar(-max_value_of_J*500, max_value_of_J*500, 0, GRB_CONTINUOUS);
                            subcluster_variable_quadratic[make_tuple(pre_number,sub_number)] = m_quadratic.addVar(-max_value_of_J*500, max_value_of_J*500,0, GRB_CONTINUOUS);
                            sub_number++;
                        }
                        prenum_to_subnum[pre_number]=sub_number;
                    }
                    else {
                        long long pre_number = it1->first;
                        set<set<tuple<int,int,int,int,int> > > equivalent_set_here = it1->second;
                        for (set<set<tuple<int,int,int,int,int> > >::iterator it2=equivalent_set_here.begin(); it2 != equivalent_set_here.end(); it2++) {
                            set<tuple<int,int,int,int,int> > subset_here = *it2;
                            if (warm_restart==true) {
                                J_other_fixed[subset_here]=warm_start_J[subset_here];
                            }
                            else{
                                J_other_fixed[subset_here]=value_of_equivalent_sets[pre_number]/index_to_equivalent_sets[pre_number].size();
                            }
                        }
                    }
                }
                
                
                
                m.update();
                m_quadratic.update();
                
                
                
                
                for (map<long long, long long>::iterator it1=prenum_to_subnum.begin(); it1!=prenum_to_subnum.end(); it1++) {
                    long long prenumber=it1->first;
                    GRBLinExpr sum = 0;
                    GRBLinExpr sum_quadratic = 0;
                    for (long long subnumber = 0; subnumber < prenum_to_subnum[prenumber]; subnumber++) {
                        sum += subcluster_variable[make_tuple(prenumber,subnumber)];
                        assert(subcluster_variable_quadratic.count(make_tuple(prenumber,subnumber))==1);
                        sum_quadratic += subcluster_variable_quadratic[make_tuple(prenumber,subnumber)];
                    }
                    m.addConstr(sum == value_of_equivalent_sets[prenumber]);
                    m_quadratic.addConstr(sum_quadratic == value_of_equivalent_sets[prenumber]);
                }
                
                
                
                map<set<tuple<int,int,int,int,int> >, double> J_input;
                //                double proximal_factor=1e-8;
                //                GRBQuadExpr quadratic_objective = mu_quadratic;
                GRBQuadExpr quadratic_objective = 0;
                for (map<long long, long long>::iterator it1=prenum_to_subnum.begin(); it1!=prenum_to_subnum.end(); it1++) {
                    long long prenumber=it1->first;
                    for (long long subnumber = 0; subnumber < prenum_to_subnum[prenumber]; subnumber++) {
                        //                        double prefactor = proximal_factor;
                        if (warm_restart==true) {
                            quadratic_objective += (-1)*(subcluster_variable_quadratic[make_tuple(prenumber,subnumber)]-warm_start_J[index_to_subcluster[make_tuple(prenumber,subnumber)]])*(subcluster_variable_quadratic[make_tuple(prenumber,subnumber)]-warm_start_J[index_to_subcluster[make_tuple(prenumber,subnumber)]]);
                        }
                        else
                        {
                            quadratic_objective += (-1)*subcluster_variable_quadratic[make_tuple(prenumber,subnumber)]*
                            subcluster_variable_quadratic[make_tuple(prenumber,subnumber)];
                        }
                    }
                }
                
                m_quadratic.setObjective(quadratic_objective*1e-10);
                double lowerbound_memory=-1e50;
                
                if (warm_restart==true) {
                    
                    //                    cout<<"\n let's check what is warm_start_J";
                    //                    printmapfromsets(warm_start_J);
                    
                    for (set<spin_no_periodic_struct>::iterator  it1=warm_start_spin_structs.begin() ; it1!=warm_start_spin_structs.end(); it1++) {
                        spin_no_periodic_struct spin_struct_now=*it1;
                        map<set<tuple<int,int,int,int,int> >, double> cluster_type_here;
                        double energy;
                        
                        
                        calculate_cluster_type_and_energy_no_periodic(warm_start_J, energy, cluster_type_here, spin_struct_now.spin);
                        
                        GRBLinExpr subcluster_sum = 0;
                        
                        
                        
                        for (map<set<tuple<int,int,int,int,int> >, double>::iterator it2 = cluster_type_here.begin(); it2 != cluster_type_here.end(); it2++) {
                            if (revert_index_to_subcluster.count(it2->first)==1) {
                                subcluster_sum += subcluster_variable[revert_index_to_subcluster[it2->first]]*it2->second;
                                //                                cout<<"\n debug #572164882 ";
                                //                                cout<<"\n let's firstly check what spin_struct_now is here";
                                //                                printblock(spin_struct_now.spin);
                                //                                cout<<"\n then check what is it2->first";
                                //                                printvector(it2->first);
                                //                                cout<<"\n then check what is it2->second";
                                //                                cout<<"\n "<<it2->second<<endl;
                            }
                            else{
                                subcluster_sum += J_other_fixed[it2->first]*it2->second;
                            }
                        }
                        
                        spin_struct_to_Constr[spin_struct_now]=m.addConstr(mu <= subcluster_sum);
                        
                        
                        
                        GRBLinExpr sum1_quadratic = 0;
                        
                        //                        cout<<"\n debug 104 here what is cluster type:"<<endl;
                        //                        printmapfromsets(cluster_type_here);
                        for (map<set<tuple<int,int,int,int,int> >, double>::iterator it2 = cluster_type_here.begin(); it2 != cluster_type_here.end(); it2++) {
                            if (revert_index_to_subcluster.count(it2->first)==1) {
                                sum1_quadratic += subcluster_variable_quadratic[revert_index_to_subcluster[it2->first]]*it2->second;
                            }
                            else {
                                sum1_quadratic+=J_other_fixed[it2->first]*it2->second;
                            }
                        }
                        m_quadratic.addConstr(mu_quadratic <= sum1_quadratic);
                        
                        
                        spin_collections_vector.push_back(spin_struct_now.spin);
                    }
                    
                    if (obscenely_verbose)
                        cout<<"\n debug here check if upper bound of lowerbound is correct";
                    m.optimize();
                    
                }
                
                
                
                GRBConstr lambda_constraint;
                bool decide_whether_to_abandon_QP=false;
                for (long long infinite_loop = 0; infinite_loop < 500000; infinite_loop++) {
                    
                    //                    cout<<"\n\n debug here do quadratic programming optimization\n"<<endl;
                    m_quadratic.optimize();
                    J_input.clear();
                    J_input=J_other_fixed;
                    
                    for (map<tuple<long long,long long>, GRBVar>::iterator it1 = subcluster_variable_quadratic.begin(); it1 != subcluster_variable_quadratic.end(); it1++) {
                        
                        tuple<long long,long long> tuple_here = it1->first;
                        GRBVar variable_here=it1->second;
                        try {
                            J_input[index_to_subcluster[tuple_here]] = variable_here.get(GRB_DoubleAttr_X);
                        } catch (GRBException a) {
                            lower_bound = lowerbound_memory;
                            J_for_proof = J_input_best_memory;
                            if (obscenely_verbose){
                                cout << "LB - Loop #" << infinite_loop <<": Either a numerical error occurred or found tightest LB: " << lower_bound << endl;
                                //                    cout << "J_for_proof: ";
                                //                    printmapfromsets(J_for_proof);
                                cout << "\n" << endl;
                            }
                            decide_whether_to_abandon_QP=true;
                            if (obscenely_verbose==true) {
                                cout<<"\n we decide to abandon QP now \n";
                            }
                            goto OUTSIDE_THE_LOOP_NUM9726372_1;
                        }
                        
                        
                    }
                    
                    if (pseudo_mode) {
                        get_lowerbound_sub(1, 0, 1, 0, 0, 1,
                                           J_input,
                                           x_range, y_range, z_range, componentnumber,
                                           spin, lower_bound,
                                           clustertypeperiodic, lowerbound_memory, id,
                                           true, false, false,
                                           obscenely_verbose, global_parameters);
                    }else if (basic_exact_mode) {
                        get_lowerbound_sub(1, 0, 1, 0, 0, 1,
                                           J_input,
                                           x_range, y_range, z_range, componentnumber,
                                           spin, lower_bound,
                                           clustertypeperiodic, lowerbound_memory, id,
                                           false, true, false,
                                           obscenely_verbose, global_parameters);
                    }
                    
                    
                    
                    //                    check if spin already in spin_collections_vector, if not add it there and add the corresponding constraints;
                    map<tuple<int,int,int,int>,int> spin_lower_bound=spin;
                    bool inside_indicator=false;
                    for (int i=0; i<spin_collections_vector.size(); i++) {
                        map<tuple<int,int,int,int>,int> spin_now= spin_collections_vector[i];
                        if (spin_now==spin_lower_bound) {
                            inside_indicator=true;
                            break;
                        }
                    }
                    
                    //                    cout<<"\n inside indicator true or false: "<<inside_indicator<<"\n";
                    
                    if (inside_indicator==false) {
                        spin_collections_vector.push_back(spin_lower_bound);
                    }
                    
                    //                    end of check if spin already in spin_collections_vector, if not add it there and add the corresponding constraints;
                    
                    map<set<tuple<int,int,int,int,int> >, double> cluster_type_now = clustertypeperiodic[make_tuple(1,0,1,0,0,1)];
                    
                    if (lowerbound_memory <= lower_bound) {
                        lowerbound_memory = lower_bound;
                        J_input_best_memory = J_input;
                        J_for_proof = J_input_best_memory;
                    }
                    
                    
                    //                    cout<<"\n debug here #879612164799965 check if inside indicator is active, inside_indicator="<<inside_indicator<<endl;
                    
                    
                    if (inside_indicator==false) {
                        
                        GRBLinExpr subcluster_sum = 0;
                        for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1 = cluster_type_now.begin(); it1 != cluster_type_now.end(); it1++) {
                            if (revert_index_to_subcluster.count(it1->first)==1) {
                                subcluster_sum += subcluster_variable[revert_index_to_subcluster[it1->first]]*it1->second;
                            }
                            else{
                                subcluster_sum += J_other_fixed[it1->first]*it1->second;
                            }
                        }
                        spin_no_periodic_struct spin_struct_temp;
                        spin_struct_temp.spin=spin_lower_bound;
                        spin_struct_to_Constr[spin_struct_temp]=m.addConstr(mu <= subcluster_sum);
                        
                        GRBLinExpr sum1_quadratic = 0;
                        for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1 = cluster_type_now.begin(); it1 != cluster_type_now.end(); it1++) {
                            if (revert_index_to_subcluster.count(it1->first)==1) {
                                sum1_quadratic += subcluster_variable_quadratic[revert_index_to_subcluster[it1->first]]*it1->second;
                            }
                            else {
                                sum1_quadratic+=J_other_fixed[it1->first]*it1->second;
                            }
                        }
                        m_quadratic.addConstr(mu_quadratic <= sum1_quadratic);
                        
                    }
                    
                    GRBLinExpr  m_objective=mu;
                    m.setObjective(m_objective);
                    m.optimize();
                    double upperbound_of_lowerbound = m.get(GRB_DoubleAttr_ObjVal);
                    
                    if (infinite_loop >= 1) {
                        m_quadratic.remove(lambda_constraint);
                    }
                    lambda_constraint = m_quadratic.addConstr(mu_quadratic >= (0.5*lowerbound_memory + 0.5*upperbound_of_lowerbound));
                    
                    
                    
                    
                    quadratic_objective = 0;
                    
                    for (map<long long, long long>::iterator it1=prenum_to_subnum.begin(); it1!=prenum_to_subnum.end(); it1++) {
                        long long prenumber=it1->first;
                        for (long long subnumber = 0; subnumber < prenum_to_subnum[prenumber]; subnumber++) {
                            quadratic_objective -= (subcluster_variable_quadratic[make_tuple(prenumber,subnumber)] -
                                                    J_input_best_memory[index_to_subcluster[make_tuple(prenumber,subnumber)]]) *
                            (subcluster_variable_quadratic[make_tuple(prenumber,subnumber)] -
                             J_input_best_memory[index_to_subcluster[make_tuple(prenumber,subnumber)]]);
                            
                            //                            subcluster_variable_quadratic[make_tuple(prenumber,subnumber)].set(GRB_DoubleAttr_Start,J_input_best_memory[index_to_subcluster[make_tuple(prenumber,subnumber)]]);
                            
                            subcluster_variable_quadratic[make_tuple(prenumber,subnumber)].set(GRB_DoubleAttr_Start,subcluster_variable[make_tuple(prenumber,subnumber)].get(GRB_DoubleAttr_X));
                        }
                    }
                    
                    m_quadratic.setObjective(quadratic_objective*1e-10);
                    
                    if (obscenely_verbose){
                        cout << "\nModel Number: "<<model_update_loops;
                        cout << "\nLB - Loop #" << infinite_loop << endl;
                        cout << "\tBound of LB: " << upperbound_of_lowerbound << endl;
                        cout << "\tCurrent LB: " << lower_bound << endl;
                        cout << "\tBest LB: " << lowerbound_memory<< endl;
                        cout << "\n" << endl;
                    }
                    
                    //                    assert(lower_bound<upperbound_of_lowerbound+max_value_of_J*1e-4);
                    if (upperbound_of_lowerbound+max_value_of_J*1e-6<lowerbound_memory) {
                        if (obscenely_verbose) {
                            cout<<"\n I am coming to a very weird scenario that upperbound_of_lowerbound+max_value_of_J*1e-6<lowerbound_memory";
                        }
                        lowerbound_memory = lower_bound;
                        J_input_best_memory = J_input;
                        J_for_proof = J_input_best_memory;
                    }
                    
                    
                    
                    if ((fabs(lowerbound_memory-upperbound_of_lowerbound)/fabs(max_value_of_J)<1e-3)||lower_bound>upperbound_of_lowerbound+max_value_of_J*1e-4) {
                        J_for_proof = J_input_best_memory;
                        lower_bound = lowerbound_memory;
                        if (obscenely_verbose){
                            cout << "\n Model number: "<< model_update_loops <<" LB - Loop with QP #" << infinite_loop << ": Found tightest lower bound: " << lower_bound << endl;
                            cout << "J_for_proof:";
                            //                            printmapfromsets(J_for_proof);
                            cout << "\n" << endl;
                            
                        }
                        
                        //                        for (map<spin_no_periodic_struct, GRBConstr>::iterator it3=spin_struct_to_Constr.begin(); it3!=spin_struct_to_Constr.end(); it3++) {
                        //                            GRBConstr constraint_now=it3->second;
                        //                            spin_no_periodic_struct spin_struct_now=it3->first;
                        //
                        //                            PI_value[it3->first]=constraint_now.get(GRB_DoubleAttr_Pi);
                        ////                            cout<<"\n debug here #78794613 calculate PI value ";
                        ////                            printblock(spin_struct_now.spin);
                        ////                            cout<<"\n PI value is: "<<PI_value[it3->first]<<endl;
                        //                        }
                        
                        decide_whether_to_abandon_QP=true;
                        if (obscenely_verbose==true) {
                            cout<<"\n we decide to abandon QP now \n";
                        }
                        
                        
                        
                        break;
                    }
                    
                    
                    
                }
            OUTSIDE_THE_LOOP_NUM9726372_1:
                if (decide_whether_to_abandon_QP==true)
                    for (long long infinite_loop_without_QP = 0; infinite_loop_without_QP < 500000; infinite_loop_without_QP++)
                    {
                        
                        //                    cout<<"\n\n debug here do quadratic programming optimization\n"<<endl;
                        // start without QP
                        if (obscenely_verbose) {
                            cout<<"\nnow we abandon QP!\n";
                        }
                        m.optimize();
                        J_input.clear();
                        J_input=J_other_fixed;
                        
                        for (map<tuple<long long,long long>, GRBVar>::iterator it1 = subcluster_variable.begin(); it1 != subcluster_variable.end(); it1++) {
                            
                            tuple<long long,long long> tuple_here = it1->first;
                            GRBVar variable_here=it1->second;
                            J_input[index_to_subcluster[tuple_here]] = variable_here.get(GRB_DoubleAttr_X);
                            
                            
                            
                        }
                        
                        if (pseudo_mode) {
                            get_lowerbound_sub(1, 0, 1, 0, 0, 1,
                                               J_input,
                                               x_range, y_range, z_range, componentnumber,
                                               spin, lower_bound,
                                               clustertypeperiodic, lowerbound_memory, id,
                                               true, false, false,
                                               obscenely_verbose, global_parameters);
                        }else if (basic_exact_mode) {
                            get_lowerbound_sub(1, 0, 1, 0, 0, 1,
                                               J_input,
                                               x_range, y_range, z_range, componentnumber,
                                               spin, lower_bound,
                                               clustertypeperiodic, lowerbound_memory, id,
                                               false, true, false,
                                               obscenely_verbose, global_parameters);
                        }
                        
                        
                        
                        //                    check if spin already in spin_collections_vector, if not add it there and add the corresponding constraints;
                        map<tuple<int,int,int,int>,int> spin_lower_bound=spin;
                        bool inside_indicator=false;
                        for (int i=0; i<spin_collections_vector.size(); i++) {
                            map<tuple<int,int,int,int>,int> spin_now= spin_collections_vector[i];
                            if (spin_now==spin_lower_bound) {
                                inside_indicator=true;
                                break;
                            }
                        }
                        
                        //                    cout<<"\n inside indicator true or false: "<<inside_indicator<<"\n";
                        
                        if (inside_indicator==false) {
                            spin_collections_vector.push_back(spin_lower_bound);
                        }
                        
                        //                    end of check if spin already in spin_collections_vector, if not add it there and add the corresponding constraints;
                        
                        map<set<tuple<int,int,int,int,int> >, double> cluster_type_now = clustertypeperiodic[make_tuple(1,0,1,0,0,1)];
                        
                        if (lowerbound_memory <= lower_bound) {
                            lowerbound_memory = lower_bound;
                            J_input_best_memory = J_input;
                            J_for_proof = J_input_best_memory;
                        }
                        

                        
                        
                        
                        
                        //                    cout<<"\n debug here #879612164799965 check if inside indicator is active, inside_indicator="<<inside_indicator<<endl;
                        
                        
                        if (inside_indicator==false) {
                            
                            GRBLinExpr subcluster_sum = 0;
                            for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1 = cluster_type_now.begin(); it1 != cluster_type_now.end(); it1++) {
                                if (revert_index_to_subcluster.count(it1->first)==1) {
                                    subcluster_sum += subcluster_variable[revert_index_to_subcluster[it1->first]]*it1->second;
                                }
                                else{
                                    subcluster_sum += J_other_fixed[it1->first]*it1->second;
                                }
                            }
                            spin_no_periodic_struct spin_struct_temp;
                            spin_struct_temp.spin=spin_lower_bound;
                            spin_struct_to_Constr[spin_struct_temp]=m.addConstr(mu <= subcluster_sum);
                            
                            GRBLinExpr sum1_quadratic = 0;
                            for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1 = cluster_type_now.begin(); it1 != cluster_type_now.end(); it1++) {
                                if (revert_index_to_subcluster.count(it1->first)==1) {
                                    sum1_quadratic += subcluster_variable_quadratic[revert_index_to_subcluster[it1->first]]*it1->second;
                                }
                                else {
                                    sum1_quadratic+=J_other_fixed[it1->first]*it1->second;
                                }
                            }
                            m_quadratic.addConstr(mu_quadratic <= sum1_quadratic);
                            
                        }
                        
                        GRBLinExpr  m_objective=mu;
                        m.setObjective(m_objective);
                        m.optimize();
                        double upperbound_of_lowerbound = m.get(GRB_DoubleAttr_ObjVal);
                        
                        
                        
                        if (obscenely_verbose){
                            cout << "\nModel Number: "<<model_update_loops;
                            cout << "\nLB without QP - Loop #" <<  infinite_loop_without_QP << endl;
                            cout << "\tBound of LB: " << upperbound_of_lowerbound << endl;
                            cout << "\tCurrent LB: " << lower_bound << endl;
                            cout << "\tBest LB: " << lowerbound_memory<< endl;
                            cout << "\n" << endl;
                        }
                        
                        //                        assert(lower_bound<upperbound_of_lowerbound+max_value_of_J*1e-4);
                        
                        if (upperbound_of_lowerbound+max_value_of_J*1e-6<lowerbound_memory) {
                            if (obscenely_verbose) {
                                cout<<"\n I am coming to a very weird scenario that upperbound_of_lowerbound+max_value_of_J*1e-6<lowerbound_memory";
                            }
                            lowerbound_memory = lower_bound;
                            J_input_best_memory = J_input;
                            J_for_proof = J_input_best_memory;
                        }
                        
                        
                        
                        if ((fabs(lowerbound_memory-upperbound_of_lowerbound)/fabs(max_value_of_J)<1e-7)||lower_bound>upperbound_of_lowerbound-max_value_of_J*1e-7) {
                            J_for_proof = J_input_best_memory;
                            lower_bound = lowerbound_memory;
                            if (obscenely_verbose){
                                cout << "\n Model number: "<< model_update_loops <<" LB without QP - Loop #" << infinite_loop_without_QP << ": Found tightest lower bound: " << lower_bound << endl;
                                cout << "J_for_proof:";
                                //                            printmapfromsets(J_for_proof);
                                cout << "\n" << endl;
                                
                            }
                            
                            for (map<spin_no_periodic_struct, GRBConstr>::iterator it3=spin_struct_to_Constr.begin(); it3!=spin_struct_to_Constr.end(); it3++) {
                                GRBConstr constraint_now=it3->second;
                                spin_no_periodic_struct spin_struct_now=it3->first;
                                
                                PI_value[it3->first]=constraint_now.get(GRB_DoubleAttr_Pi);
                                //                            cout<<"\n debug here #78794613 calculate PI value ";
                                //                            printblock(spin_struct_now.spin);
                                //                            cout<<"\n PI value is: "<<PI_value[it3->first]<<endl;
                            }
                            
                            
                            break;
                        }
                        
                    }
                1;
                
            }
            
            
            //            build up lowest_energy_spin_struct_set
            double min=1e30;
            set<spin_no_periodic_struct> lowest_energy_spin_struct_set;
            {
                double min_energy=3e50;
                vector<double> energy_vector;
                
                
                for (int i1=0; i1<spin_collections_vector.size(); i1++) {
                    map<tuple<int,int,int,int>,int> spin_now=spin_collections_vector[i1];
                    map<set<tuple<int,int,int,int,int> >, double> cluster_type_here;
                    double energy;
                    //                calculate_cluster_type_and_energy_without_periodicity;
                    calculate_cluster_type_and_energy_no_periodic(J_for_proof, energy, cluster_type_here, spin_now);
                    
                    energy_vector.push_back(energy);
                }
                
                
                min=1e30;
                for (int i1=0; i1<energy_vector.size(); i1++) {
                    if (min>energy_vector[i1]) {
                        min=energy_vector[i1];
                    }
                }
                
                
                //            cout<<"\ndebug here, what is min: "<<min;
                vector<int> min_list;
                for (int i1=0; i1<energy_vector.size(); i1++) {
                    //                assert(min<1);
                    if (energy_vector[i1]<min+1e-5*max_value_of_J) {
                        min_list.push_back(i1);
                    }
                }
                
                
                
                for (int i1=0; i1<min_list.size(); i1++) {
                    map<tuple<int,int,int,int>,int> spin_temp=spin_collections_vector[min_list[i1]];
                    spin_no_periodic_struct spin_no_periodic_struct_temp(spin_temp);
                    lowest_energy_spin_struct_set.insert(spin_no_periodic_struct_temp);
                }
                
                
                
            }
            
            set<tuple<int,int,int,int,int> > prototype_set_here;
            bool shall_i_return=false;
            
            //            graph model initiated!
            //            build up models to check connectivity
            //            build up tree
            
            graph_model_and_prototype_set_construction(obscenely_verbose, lowest_energy_spin_struct_set, x_range, y_range, z_range, shall_i_return, model_update_loops, warm_restart, warm_start_spin_structs, warm_start_J, J_input_best_memory, index_to_equivalent_sets, PI_value, revert_index_to_equivalent_sets, prototype_set_here, significant_index);
            set<set<tuple<int,int,int,int,int> > > prototype_set_here_fullset;
            if (prototype_set_here.size()>0) {
                prototype_set_here_fullset.insert(prototype_set_here);
            }
            
            
            if (shall_i_return==true) {
                return;
            }
            
            
            
            
            if (prototype_set_here_fullset.size()==0) {
                //                exhaust to find all the lowest energy states;
                
                {
                    bool new_lowest_states=true;
                    while (new_lowest_states==true) {
                        
                        //                    cout<<"\n now inside debug 44133564136";
                        
                        map<set<tuple<int,int,int,int,int> >, double> J_to_extract_new_GS =J_input_best_memory;
                        for (set<spin_no_periodic_struct>::iterator it1=lowest_energy_spin_struct_set.begin(); it1!=lowest_energy_spin_struct_set.end(); it1++) {
                            map<tuple<int,int,int,int>,int> spin_here=(*it1).spin;
                            
                            set<tuple<int,int,int,int,int> > prototype_set_for_finding_new_lowest_spin;
                            for (map<tuple<int,int,int,int>,int>::iterator it3=spin_here.begin(); it3!=spin_here.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_for_finding_new_lowest_spin.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            
                            J_to_extract_new_GS[prototype_set_for_finding_new_lowest_spin]+=50*max_value_of_J;
                        }
                        
                        
                        map<tuple<int,int,int,int>,int> spin_out;
                        double energy_lowe_bound_dump;
                        map<tuple<int, int, int, int, int, int>, map<set<tuple<int, int, int, int, int> >, double> > clustertypeperiodic_dump;
                        double LB_so_far_dump=-1e20;
                        
                        get_lowerbound_sub(a0, a1, a2, a3, a4, a5, J_to_extract_new_GS, x_range, y_range, z_range, componentnumber, spin_out, energy_lowe_bound_dump, clustertypeperiodic_dump, LB_so_far_dump, id, 0, 0, 1, obscenely_verbose, global_parameters);
                        
                        if(obscenely_verbose)
                            cout<<"\n min energy is "<<min <<" energy_lowe_bound_dump is:"<<energy_lowe_bound_dump;
                        new_lowest_states=false;
                        if (energy_lowe_bound_dump<min+1e-5*max_value_of_J) {
                            if(obscenely_verbose)
                                cout<<"\nfound new low energy states";
                            new_lowest_states=true;
                            spin_no_periodic_struct new_spin_states(spin_out);
                            lowest_energy_spin_struct_set.insert(new_spin_states);
                        }
                        
                        
                    }
                }
                
                cout<<"\n now we are in a very difficult situation and have to add everything in";
                
                
                for (set<spin_no_periodic_struct>::iterator it1=lowest_energy_spin_struct_set.begin(); it1!=lowest_energy_spin_struct_set.end(); it1++) {
                    
                    spin_no_periodic_struct spin_struct_now=*it1;
                    
                    //extract prototype set here
                    {
                        if (x_range>1)
                        {
                            
                            //                        cout<<"\n what is the initial spin\n";
                            //                        printblock(spin_struct_now.spin);
                            spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_x_pos();
                            //                        cout<<"\n what is the truncated spin\n";
                            //                        printblock(truncated_initial_spin_struct_now.spin);
                            //firstly_enforce_we_do_not_loop_back
                            
                            
                            bool already_in_model=false;
                            set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                            construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                            if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                                if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                                    already_in_model=true;
                                }
                            }
                            //end of enforce
                            
                            
                            
                            prototype_set_here.clear();
                            if (already_in_model==false) {
                                spin_no_periodic_struct to_compare_spin_struct;
                                bool inside_some=false;
                                
                                if (inside_some==false) {
                                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                        
                                        tuple<int,int,int,int> initial_tuple=it3->first;
                                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                        int initial_type=it3->second;
                                        prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                                    }
                                    
                                    prototype_set_here_fullset.insert(prototype_set_here);
                                }
                            }
                            
                            
                        }
                        if (x_range>1)
                        {
                            spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_x_neg();
                            
                            //firstly_enforce_we_do_not_loop_back
                            bool already_in_model=false;
                            set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                            construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                            if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                                if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                                    already_in_model=true;
                                }
                            }
                            //end of enforce
                            prototype_set_here.clear();
                            if (already_in_model==false) {
                                spin_no_periodic_struct to_compare_spin_struct;
                                bool inside_some=false;
                                
                                if (inside_some==false) {
                                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                        
                                        tuple<int,int,int,int> initial_tuple=it3->first;
                                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                        int initial_type=it3->second;
                                        prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                                    }
                                    prototype_set_here_fullset.insert(prototype_set_here);
                                }
                            }
                            
                            
                        }
                        if (y_range>1)
                        {
                            spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_y_pos();
                            
                            //firstly_enforce_we_do_not_loop_back
                            bool already_in_model=false;
                            set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                            construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                            if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                                if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                                    already_in_model=true;
                                }
                            }
                            //end of enforce
                            prototype_set_here.clear();
                            if (already_in_model==false) {
                                spin_no_periodic_struct to_compare_spin_struct;
                                bool inside_some=false;
                                if (inside_some==false) {
                                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                        
                                        tuple<int,int,int,int> initial_tuple=it3->first;
                                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                        int initial_type=it3->second;
                                        prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                                    }
                                    prototype_set_here_fullset.insert(prototype_set_here);
                                }
                            }
                            
                            
                        }
                        if (y_range>1)
                        {
                            spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_y_neg();
                            
                            //firstly_enforce_we_do_not_loop_back
                            bool already_in_model=false;
                            set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                            construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                            if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                                if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                                    already_in_model=true;
                                }
                            }
                            //end of enforce
                            prototype_set_here.clear();
                            if (already_in_model==false) {
                                spin_no_periodic_struct to_compare_spin_struct;
                                bool inside_some=false;
                                if (inside_some==false) {
                                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                        
                                        tuple<int,int,int,int> initial_tuple=it3->first;
                                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                        int initial_type=it3->second;
                                        prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                                    }
                                    prototype_set_here_fullset.insert(prototype_set_here);
                                }
                            }
                            
                            
                        }
                        if (z_range>1)
                        {
                            spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_z_pos();
                            
                            //firstly_enforce_we_do_not_loop_back
                            bool already_in_model=false;
                            set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                            construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                            if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                                if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                                    already_in_model=true;
                                }
                            }
                            //end of enforce
                            
                            prototype_set_here.clear();
                            if (already_in_model==false) {
                                spin_no_periodic_struct to_compare_spin_struct;
                                bool inside_some=false;
                                if (inside_some==false) {
                                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                        
                                        tuple<int,int,int,int> initial_tuple=it3->first;
                                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                        int initial_type=it3->second;
                                        prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                                    }
                                    prototype_set_here_fullset.insert(prototype_set_here);
                                }
                            }
                            
                            
                        }
                        if (z_range>1)
                        {
                            spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_z_neg();
                            
                            //firstly_enforce_we_do_not_loop_back
                            bool already_in_model=false;
                            set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                            construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                            if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                                if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                                    already_in_model=true;
                                }
                            }
                            //end of enforce
                            
                            prototype_set_here.clear();
                            if (already_in_model==false) {
                                spin_no_periodic_struct to_compare_spin_struct;
                                bool inside_some=false;
                                if (inside_some==false) {
                                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                        
                                        tuple<int,int,int,int> initial_tuple=it3->first;
                                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                        int initial_type=it3->second;
                                        prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                                    }
                                    prototype_set_here_fullset.insert(prototype_set_here);
                                }
                            }
                            
                            
                            
                        }
                    }
                    
                }
                
            }
            
            
            
            
            
            if (prototype_set_here_fullset.size()==0) {
                if (obscenely_verbose) {
                    cout<<"\nGraph compatability cannot be achieved";
                }
                break;

            }
            
            
            assert(prototype_set_here_fullset.size()>0);
            
            warm_start_spin_structs=lowest_energy_spin_struct_set;
            
            
            for (set<set<tuple<int,int,int,int,int> > >::iterator it1=prototype_set_here_fullset.begin(); it1!=prototype_set_here_fullset.end(); it1++) {
                
                //                set<tuple<int,int,int,int,int> > prototype_set=prototype_set_here;
                set<tuple<int,int,int,int,int> > prototype_set=*it1;
                
                //                assert no empty subset
                assert(prototype_set.size()>0);
                
                
                set<set<tuple<int,int,int,int,int> > > set_of_equivalent;
                construct_set_of_equivalent_from_prototype_set(set_of_equivalent, prototype_set, x_range, y_range, z_range);
                
                if (obscenely_verbose) {
                    cout<<"\ndoing very painful debugging here #14672136464115\n";
                    cout<<"what is prototype_set\n";
                    printvector(prototype_set);
                    cout<<"\nwhat is set_of_equivalent";
                    printvectorofvector(set_of_equivalent);
                    cout<<"\nwhat is revert_index_to_equivalent_sets.count(set_of_equivalent):"<<revert_index_to_equivalent_sets.count(set_of_equivalent);
                    if (revert_index_to_equivalent_sets.count(set_of_equivalent)>0.1) {
                        cout<<"\nwhat is significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent]) "<<significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent]);
                    }
                    
                    
                    
                }
                
                
                if (revert_index_to_equivalent_sets.count(set_of_equivalent) == 0) {
                    index_to_equivalent_sets[equivalent_sets_count] = set_of_equivalent;
                    revert_index_to_equivalent_sets[set_of_equivalent] = equivalent_sets_count;
                    value_of_equivalent_sets[equivalent_sets_count] = 0;
                    significant_index.insert(equivalent_sets_count);
                    assert(warm_start_J.count(prototype_set)==0);
                    warm_start_J[prototype_set]=0;
                    equivalent_sets_count++;
                } else {
                    long long number = revert_index_to_equivalent_sets[set_of_equivalent];
                    significant_index.insert(number);
                }
                
                
            }
            
            
            //            if(index_to_equivalent_sets.count(equivalent_sets_count)==0)
            //            {
            //                    cout<<"\nthis is very weird, return";
            //                return;
            //            }
            
            assert(index_to_equivalent_sets.count(equivalent_sets_count)==0);
            
            
            //            construct graph ;
            
            if (model_update_loops+1==model_number_max) {
                cout<<"\ngraph compatability not achieved at model_number_max= "<<model_number_max<<endl;
            }
            
            
            
            
        }//model loop end
        
    }
    else if (use_weighted_dual_average){
        // Define optimization range (based on loopnumber and dimensionality)
        x_range += a0-1;
        y_range += a2-1;
        z_range += a5-1;
        
        // Set up Gurobi optimizer
        GRBEnv env = GRBEnv();
        GRBModel m = GRBModel(env);
        if (!obscenely_verbose){
            m.getEnv().set(GRB_IntParam_OutputFlag, 0);
        }
        m.update();
        m.set(GRB_IntAttr_ModelSense,-1);
        
        long long equivalent_sets_count = 0;
        map<long long, set<set<tuple<int,int,int,int,int> > > > index_to_equivalent_sets;
        map<set<set<tuple<int,int,int,int,int> > >,long long> revert_index_to_equivalent_sets;
        map<long long, double> value_of_equivalent_sets;
        double max_value_of_J = 0;
        
        for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1 = J.begin(); it1 != J.end(); it1++) {
            set<tuple<int,int,int,int,int> > prototype_set=it1->first;
            double value = it1->second;
            
            // Set some limits on x, y, z range
            int x_max = -1e5;
            int x_min = 1e5;
            int y_max = -1e5;
            int y_min = 1e5;
            int z_max = -1e5;
            int z_min = 1e5;
            for (set<tuple<int,int,int,int,int> >::iterator it2=prototype_set.begin(); it2!=prototype_set.end(); it2++) {
                tuple<int,int,int,int,int> tuple_element=*it2;
                int x_position=tuple_element.get<0>();
                int y_position=tuple_element.get<1>();
                int z_position=tuple_element.get<2>();
                if (x_min > x_position) { x_min = x_position; }
                if (x_max < x_position) { x_max = x_position; }
                if (y_min > y_position) { y_min = y_position; }
                if (y_max < y_position) { y_max = y_position; }
                if (z_min > z_position) { z_min = z_position; }
                if (z_max < z_position) { z_max = z_position; }
            }
            
            // Find equivalent sets of ECIs
            set<set<tuple<int,int,int,int,int> > > set_of_equivalent;
            for (int translation_z = 1-z_min; translation_z <= z_range-z_max; translation_z++) {
                for (int translation_y = 1-y_min; translation_y <= y_range-y_max; translation_y++) {
                    for (int translation_x = 1-x_min; translation_x <= x_range-x_max; translation_x++) {
                        set<tuple<int,int,int,int,int> > temp_equivalent_set;
                        
                        for (set<tuple<int,int,int,int,int> >::iterator it2=prototype_set.begin(); it2!=prototype_set.end(); it2++) {
                            tuple<int,int,int,int,int> tuple_element=*it2;
                            tuple<int,int,int,int,int> new_tuple=make_tuple(tuple_element.get<0>() +
                                                                            translation_x,tuple_element.get<1>() +
                                                                            translation_y,tuple_element.get<2>() +
                                                                            translation_z,tuple_element.get<3>(),
                                                                            tuple_element.get<4>());
                            temp_equivalent_set.insert(new_tuple);
                        }
                        set_of_equivalent.insert(temp_equivalent_set);
                    }
                }
            }
            
            if (revert_index_to_equivalent_sets.count(set_of_equivalent) == 0) {
                index_to_equivalent_sets[equivalent_sets_count] = set_of_equivalent;
                revert_index_to_equivalent_sets[set_of_equivalent] = equivalent_sets_count;
                value_of_equivalent_sets[equivalent_sets_count] = value;
                equivalent_sets_count++;
            } else {
                long long number = revert_index_to_equivalent_sets[set_of_equivalent];
                value_of_equivalent_sets[number] += value;
            }
        }
        
        for (map<long long, double>::iterator it1 = value_of_equivalent_sets.begin(); it1 != value_of_equivalent_sets.end(); it1++) {
            if (max_value_of_J < abs(it1->second)) { max_value_of_J = abs(it1->second); }
        }
        
        map<tuple<long long,long long>, set<tuple<int,int,int,int,int> > > index_to_subcluster;
        map<set<tuple<int,int,int,int,int> >, tuple<long long,long long> > revert_index_to_subcluster;
        map<tuple<long long,long long>, double> subcluster_value;
        map<tuple<long long,long long>, GRBVar> subcluster_variable;
        map<long long,long long> prenum_to_subnum; //Last one does not count
        
        //WX 2015-03-17 begin
        
        long long count_dimension=0;
        vector<tuple<double,long long> > vector_pair_for_value_of_equivalent_sets;
        set<long long> significant_index;
        map<set<tuple<int,int,int,int,int> >, double> J_other_fixed;
        
        for (map<long long, double>::iterator it1=value_of_equivalent_sets.begin(); it1!=value_of_equivalent_sets.end();it1++ ) {
            vector_pair_for_value_of_equivalent_sets.push_back(make_tuple(fabs(it1->second),it1->first));
        }
        
        sort(vector_pair_for_value_of_equivalent_sets.rbegin(), vector_pair_for_value_of_equivalent_sets.rend());
        
        for (vector<tuple<double,long long> >::iterator it1=vector_pair_for_value_of_equivalent_sets.begin(); it1!=vector_pair_for_value_of_equivalent_sets.end(); it1++) {
            tuple<double,long long> tuple_now=*it1;
            long long index=tuple_now.get<1>();
            long long size=index_to_equivalent_sets[index].size();
            if (count_dimension<=limit_dimension) {
                significant_index.insert(index);
                count_dimension+=size;
            }
            else if (count_dimension>limit_dimension){
                break;
            }
        }
        
        
        
        //WX 2015-03-17 end
        
        
        
        for (map<long long, set<set<tuple<int,int,int,int,int> > > >::iterator it1 = index_to_equivalent_sets.begin(); it1 != index_to_equivalent_sets.end(); it1++) {
            if (significant_index.count(it1->first)) {
                long long sub_number = 0;
                long long pre_number = it1->first;
                set<set<tuple<int,int,int,int,int> > > equivalent_set_here = it1->second;
                for (set<set<tuple<int,int,int,int,int> > >::iterator it2=equivalent_set_here.begin(); it2 != equivalent_set_here.end(); it2++) {
                    set<tuple<int,int,int,int,int> > subset_here = *it2;
                    index_to_subcluster[make_tuple(pre_number,sub_number)] = subset_here;
                    revert_index_to_subcluster[subset_here]=make_tuple(pre_number, sub_number);
                    subcluster_variable[make_tuple(pre_number,sub_number)] = m.addVar(-max_value_of_J, max_value_of_J, 0, GRB_CONTINUOUS);
                    
                    sub_number++;
                }
                prenum_to_subnum[pre_number]=sub_number;
            }
            else {
                long long pre_number = it1->first;
                set<set<tuple<int,int,int,int,int> > > equivalent_set_here = it1->second;
                for (set<set<tuple<int,int,int,int,int> > >::iterator it2=equivalent_set_here.begin(); it2 != equivalent_set_here.end(); it2++) {
                    set<tuple<int,int,int,int,int> > subset_here = *it2;
                    J_other_fixed[subset_here]=value_of_equivalent_sets[pre_number]/index_to_equivalent_sets[pre_number].size();
                }
            }
        }
        
        m.update();
        
        
        
        for (map<long long, long long>::iterator it1=prenum_to_subnum.begin(); it1!=prenum_to_subnum.end(); it1++) {
            long long prenumber=it1->first;
            GRBLinExpr sum = 0;
            for (long long subnumber = 0; subnumber < prenum_to_subnum[prenumber]; subnumber++) {
                sum += subcluster_variable[make_tuple(prenumber,subnumber)];
            }
            m.addConstr(sum == value_of_equivalent_sets[prenumber]);
        }
        
        
        
        map<set<tuple<int,int,int,int,int> >, double> J_input;
        double proximal_factor=1e-8;
        GRBQuadExpr objective = 0;
        for (map<long long, long long>::iterator it1=prenum_to_subnum.begin(); it1!=prenum_to_subnum.end(); it1++) {
            long long prenumber=it1->first;
            for (long long subnumber = 0; subnumber < prenum_to_subnum[prenumber]; subnumber++) {
                double prefactor = proximal_factor;
                objective += (-1)*prefactor*subcluster_variable[make_tuple(prenumber,subnumber)]*
                subcluster_variable[make_tuple(prenumber,subnumber)];
                
            }
            
        }
        
        m.setObjective(objective);
        
        //        m_quadratic.setObjective(quadratic_objective);
        
        double lowerbound_memory=-1e20;
        map<set<tuple<int,int,int,int,int> >, double> J_input_best_memory;
        double upperbound_of_lowerbound=1e20;
        
        map<long long,double> big_S_k,beta_hat,beta,lambda_k;
        map<long long,map<tuple<long long,long long>,double > > x_k,s_k,s_hat_k,g_k,x_hat_k;
        double D=(max_value_of_J)*(max_value_of_J);
        double L=sqrt(count_dimension);
        double gamma=L/sqrt(D);
        
        beta_hat[0]=1;
        beta_hat[1]=1;
        big_S_k[-1]=0;
        
        GRBQuadExpr d_x;
        map<long long, GRBLinExpr> gap_function_objective_k;
        GRBModel gap_function=m;
        
        for (long long infinite_loop = 0; infinite_loop < 1000; infinite_loop++) {
            
            
            
            m.optimize();
            
            lambda_k[infinite_loop]=1;
            big_S_k[infinite_loop]=big_S_k[infinite_loop-1]+lambda_k[infinite_loop];
            
            if (infinite_loop>0.1) {
                beta_hat[infinite_loop+1]=beta_hat[infinite_loop]+1/beta_hat[infinite_loop];
            }
            
            
            for (map<tuple<long long,long long>, GRBVar>::iterator it1=subcluster_variable.begin(); it1!=subcluster_variable.end(); it1++) {
                tuple<long long,long long> index=it1->first;
                GRBVar variable=it1->second;
                x_k[infinite_loop][index]=variable.get(GRB_DoubleAttr_X);
            }
            
            if (infinite_loop==0) {
                for (map<tuple<long long,long long>,double >::iterator it1=x_k[infinite_loop].begin(); it1!=x_k[infinite_loop].end(); it1++) {
                    tuple<long long,long long> tuple_now=it1->first;
                    double value=it1->second;
                    x_hat_k[infinite_loop][tuple_now]=lambda_k[infinite_loop]*value/big_S_k[infinite_loop];
                }
            }
            else{
                for (map<tuple<long long,long long>,double >::iterator it1=x_k[infinite_loop].begin(); it1!=x_k[infinite_loop].end(); it1++) {
                    tuple<long long,long long> tuple_now=it1->first;
                    double value=it1->second;
                    x_hat_k[infinite_loop][tuple_now]=(lambda_k[infinite_loop]*value+big_S_k[infinite_loop-1]*x_hat_k[infinite_loop-1][tuple_now])/big_S_k[infinite_loop];
                }
            }
            
            
            
            
            J_input.clear();
            J_input=J_other_fixed;
            
            for (map<tuple<long long,long long>, GRBVar>::iterator it1 = subcluster_variable.begin(); it1 != subcluster_variable.end(); it1++) {
                
                tuple<long long,long long> tuple_here = it1->first;
                GRBVar variable_here=it1->second;
                try {
                    J_input[index_to_subcluster[tuple_here]] = variable_here.get(GRB_DoubleAttr_X);
                } catch (GRBException a) {
                    lower_bound = lowerbound_memory;
                    J_for_proof = J_input_best_memory;
                    if (obscenely_verbose){
                        cout << "LB - Loop #" << infinite_loop <<": Either a numerical error occurred or found tightest LB: " << lower_bound << endl;
                        //                    cout << "J_for_proof: ";
                        //                    printmapfromsets(J_for_proof);
                        cout << "\n" << endl;
                    }
                    return;
                }
            }
            
            if (pseudo_mode) {
                get_lowerbound_sub(1, 0, 1, 0, 0, 1,
                                   J_input,
                                   x_range, y_range, z_range, componentnumber,
                                   spin, lower_bound,
                                   clustertypeperiodic, lowerbound_memory, id,
                                   true, false, false,
                                   obscenely_verbose, global_parameters );
            }else if (basic_exact_mode) {
                get_lowerbound_sub(1, 0, 1, 0, 0, 1,
                                   J_input,
                                   x_range, y_range, z_range, componentnumber,
                                   spin, lower_bound,
                                   clustertypeperiodic, lowerbound_memory, id,
                                   false, true, false,
                                   obscenely_verbose, global_parameters);
            }
            
            map<set<tuple<int,int,int,int,int> >, double> cluster_type_now = clustertypeperiodic[make_tuple(1,0,1,0,0,1)];
            
            if (lowerbound_memory <= lower_bound) {
                lowerbound_memory = lower_bound;
                J_input_best_memory = J_input;
                J_for_proof = J_input_best_memory;
            }
            
            for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1 = cluster_type_now.begin(); it1 != cluster_type_now.end(); it1++) {
                if (revert_index_to_subcluster.count(it1->first)==1) {
                    g_k[infinite_loop][revert_index_to_subcluster[it1->first]]=it1->second;
                    s_k[infinite_loop+1][revert_index_to_subcluster[it1->first]]=s_k[infinite_loop][revert_index_to_subcluster[it1->first]]+lambda_k[infinite_loop]*g_k[infinite_loop][revert_index_to_subcluster[it1->first]];
                }
            }
            
            
            for (map<tuple<long long,long long>,double >::iterator it1=s_k[infinite_loop].begin(); it1!=s_k[infinite_loop].end(); it1++) {
                tuple<long long,long long> tuple_now=it1->first;
                double value=it1->second;
                s_hat_k[infinite_loop][tuple_now]=s_k[infinite_loop][tuple_now]/big_S_k[infinite_loop];
            }
            
            
            
            
            if (infinite_loop==0) {
                d_x=0;
                for (map<tuple<long long,long long>, GRBVar>::iterator it1=subcluster_variable.begin(); it1!=subcluster_variable.end(); it1++) {
                    tuple<long long,long long> tuple_now=it1->first;
                    GRBVar variable_now=it1->second;
                    d_x-=(variable_now-x_k[0][tuple_now])*(variable_now-x_k[0][tuple_now]);
                }
                
            }
            
            
            
            
            
            //            double gap_function_value=0;
            //
            //            if (infinite_loop==0) {
            //                gap_function.addQConstr(d_x>=-D);
            //
            //            }
            //            //                gap_function_objective+=lambda_k[i]*<g_k[i],x-x_k>;
            //            if (infinite_loop==0) {
            //                gap_function_objective_k[infinite_loop]=0;
            //
            //            }
            //            else {
            //                gap_function_objective_k[infinite_loop]=gap_function_objective_k[infinite_loop-1];
            //            }
            //
            //            for (map<tuple<long long,long long>,double >::iterator it1=g_k[infinite_loop].begin(); it1!=g_k[infinite_loop].end(); it1++) {
            //                tuple<long long,long long> tuple_now=it1->first;
            //                double value=it1->second;
            //                gap_function_objective_k[infinite_loop]+=lambda_k[infinite_loop]*g_k[infinite_loop][tuple_now]*(subcluster_variable[tuple_now]-x_k[infinite_loop][tuple_now]);
            //            }
            //
            //            gap_function.setObjective(gap_function_objective_k[infinite_loop]);
            //            gap_function.optimize();
            //            gap_function_value=gap_function.get(GRB_DoubleAttr_ObjVal)/big_S_k[infinite_loop];
            
            
            
            
            
            
            
            
            
            
            GRBQuadExpr objective=0;
            
            for (map<tuple<long long,long long>,double > ::iterator it1=s_k[infinite_loop+1].begin(); it1!=s_k[infinite_loop+1].end(); it1++) {
                
                tuple<long long,long long> tuple_now=it1->first;
                double value=it1->second;
                objective+=subcluster_variable[tuple_now]*value;
                
            }
            
            beta[infinite_loop+1]=gamma*beta_hat[infinite_loop+1];
            objective+=beta[infinite_loop+1]*d_x;
            m.setObjective(objective);
            
            
            
            
            if (obscenely_verbose){
                cout << "\nLB - Loop #" << infinite_loop << endl;
                //                cout << "\tBound of LB: " << upperbound_of_lowerbound << endl;
                cout << "\tCurrent LB: " << lower_bound << endl;
                cout << "\tBest LB: " << lowerbound_memory<< endl;
                //                cout << "\tGap_function value: "<<gap_function_value<<endl;
                
                //                cout << "\tProximity factor:" << proximal_factor << endl;
                
                cout << "\n" << endl;
            }
            
            if (lower_bound + 0.1 > upperbound_of_lowerbound) {
                J_for_proof = J_input_best_memory;
                lower_bound = lowerbound_memory;
                if (obscenely_verbose){
                    cout << "LB - Loop #" << infinite_loop << ": Found tightest lower bound: " << lower_bound << endl;
                    cout << "J_for_proof:";
                    //                    printmapfromsets(J_for_proof);
                    cout << "\n" << endl;
                }
                break;
            }
            
            
        }
        lower_bound=lowerbound_memory;
    }
  
    
    
}




void get_lowerbound_sub(int a0, int a1, int a2, int a3, int a4, int a5,
                        map<set<tuple<int,int,int,int,int> >, double> &J,
                        int x_range, int y_range, int z_range, map<int,int> componentnumber,
                        map<tuple<int,int,int,int>,int> &spin,
                        double &lower_bound,
                        map<tuple<int,int,int,int,int,int>, map<set<tuple<int,int,int,int,int> >, double> > &clustertypeperiodic,
                        double LB_so_far,
                        std::string id,
                        bool stop_at_first_iteration,
                        bool stop_when_upperbound_smaller_than_predicted_lower_bound,
                        bool stop_till_exact_result,
                        bool obscenely_verbose,solver_variable &global_parameters)
{
    
    srand(time(NULL));
    
    if (obscenely_verbose) {
        cout<<"\n I am here #48761354, I want to examine what J is";
        printmapfromsets(J);
    }
    
    
    
    bool is_zero_involved=false;
    for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1=J.begin(); it1!=J.end(); it1++) {
        set<tuple<int,int,int,int,int> > temp_set=it1->first;
        for (set<tuple<int,int,int,int,int> >::iterator it2=temp_set.begin(); it2!=temp_set.end(); it2++) {
            tuple<int,int,int,int,int> temp_tuple=*it2;
            int type=temp_tuple.get<4>();
            if (type==0&&it1->second!=0)
            {
                is_zero_involved=true;
                break;
            }
            
        }
        if (is_zero_involved==true) {
            break;
        }
    }
    
    if (obscenely_verbose)
        cout<<"\n is_zero_involved: "<<is_zero_involved<<endl;
    
    
    if(is_zero_involved==false)
    {
        
        
        
        
        if (obscenely_verbose) {
            cout << "\nCalculating LB with a0 = " << a0;
            cout << ", a1 = " << a1;
            cout << ", a2 = " << a2;
            cout << ", a3 = " << a3;
            cout << ", a4 = " << a4;
            cout << ", a5 = " << a5 << endl;
        }
        
        set< set<tuple<int,int,int,int,int> > > nonzerolist;
        map<set<tuple<int,int,int,int,int> >, long long> J_integral;
        long long softer_count_variable = 1;
        map<tuple<int,int,int,int,int>, long long> softer_variable_encode;
        map<long long,tuple<int,int,int,int,int> > softer_variable_decode;
        map<set<long long>, long long> softer_maxsat_model;
        set<set<long long> > softer_maxsat_model_hard;
        long long softer_max_element = 1e18;
        map<long long, int> soft_result;
        long long overall_constant = 0;
        
        for(map<set<tuple<int,int,int,int,int> >, double>::iterator it = J.begin(); it != J.end(); it++){
            nonzerolist.insert(it->first);
            J_integral[it->first] = round(it->second);
        }
        
        for (int i = 1; i < a0 + x_range; i++) {
            for(int j = 1; j < a2 + y_range; j++) {
                for (int k = 1; k < a5 + z_range; k++) {
                    for (map<int, int>::iterator cit = componentnumber.begin(); cit != componentnumber.end(); cit++) {
                        int position = cit->first;
                        int element_count = cit->second;
                        for (int z = 0; z < element_count; z++) {
                            if (z > 0.1) {
                                softer_variable_encode[make_tuple(i,j,k,position,z)] = softer_count_variable;
                                softer_variable_decode[softer_count_variable] = make_tuple(i,j,k,position,z);
                                softer_count_variable++;
                            }
                        }
                    }
                }
            }
        }
        
        
        
        for (int i = 1; i < a0 + x_range; i++) {
            for (int j = 1; j < a2 + y_range; j++) {
                for (int k = 1; k < a5 + z_range; k++) {
                    for (map<int, int>::iterator cit = componentnumber.begin(); cit != componentnumber.end(); cit++) {
                        int position = cit->first;
                        int element_count = cit->second;
                        for (int temp5 = 0; temp5 < element_count - 1; temp5++) {
                            for (int temp6 = temp5 + 1; temp6 < element_count; temp6++) {
                                set<long long> set_temp;
                                if (temp5 > 0.1) {
                                    set<long long> softer_set_temp;
                                    softer_set_temp.insert(-softer_variable_encode[make_tuple(i, j, k, position, temp5)]);
                                    softer_set_temp.insert(-softer_variable_encode[make_tuple(i, j, k, position, temp6)]);
                                    softer_maxsat_model_hard.insert(softer_set_temp);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        
        for (int i = 1; i < a0 + 1; i++) {
            for (int j = 1; j < a2 + 1; j++) {
                for (int k = 1; k < a5 + 1; k++) {
                    for(set< set<tuple<int,int,int,int,int> > >::iterator itv = nonzerolist.begin();
                        itv != nonzerolist.end(); ++itv){
                        tuple<int,int,int,int,int> first, second;
                        set<tuple<int,int,int,int,int> > tuple_set = *itv;
                        if (J_integral[tuple_set] > 0.1) {
                            set<long long> set_of_terms;
                            set<long long> softer_set_of_terms;
                            for(set< tuple<int,int,int,int,int> >::iterator itv1 = tuple_set.begin();
                                itv1 != tuple_set.end(); itv1++){
                                tuple<int,int,int,int,int> tuple_here = *itv1;
                                int x0 = tuple_here.get<0>();
                                int y0 = tuple_here.get<1>();
                                int z0 = tuple_here.get<2>();
                                int position = tuple_here.get<3>();
                                int mu = tuple_here.get<4>();
                                softer_set_of_terms.insert(-softer_variable_encode[make_tuple(x0 + i - 1,
                                                                                              y0 + j - 1,
                                                                                              z0 + k - 1,
                                                                                              position, mu)]);
                            }
                            softer_maxsat_model[softer_set_of_terms] += J_integral[tuple_set];
                        }
                        
                        if (J_integral[tuple_set] < -0.1) {
                            for(set<tuple<int,int,int,int,int> >::iterator itv1 = tuple_set.begin();
                                itv1 != tuple_set.end(); itv1++){
                                set<long long> set_of_terms;
                                set<long long> softer_set_of_terms;
                                tuple<int,int,int,int,int> tuple_here = *itv1;
                                int x0 = tuple_here.get<0>();
                                int y0 = tuple_here.get<1>();
                                int z0 = tuple_here.get<2>();
                                int position = tuple_here.get<3>();
                                int mu = tuple_here.get<4>();
                                softer_set_of_terms.insert(softer_variable_encode[make_tuple(x0 + i - 1,
                                                                                             y0 + j - 1,
                                                                                             z0 + k - 1,
                                                                                             position, mu)]);
                                itv1++;
                                for (set<tuple<int,int,int,int,int> >::iterator itv2 = itv1;
                                     itv2 != tuple_set.end(); itv2++) {
                                    tuple<int,int,int,int,int> tuple_here = *itv2;
                                    int x0 = tuple_here.get<0>();
                                    int y0 = tuple_here.get<1>();
                                    int z0 = tuple_here.get<2>();
                                    int position = tuple_here.get<3>();
                                    int mu = tuple_here.get<4>();
                                    softer_set_of_terms.insert(-softer_variable_encode[make_tuple(x0 + i - 1,
                                                                                                  y0 + j - 1,
                                                                                                  z0 + k - 1,
                                                                                                  position, mu)]);
                                }
                                itv1--;
                                softer_maxsat_model[softer_set_of_terms]-=J_integral[tuple_set];
                            }
                            overall_constant -= J_integral[tuple_set];
                        }
                    }
                }
            }
        }
        
        
        // Construct input file for MAXSAT optimizer
        ofstream softer_input;
        string softer_filename=id + "_" + to_string_is(a0) + "." + to_string_is(a2) + "." + to_string_is(a5);
        softer_filename += "_softer_lowerbound.wcnf";
        softer_input.open(softer_filename.c_str());
        softer_input << "p wcnf " << softer_variable_encode.size() << " "
        << softer_maxsat_model.size() + softer_maxsat_model_hard.size() << " " << softer_max_element;
        
        
        for (map<set<long long>, long long>::iterator it1 = softer_maxsat_model.begin();
             it1 != softer_maxsat_model.end(); it1++) {
            set<long long> the_set = it1->first;
            softer_input << "\n" << (it1->second);
            for (set<long long>::iterator it2 = the_set.begin(); it2 != the_set.end(); it2++) {
                softer_input << " " << (*it2);
            }
            softer_input << " 0";
        }
        
        
        for (set<set<long long> >::iterator it1 = softer_maxsat_model_hard.begin();
             it1 != softer_maxsat_model_hard.end(); it1++) {
            set<long long> the_set = *it1;
            softer_input << "\n" << softer_max_element;
            for (set<long long>::iterator it2 = the_set.begin(); it2 != the_set.end(); it2++) {
                softer_input << " " << *it2;
            }
            softer_input << " 0";
        }
        
        
        softer_input.close();
        
        long long criteria = overall_constant + round(LB_so_far * a0 * a2 * a5) - 100;
        if (stop_when_upperbound_smaller_than_predicted_lower_bound) {
            criteria = (overall_constant + round(LB_so_far * a0 * a2 * a5) - 100)*0.95;
        }
        if (criteria < 0) { criteria = 0; }
        if (stop_at_first_iteration) { criteria = 1e18; }
        if (stop_till_exact_result) { criteria = 0; }
        
        
        //use stop_till_exact_result for robustness
        criteria = 0;
        
        string criteria_string = lexical_cast<string>(criteria);
        string result_now = exec("./CCLS_to_akmaxsat_LB " + softer_filename + " " + criteria_string);
        std::remove(softer_filename.c_str());
        
        vector<string> result_pieces;
        split(result_now, '\n', result_pieces);
        string s_line;
        vector<string> v_line;
        
        //long long o_value;
        for (vector<string>::iterator it=result_pieces.begin(); it != result_pieces.end(); it++) {
            string temp_string = *it;
            vector<string> temp_segment;
            split(temp_string, ' ', temp_segment);
            string first_segment = *temp_segment.begin();
            if (first_segment == "s") {
                s_line = temp_string;
            }
            //if (first_segment=="o") {
            //vector<string>::iterator it1=temp_segment.begin();
            //o_value=boost::lexical_cast<long long>(*(++it1));
            //}
            if (first_segment == "v") {
                v_line.push_back(temp_string);
                for (vector<string>::iterator it1 = temp_segment.begin(); it1 != temp_segment.end(); it1++) {
                    if (it1!=temp_segment.begin()) {
                        long long number_now = boost::lexical_cast<long long>(*it1);
                        if (number_now < -0.1) {
                            soft_result[-number_now] = 0;
                        }
                        else if (number_now > 0.1){
                            soft_result[number_now] = 1;
                        }
                    }
                }
            }
        }
        
        map<tuple<int, int,int,int,int>, int> s_result;
        for (map<long long, int>::iterator it = soft_result.begin(); it != soft_result.end(); it++) {
            s_result[softer_variable_decode[it->first]] = it->second;
        }
        
        for (int i = 1; i < a0 + x_range; i++){
            for (int j = 1; j < a2 + y_range; j++) {
                for (int k = 1; k < a5 + z_range; k++) {
                    for (map<int, int>::iterator cit = componentnumber.begin();
                         cit != componentnumber.end(); cit++) {
                        int position = cit->first;
                        int element_count = cit->second;
                        for (int mu = 0; mu < element_count; mu++) {
                            if (s_result[make_tuple(i, j, k, position, mu)] == 1) {
                                spin[make_tuple(i, j, k, position)] = mu;
                                break;
                            }
                            spin[make_tuple(i, j, k, position)] = 0;
                        }
                    }
                    
                }
            }
        }
        
        
        
        if (obscenely_verbose){
            cout << "\nBlock result is: ";
            printblock(spin);
            //            printmap(spin);
            cout << endl;
        }
        
        for (set< set<tuple<int,int,int,int,int> > >::iterator it = nonzerolist.begin();
             it != nonzerolist.end(); it++) {
            set<tuple<int,int,int,int,int> > thisset = *it;
            clustertypeperiodic[make_tuple(a0,a1,a2,a3,a4,a5)][*it] = 0;
            for (int i = 1; i < a0 + 1; i++) {
                for (int j = 1; j < a2 + 1; j++) {
                    for (int k = 1; k < a5 + 1; k++) {
                        int all_matched_indicator = 0;
                        for (set<tuple<int,int,int,int,int> >::iterator it1 = thisset.begin();
                             it1 != thisset.end(); it1++) {
                            tuple<int,int,int,int,int> temp_tuple = *it1;
                            int x = temp_tuple.get<0>();
                            int y = temp_tuple.get<1>();
                            int z = temp_tuple.get<2>();
                            int position = temp_tuple.get<3>();
                            int var = temp_tuple.get<4>();
                            if (spin[make_tuple(i + x - 1, j + y - 1, k + z - 1, position)] != var) {
                                break;
                            }
                            it1++;
                            if (it1 == thisset.end()) {
                                all_matched_indicator = 1;
                                break;
                            }
                            it1--;
                        }
                        
                        if (all_matched_indicator == 1) {
                            clustertypeperiodic[make_tuple(a0, a1, a2, a3, a4, a5)][*it] += 1.0/(a0 * a2 * a5);
                        }
                    }
                }
            }
        }
        
        map<set<tuple<int,int,int,int,int> >, double> cluster_type_here =
        clustertypeperiodic[make_tuple(a0, a1, a2, a3, a4, a5)];
        
        lower_bound = 0;
        
        for (map<set<tuple<int,int,int,int,int> >, double>::iterator it = cluster_type_here.begin();
             it != cluster_type_here.end(); it++) {
            lower_bound += J[it->first]*it->second;
        }
        
        if (obscenely_verbose){
            cout << "\n(" << id <<") LB average energy: " << lower_bound << endl;;
        }
        
    }
    else if (is_zero_involved==true)
    {
        
        
        
        
        if (obscenely_verbose) {
            cout << "\nCalculating LB with a0 = " << a0;
            cout << ", a1 = " << a1;
            cout << ", a2 = " << a2;
            cout << ", a3 = " << a3;
            cout << ", a4 = " << a4;
            cout << ", a5 = " << a5 << endl;
        }
        
        set< set<tuple<int,int,int,int,int> > > nonzerolist;
        map<set<tuple<int,int,int,int,int> >, long long> J_integral;
        long long softer_count_variable = 1;
        map<tuple<int,int,int,int,int>, long long> softer_variable_encode;
        map<long long,tuple<int,int,int,int,int> > softer_variable_decode;
        map<set<long long>, long long> softer_maxsat_model;
        set<set<long long> > softer_maxsat_model_hard;
        long long softer_max_element = 1e18;
        map<long long, int> soft_result;
        long long overall_constant = 0;
        
        for(map<set<tuple<int,int,int,int,int> >, double>::iterator it = J.begin(); it != J.end(); it++){
            nonzerolist.insert(it->first);
            J_integral[it->first] = round(it->second);
        }
        
        for (int i = 1; i < a0 + x_range; i++) {
            for(int j = 1; j < a2 + y_range; j++) {
                for (int k = 1; k < a5 + z_range; k++) {
                    for (map<int, int>::iterator cit = componentnumber.begin(); cit != componentnumber.end(); cit++) {
                        int position = cit->first;
                        int element_count = cit->second;
                        for (int z = 0; z < element_count; z++) {
                            //                            if (z > 0.1) {
                            softer_variable_encode[make_tuple(i,j,k,position,z)] = softer_count_variable;
                            softer_variable_decode[softer_count_variable] = make_tuple(i,j,k,position,z);
                            softer_count_variable++;
                            //                            }
                        }
                    }
                }
            }
        }
        
        
        
        for (int i = 1; i < a0 + x_range; i++) {
            for (int j = 1; j < a2 + y_range; j++) {
                for (int k = 1; k < a5 + z_range; k++) {
                    for (map<int, int>::iterator cit = componentnumber.begin(); cit != componentnumber.end(); cit++) {
                        int position = cit->first;
                        int element_count = cit->second;
                        for (int temp5 = 0; temp5 < element_count - 1; temp5++) {
                            for (int temp6 = temp5 + 1; temp6 < element_count; temp6++) {
                                set<long long> set_temp;
                                //                                if (temp5 > 0.1) {
                                set<long long> softer_set_temp;
                                softer_set_temp.insert(-softer_variable_encode[make_tuple(i, j, k, position, temp5)]);
                                softer_set_temp.insert(-softer_variable_encode[make_tuple(i, j, k, position, temp6)]);
                                softer_maxsat_model_hard.insert(softer_set_temp);
                                //                                }
                            }
                        }
                        
                        set<long long> set_temp1;
                        for (int temp5 = 0; temp5 < element_count; temp5++) {
                            set_temp1.insert(softer_variable_encode[make_tuple(i, j, k, position, temp5)]);
                        }
                        softer_maxsat_model_hard.insert(set_temp1);
                    }
                }
            }
        }
        
        for (int i = 1; i < a0 + 1; i++) {
            for (int j = 1; j < a2 + 1; j++) {
                for (int k = 1; k < a5 + 1; k++) {
                    for(set< set<tuple<int,int,int,int,int> > >::iterator itv = nonzerolist.begin();
                        itv != nonzerolist.end(); ++itv){
                        tuple<int,int,int,int,int> first, second;
                        set<tuple<int,int,int,int,int> > tuple_set = *itv;
                        if (J_integral[tuple_set] > 0.1) {
                            set<long long> set_of_terms;
                            set<long long> softer_set_of_terms;
                            for(set< tuple<int,int,int,int,int> >::iterator itv1 = tuple_set.begin();
                                itv1 != tuple_set.end(); itv1++){
                                tuple<int,int,int,int,int> tuple_here = *itv1;
                                int x0 = tuple_here.get<0>();
                                int y0 = tuple_here.get<1>();
                                int z0 = tuple_here.get<2>();
                                int position = tuple_here.get<3>();
                                int mu = tuple_here.get<4>();
                                softer_set_of_terms.insert(-softer_variable_encode[make_tuple(x0 + i - 1,
                                                                                              y0 + j - 1,
                                                                                              z0 + k - 1,
                                                                                              position, mu)]);
                            }
                            softer_maxsat_model[softer_set_of_terms] += J_integral[tuple_set];
                        }
                        
                        if (J_integral[tuple_set] < -0.1) {
                            for(set<tuple<int,int,int,int,int> >::iterator itv1 = tuple_set.begin();
                                itv1 != tuple_set.end(); itv1++){
                                set<long long> set_of_terms;
                                set<long long> softer_set_of_terms;
                                tuple<int,int,int,int,int> tuple_here = *itv1;
                                int x0 = tuple_here.get<0>();
                                int y0 = tuple_here.get<1>();
                                int z0 = tuple_here.get<2>();
                                int position = tuple_here.get<3>();
                                int mu = tuple_here.get<4>();
                                softer_set_of_terms.insert(softer_variable_encode[make_tuple(x0 + i - 1,
                                                                                             y0 + j - 1,
                                                                                             z0 + k - 1,
                                                                                             position, mu)]);
                                itv1++;
                                for (set<tuple<int,int,int,int,int> >::iterator itv2 = itv1;
                                     itv2 != tuple_set.end(); itv2++) {
                                    tuple<int,int,int,int,int> tuple_here = *itv2;
                                    int x0 = tuple_here.get<0>();
                                    int y0 = tuple_here.get<1>();
                                    int z0 = tuple_here.get<2>();
                                    int position = tuple_here.get<3>();
                                    int mu = tuple_here.get<4>();
                                    softer_set_of_terms.insert(-softer_variable_encode[make_tuple(x0 + i - 1,
                                                                                                  y0 + j - 1,
                                                                                                  z0 + k - 1,
                                                                                                  position, mu)]);
                                }
                                itv1--;
                                softer_maxsat_model[softer_set_of_terms]-=J_integral[tuple_set];
                            }
                            overall_constant -= J_integral[tuple_set];
                        }
                    }
                }
            }
        }
        
        
        // Construct input file for MAXSAT optimizer
        ofstream softer_input;
        string softer_filename=id + "_" + to_string_is(a0) + "." + to_string_is(a2) + "." + to_string_is(a5);
        
        
        
        
        softer_filename += "_softer_lowerbound_.wcnf";
        softer_input.open(softer_filename.c_str());
        softer_input << "p wcnf " << softer_variable_encode.size() << " "
        << softer_maxsat_model.size() + softer_maxsat_model_hard.size() << " " << softer_max_element;
        
        
        for (map<set<long long>, long long>::iterator it1 = softer_maxsat_model.begin();
             it1 != softer_maxsat_model.end(); it1++) {
            set<long long> the_set = it1->first;
            softer_input << "\n" << (it1->second);
            for (set<long long>::iterator it2 = the_set.begin(); it2 != the_set.end(); it2++) {
                softer_input << " " << (*it2);
            }
            softer_input << " 0";
        }
        
        
        for (set<set<long long> >::iterator it1 = softer_maxsat_model_hard.begin();
             it1 != softer_maxsat_model_hard.end(); it1++) {
            set<long long> the_set = *it1;
            softer_input << "\n" << softer_max_element;
            for (set<long long>::iterator it2 = the_set.begin(); it2 != the_set.end(); it2++) {
                softer_input << " " << *it2;
            }
            softer_input << " 0";
        }
        
        
        softer_input.close();
        
        long long criteria = overall_constant + round(LB_so_far * a0 * a2 * a5) - 100;
        if (stop_when_upperbound_smaller_than_predicted_lower_bound) {
            criteria = (overall_constant + round(LB_so_far * a0 * a2 * a5) - 100)*0.95;

        }
        if (criteria < 0) { criteria = 0; }
        if (stop_at_first_iteration) { criteria = 1e18; }
        if (stop_till_exact_result) { criteria = 0; }
        
        //use stop_till_exact_result for robustness
        criteria = 0;
        
        string criteria_string = lexical_cast<string>(criteria);
        
        string result_now = exec("./CCLS_to_akmaxsat_LB " + softer_filename + " " + criteria_string);
        std::remove(softer_filename.c_str());
        //        cout<<"\nresult_now is:\n"<<result_now;
        
        vector<string> result_pieces;
        split(result_now, '\n', result_pieces);
        string s_line;
        vector<string> v_line;
        
        //long long o_value;
        for (vector<string>::iterator it=result_pieces.begin(); it != result_pieces.end(); it++) {
            string temp_string = *it;
            vector<string> temp_segment;
            split(temp_string, ' ', temp_segment);
            string first_segment = *temp_segment.begin();
            if (first_segment == "s") {
                s_line = temp_string;
            }
            //if (first_segment=="o") {
            //vector<string>::iterator it1=temp_segment.begin();
            //o_value=boost::lexical_cast<long long>(*(++it1));
            //}
            if (first_segment == "v") {
                v_line.push_back(temp_string);
                for (vector<string>::iterator it1 = temp_segment.begin(); it1 != temp_segment.end(); it1++) {
                    if (it1!=temp_segment.begin()) {
                        long long number_now = boost::lexical_cast<long long>(*it1);
                        if (number_now < -0.1) {
                            soft_result[-number_now] = 0;
                        }
                        else if (number_now > 0.1){
                            soft_result[number_now] = 1;
                        }
                    }
                }
            }
        }
        
        map<tuple<int, int,int,int,int>, int> s_result;
        for (map<long long, int>::iterator it = soft_result.begin(); it != soft_result.end(); it++) {
            s_result[softer_variable_decode[it->first]] = it->second;
        }
        
        for (int i = 1; i < a0 + x_range; i++){
            for (int j = 1; j < a2 + y_range; j++) {
                for (int k = 1; k < a5 + z_range; k++) {
                    for (map<int, int>::iterator cit = componentnumber.begin();
                         cit != componentnumber.end(); cit++) {
                        int position = cit->first;
                        int element_count = cit->second;
                        for (int mu = 0; mu < element_count; mu++) {
                            if (s_result[make_tuple(i, j, k, position, mu)] == 1) {
                                spin[make_tuple(i, j, k, position)] = mu;
                                break;
                            }
                        }
                        assert(spin.count(make_tuple(i, j, k, position))==1);
                    }
                    
                }
            }
        }
        
        
        
        if (obscenely_verbose){
            cout << "\nBlock result is: ";
            printblock(spin);
            cout << endl;
        }
        
        for (set< set<tuple<int,int,int,int,int> > >::iterator it = nonzerolist.begin();
             it != nonzerolist.end(); it++) {
            set<tuple<int,int,int,int,int> > thisset = *it;
            clustertypeperiodic[make_tuple(a0,a1,a2,a3,a4,a5)][*it] = 0;
            for (int i = 1; i < a0 + 1; i++) {
                for (int j = 1; j < a2 + 1; j++) {
                    for (int k = 1; k < a5 + 1; k++) {
                        int all_matched_indicator = 0;
                        for (set<tuple<int,int,int,int,int> >::iterator it1 = thisset.begin();
                             it1 != thisset.end(); it1++) {
                            tuple<int,int,int,int,int> temp_tuple = *it1;
                            int x = temp_tuple.get<0>();
                            int y = temp_tuple.get<1>();
                            int z = temp_tuple.get<2>();
                            int position = temp_tuple.get<3>();
                            int var = temp_tuple.get<4>();
                            if (spin[make_tuple(i + x - 1, j + y - 1, k + z - 1, position)] != var) {
                                break;
                            }
                            it1++;
                            if (it1 == thisset.end()) {
                                all_matched_indicator = 1;
                                break;
                            }
                            it1--;
                        }
                        
                        if (all_matched_indicator == 1) {
                            clustertypeperiodic[make_tuple(a0, a1, a2, a3, a4, a5)][*it] += 1.0/(a0 * a2 * a5);
                        }
                    }
                }
            }
        }
        
        map<set<tuple<int,int,int,int,int> >, double> cluster_type_here =
        clustertypeperiodic[make_tuple(a0, a1, a2, a3, a4, a5)];
        
        lower_bound = 0;
        
        for (map<set<tuple<int,int,int,int,int> >, double>::iterator it = cluster_type_here.begin();
             it != cluster_type_here.end(); it++) {
            lower_bound += J[it->first]*it->second;
        }
        
        if (obscenely_verbose){
            cout << "\n(" << id <<") LB average energy: " << lower_bound << endl;;
        }
        
    }
    
}







void graph_model_and_prototype_set_construction(bool obscenely_verbose,set<spin_no_periodic_struct> lowest_energy_spin_struct_set, int x_range, int y_range, int z_range, bool &shall_i_return,long long model_update_loops, bool &warm_restart,set<spin_no_periodic_struct> &warm_start_spin_structs,map<set<tuple<int,int,int,int,int> >, double> &warm_start_J, map<set<tuple<int,int,int,int,int> >, double> J_input_best_memory,map<long long, set<set<tuple<int,int,int,int,int> > > > index_to_equivalent_sets,map<spin_no_periodic_struct, double> PI_value,map<set<set<tuple<int,int,int,int,int> > >,long long> revert_index_to_equivalent_sets,  set<tuple<int,int,int,int,int> > &prototype_set_here,set<long long> significant_index)
{
    
    {
        GRBEnv env_graph = GRBEnv();
        if (!obscenely_verbose) {
            env_graph.set(GRB_IntParam_OutputFlag, 0);
        }
        GRBModel model_graph  = GRBModel(env_graph);
        model_graph.set(GRB_IntAttr_ModelSense, -1);
        
        
        
        map<spin_no_periodic_struct, GRBVar> graph_variable;
        
        //            define variable
        
        for (set<spin_no_periodic_struct>::iterator it1=lowest_energy_spin_struct_set.begin(); it1!=lowest_energy_spin_struct_set.end(); it1++) {
            graph_variable[*it1]=model_graph.addVar(0, 1, 1, GRB_CONTINUOUS);
        }
        
        //            add node in the x graph
        set<spin_no_periodic_struct> spin_struct_x_nodes_set,spin_struct_y_nodes_set,spin_struct_z_nodes_set;
        
        if (x_range>1) {
            for (set<spin_no_periodic_struct>::iterator it1=lowest_energy_spin_struct_set.begin() ; it1!=lowest_energy_spin_struct_set.end(); it1++) {
                spin_no_periodic_struct spin_struct_here=*it1,spin_struct_temp;
                spin_struct_x_nodes_set.insert(spin_struct_here.remove_x_pos());
                spin_struct_x_nodes_set.insert(spin_struct_here.remove_x_neg());
                
            }
        }
        
        if (y_range>1) {
            for (set<spin_no_periodic_struct>::iterator it1=lowest_energy_spin_struct_set.begin() ; it1!=lowest_energy_spin_struct_set.end(); it1++) {
                
                spin_no_periodic_struct spin_struct_here=*it1,spin_struct_temp;
                spin_struct_y_nodes_set.insert(spin_struct_here.remove_y_pos());
                spin_struct_y_nodes_set.insert(spin_struct_here.remove_y_neg());
                
            }
        }
        
        if (z_range>1) {
            for (set<spin_no_periodic_struct>::iterator it1=lowest_energy_spin_struct_set.begin() ; it1!=lowest_energy_spin_struct_set.end(); it1++) {
                
                spin_no_periodic_struct spin_struct_here=*it1,spin_struct_temp;
                spin_struct_z_nodes_set.insert(spin_struct_here.remove_z_pos());
                spin_struct_z_nodes_set.insert(spin_struct_here.remove_z_neg());
                
            }
        }
        
        
        model_graph.update();
        
        for (set<spin_no_periodic_struct>::iterator it1=spin_struct_x_nodes_set.begin(); it1!=spin_struct_x_nodes_set.end(); it1++) {
            spin_no_periodic_struct node_here=*it1;
            GRBLinExpr lhs=0,rhs=0;
            for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                spin_no_periodic_struct spin_struct_arcs=*it2;
                if (spin_struct_arcs.remove_x_neg()==node_here) {
                    lhs+=graph_variable[spin_struct_arcs];
                }
                if (spin_struct_arcs.remove_x_pos()==node_here) {
                    rhs+=graph_variable[spin_struct_arcs];
                }
            }
            model_graph.addConstr(lhs==rhs);
        }
        
        for (set<spin_no_periodic_struct>::iterator it1=spin_struct_y_nodes_set.begin(); it1!=spin_struct_y_nodes_set.end(); it1++) {
            spin_no_periodic_struct node_here=*it1;
            GRBLinExpr lhs=0,rhs=0;
            for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                spin_no_periodic_struct spin_struct_arcs=*it2;
                if (spin_struct_arcs.remove_y_neg()==node_here) {
                    lhs+=graph_variable[spin_struct_arcs];
                }
                if (spin_struct_arcs.remove_y_pos()==node_here) {
                    rhs+=graph_variable[spin_struct_arcs];
                }
            }
            model_graph.addConstr(lhs==rhs);
        }
        
        for (set<spin_no_periodic_struct>::iterator it1=spin_struct_z_nodes_set.begin(); it1!=spin_struct_z_nodes_set.end(); it1++) {
            spin_no_periodic_struct node_here=*it1;
            GRBLinExpr lhs=0,rhs=0;
            for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                spin_no_periodic_struct spin_struct_arcs=*it2;
                if (spin_struct_arcs.remove_z_neg()==node_here) {
                    lhs+=graph_variable[spin_struct_arcs];
                }
                if (spin_struct_arcs.remove_z_pos()==node_here) {
                    rhs+=graph_variable[spin_struct_arcs];
                }
            }
            model_graph.addConstr(lhs==rhs);
        }
        
        
        //            cout<<"\n graph model optimization here"<<endl;
        
        model_graph.optimize();
        
        if (model_graph.get(GRB_DoubleAttr_ObjVal)>0.01) {
            cout<<"\n graph model compatibile at model number: "<<model_update_loops<<" return"<<endl;
            shall_i_return=true;
        }
        
        
    }
    
    //end of graph model checking
    
    //start building for new variables
    
    
    //let's just record and resent as necessary
    
    
    //            remember to update value_of_equivalent_sets,significant index, index_to_equivalent_sets;
    
    warm_restart=true;
    warm_start_spin_structs=lowest_energy_spin_struct_set;
    warm_start_J=J_input_best_memory;
    
    
    
    //            cout<<"\n debug519423 lowest_energy_spin_struct_set.size() is "<<lowest_energy_spin_struct_set.size();
    
    if  (obscenely_verbose){
        cout<<"\ndebug here let me see what is the lowest energy spin struct sets";
        for (set<long long>::iterator it1=significant_index.begin(); it1!=significant_index.end(); it1++) {
            long long index_here=*it1;
            
            printvectorofvector(index_to_equivalent_sets[index_here]);
            
        }
    }
    
    
    //            debug lines
    for (set<spin_no_periodic_struct>::iterator it1=lowest_energy_spin_struct_set.begin(); it1!=lowest_energy_spin_struct_set.end(); it1++){
        spin_no_periodic_struct spin_struct_temp=*it1;
        if  (obscenely_verbose)
        {
            printblock(spin_struct_temp.spin);
            cout<<"\n has PI value of "<<PI_value[spin_struct_temp]<<"\n";
            double energy_temp=0;
            map<set<tuple<int,int,int,int,int> >, double> cluster_type_temp;
            calculate_cluster_type_and_energy_no_periodic(J_input_best_memory, energy_temp, cluster_type_temp, spin_struct_temp.spin);
            cout<<"has energy value of"<<energy_temp;
        }

        //                printmap(spin_struct_temp.spin);
        assert(PI_value.count(spin_struct_temp)>0);

    }
    
    //            look for prototypes set
    for (set<spin_no_periodic_struct>::iterator it1=lowest_energy_spin_struct_set.begin(); it1!=lowest_energy_spin_struct_set.end(); it1++) {
        
        spin_no_periodic_struct spin_struct_now=*it1;
        
        bool shall_i_do_it=false;
        if (PI_value.count(spin_struct_now)>0) {
            if (fabs(PI_value[spin_struct_now])>1e-5) {
                shall_i_do_it=true;
            }
        }
        
        if (shall_i_do_it==false) {
            continue;
        }
        
        
        
        //extract prototype set here
        {
            
            
            if (x_range>1)
            {
                
                //                        cout<<"\n what is the initial spin\n";
                //                        printblock(spin_struct_now.spin);
                spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_x_pos();
                //                        cout<<"\n what is the truncated spin\n";
                //                        printblock(truncated_initial_spin_struct_now.spin);
                //firstly_enforce_we_do_not_loop_back
                
                
                bool already_in_model=false;
                set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                    
                    tuple<int,int,int,int> initial_tuple=it3->first;
                    int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                    int initial_type=it3->second;
                    prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                }
                set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                    if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                        already_in_model=true;
                    }
                }
                //end of enforce
                
                
                
                
                if (already_in_model==false) {
                    spin_no_periodic_struct to_compare_spin_struct;
                    bool inside_some=false;
                    if (x_range>1) {
                        for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                            
                            spin_no_periodic_struct temp_spin_struct=*it2;
                            
                            to_compare_spin_struct=temp_spin_struct.remove_x_neg();
                            
                            //                                    cout<<"\n temp_spin_struct is";
                            //                                    printblock(temp_spin_struct.spin);
                            //                                    cout<<"\n to_compare_spin_struct is";
                            //                                    printblock(to_compare_spin_struct.spin);
                            
                            if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                inside_some=true;
                                break;
                            }
                            
                        }
                    }
                    
                    if (obscenely_verbose) {
                        cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                        printblock(truncated_initial_spin_struct_now.spin);
                        cout<<"\nwhat is to_compare_spin_struct";
                        printblock(to_compare_spin_struct.spin);
                        cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                    }

                    
                    
                    if (inside_some==false) {
                        for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                            
                            tuple<int,int,int,int> initial_tuple=it3->first;
                            int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                            int initial_type=it3->second;
                            prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                        }
                        
                        break;
                    }
                }
                
                
            }
            if (x_range>1)
            {
                spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_x_neg();
                
                //firstly_enforce_we_do_not_loop_back
                bool already_in_model=false;
                set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                    
                    tuple<int,int,int,int> initial_tuple=it3->first;
                    int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                    int initial_type=it3->second;
                    prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                }
                set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                    if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                        already_in_model=true;
                    }
                }
                //end of enforce
                
                if (already_in_model==false) {
                    spin_no_periodic_struct to_compare_spin_struct;
                    bool inside_some=false;
                    if (x_range>1) {
                        for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                            
                            spin_no_periodic_struct temp_spin_struct=*it2;
                            to_compare_spin_struct=temp_spin_struct.remove_x_pos();
                            
                            if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                inside_some=true;
                                break;
                            }
                            
                        }
                    }
                    
                    if (obscenely_verbose) {
                        cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                        printblock(truncated_initial_spin_struct_now.spin);
                        cout<<"\nwhat is to_compare_spin_struct";
                        printblock(to_compare_spin_struct.spin);
                        cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                    }

                    
                    
                    if (inside_some==false) {
                        for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                            
                            tuple<int,int,int,int> initial_tuple=it3->first;
                            int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                            int initial_type=it3->second;
                            prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                        }
                        
                        break;
                    }
                }
                
                
            }
            if (y_range>1)
            {
                spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_y_pos();
                
                //firstly_enforce_we_do_not_loop_back
                bool already_in_model=false;
                set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                    
                    tuple<int,int,int,int> initial_tuple=it3->first;
                    int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                    int initial_type=it3->second;
                    prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                }
                set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                    if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                        already_in_model=true;
                    }
                }
                //end of enforce
                
                if (already_in_model==false) {
                    spin_no_periodic_struct to_compare_spin_struct;
                    bool inside_some=false;
                    if (y_range>1) {
                        for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                            
                            spin_no_periodic_struct temp_spin_struct=*it2;
                            to_compare_spin_struct=temp_spin_struct.remove_y_neg();
                            
                            if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                inside_some=true;
                                break;
                            }
                            
                        }
                    }
                    
                    
                    //                        cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                    //                        printblock(truncated_initial_spin_struct_now.spin);
                    //                        cout<<"\nwhat is to_compare_spin_struct";
                    //                        printblock(to_compare_spin_struct.spin);
                    //                        cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                    
                    
                    if (inside_some==false) {
                        for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                            
                            tuple<int,int,int,int> initial_tuple=it3->first;
                            int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                            int initial_type=it3->second;
                            prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                        }
                        
                        break;
                    }
                }
                
                
            }
            if (y_range>1)
            {
                spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_y_neg();
                
                //firstly_enforce_we_do_not_loop_back
                bool already_in_model=false;
                set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                    
                    tuple<int,int,int,int> initial_tuple=it3->first;
                    int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                    int initial_type=it3->second;
                    prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                }
                set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                    if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                        already_in_model=true;
                    }
                }
                //end of enforce
                
                if (already_in_model==false) {
                    spin_no_periodic_struct to_compare_spin_struct;
                    bool inside_some=false;
                    if (y_range>1) {
                        for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                            
                            spin_no_periodic_struct temp_spin_struct=*it2;
                            to_compare_spin_struct=temp_spin_struct.remove_y_pos();
                            
                            if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                inside_some=true;
                                break;
                            }
                            
                        }
                    }
                    
                    
                    //                        cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                    //                        printblock(truncated_initial_spin_struct_now.spin);
                    //                        cout<<"\nwhat is to_compare_spin_struct";
                    //                        printblock(to_compare_spin_struct.spin);
                    //                        cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                    
                    
                    if (inside_some==false) {
                        for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                            
                            tuple<int,int,int,int> initial_tuple=it3->first;
                            int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                            int initial_type=it3->second;
                            prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                        }
                        
                        break;
                    }
                }
                
                
            }
            if (z_range>1)
            {
                spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_z_pos();
                
                //firstly_enforce_we_do_not_loop_back
                bool already_in_model=false;
                set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                    
                    tuple<int,int,int,int> initial_tuple=it3->first;
                    int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                    int initial_type=it3->second;
                    prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                }
                set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                    if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                        already_in_model=true;
                    }
                }
                //end of enforce
                
                if (already_in_model==false) {
                    spin_no_periodic_struct to_compare_spin_struct;
                    bool inside_some=false;
                    if (z_range>1) {
                        for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                            
                            spin_no_periodic_struct temp_spin_struct=*it2;
                            to_compare_spin_struct=temp_spin_struct.remove_z_neg();
                            
                            if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                inside_some=true;
                                break;
                            }
                            
                        }
                    }
                    
                    
                    //                        cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                    //                        printblock(truncated_initial_spin_struct_now.spin);
                    //                        cout<<"\nwhat is to_compare_spin_struct";
                    //                        printblock(to_compare_spin_struct.spin);
                    //                        cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                    
                    
                    if (inside_some==false) {
                        for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                            
                            tuple<int,int,int,int> initial_tuple=it3->first;
                            int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                            int initial_type=it3->second;
                            prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                        }
                        
                        break;
                    }
                }
                
                
            }
            if (z_range>1)
            {
                spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_z_neg();
                
                //firstly_enforce_we_do_not_loop_back
                bool already_in_model=false;
                set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                    
                    tuple<int,int,int,int> initial_tuple=it3->first;
                    int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                    int initial_type=it3->second;
                    prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                }
                set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                    if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                        already_in_model=true;
                    }
                }
                //end of enforce
                
                if (already_in_model==false) {
                    spin_no_periodic_struct to_compare_spin_struct;
                    bool inside_some=false;
                    if (z_range>1) {
                        for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                            
                            spin_no_periodic_struct temp_spin_struct=*it2;
                            to_compare_spin_struct=temp_spin_struct.remove_z_pos();
                            
                            if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                inside_some=true;
                                break;
                            }
                            
                        }
                    }
                    
                    
                    //                        cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                    //                        printblock(truncated_initial_spin_struct_now.spin);
                    //                        cout<<"\nwhat is to_compare_spin_struct";
                    //                        printblock(to_compare_spin_struct.spin);
                    //                        cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                    
                    
                    if (inside_some==false) {
                        for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                            
                            tuple<int,int,int,int> initial_tuple=it3->first;
                            int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                            int initial_type=it3->second;
                            prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                        }
                        
                        break;
                    }
                }
                
                
            }
        }
        
    }
    
    //            cout<<"\n prototype_set_here.size() is "<<prototype_set_here.size();
    //            cout<<"\n debug here 9561723412, what is the prototype set";
    //            printvector(prototype_set_here);
    //            cout<<"\n\n";
    
    
    
    if (prototype_set_here.size()==0) {
        for (set<spin_no_periodic_struct>::iterator it1=lowest_energy_spin_struct_set.begin(); it1!=lowest_energy_spin_struct_set.end(); it1++) {
            
            spin_no_periodic_struct spin_struct_now=*it1;
            
            //extract prototype set here
            {
                if (x_range>1)
                {
                    spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_x_pos();
                    
                    //firstly_enforce_we_do_not_loop_back
                    bool already_in_model=false;
                    set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                        
                        tuple<int,int,int,int> initial_tuple=it3->first;
                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                        int initial_type=it3->second;
                        prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                    }
                    set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                    construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                    if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                        if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                            already_in_model=true;
                        }
                    }
                    //end of enforce
                    
                    if (already_in_model==false) {
                        spin_no_periodic_struct to_compare_spin_struct;
                        bool inside_some=false;
                        if (x_range>1) {
                            for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                                
                                spin_no_periodic_struct temp_spin_struct=*it2;
                                to_compare_spin_struct=temp_spin_struct.remove_x_neg();
                                
                                if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                    inside_some=true;
                                    break;
                                }
                                
                            }
                        }
                        
                        if (obscenely_verbose) {
                            cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                            printblock(truncated_initial_spin_struct_now.spin);
                            cout<<"\nwhat is to_compare_spin_struct";
                            printblock(to_compare_spin_struct.spin);
                            cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                        }

                        
                        
                        if (inside_some==false) {
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            
                            break;
                        }
                    }
                    
                    
                }
                if (x_range>1)
                {
                    spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_x_neg();
                    
                    //firstly_enforce_we_do_not_loop_back
                    bool already_in_model=false;
                    set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                        
                        tuple<int,int,int,int> initial_tuple=it3->first;
                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                        int initial_type=it3->second;
                        prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                    }
                    set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                    construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                    if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                        if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                            already_in_model=true;
                        }
                    }
                    //end of enforce
                    
                    if (already_in_model==false) {
                        spin_no_periodic_struct to_compare_spin_struct;
                        bool inside_some=false;
                        if (x_range>1) {
                            for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                                
                                spin_no_periodic_struct temp_spin_struct=*it2;
                                to_compare_spin_struct=temp_spin_struct.remove_x_pos();
                                
                                if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                    inside_some=true;
                                    break;
                                }
                                
                            }
                        }
                        
                        if (obscenely_verbose) {
                            cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                            printblock(truncated_initial_spin_struct_now.spin);
                            cout<<"\nwhat is to_compare_spin_struct";
                            printblock(to_compare_spin_struct.spin);
                            cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                        }
                        

                        
                        
                        if (inside_some==false) {
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            
                            break;
                        }
                    }
                    
                    
                }
                if (y_range>1)
                {
                    spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_y_pos();
                    
                    //firstly_enforce_we_do_not_loop_back
                    bool already_in_model=false;
                    set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                        
                        tuple<int,int,int,int> initial_tuple=it3->first;
                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                        int initial_type=it3->second;
                        prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                    }
                    set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                    construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                    if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                        if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                            already_in_model=true;
                        }
                    }
                    //end of enforce
                    
                    if (already_in_model==false) {
                        spin_no_periodic_struct to_compare_spin_struct;
                        bool inside_some=false;
                        if (y_range>1) {
                            for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                                
                                spin_no_periodic_struct temp_spin_struct=*it2;
                                to_compare_spin_struct=temp_spin_struct.remove_y_neg();
                                
                                if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                    inside_some=true;
                                    break;
                                }
                                
                            }
                        }
                        
                        
                        //                        cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                        //                        printblock(truncated_initial_spin_struct_now.spin);
                        //                        cout<<"\nwhat is to_compare_spin_struct";
                        //                        printblock(to_compare_spin_struct.spin);
                        //                        cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                        
                        
                        if (inside_some==false) {
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            
                            break;
                        }
                    }
                    
                    
                }
                if (y_range>1)
                {
                    spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_y_neg();
                    
                    //firstly_enforce_we_do_not_loop_back
                    bool already_in_model=false;
                    set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                        
                        tuple<int,int,int,int> initial_tuple=it3->first;
                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                        int initial_type=it3->second;
                        prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                    }
                    set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                    construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                    if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                        if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                            already_in_model=true;
                        }
                    }
                    //end of enforce
                    
                    if (already_in_model==false) {
                        spin_no_periodic_struct to_compare_spin_struct;
                        bool inside_some=false;
                        if (y_range>1) {
                            for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                                
                                spin_no_periodic_struct temp_spin_struct=*it2;
                                to_compare_spin_struct=temp_spin_struct.remove_y_pos();
                                
                                if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                    inside_some=true;
                                    break;
                                }
                                
                            }
                        }
                        
                        
                        //                        cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                        //                        printblock(truncated_initial_spin_struct_now.spin);
                        //                        cout<<"\nwhat is to_compare_spin_struct";
                        //                        printblock(to_compare_spin_struct.spin);
                        //                        cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                        
                        
                        if (inside_some==false) {
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            
                            break;
                        }
                    }
                    
                    
                }
                if (z_range>1)
                {
                    spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_z_pos();
                    
                    //firstly_enforce_we_do_not_loop_back
                    bool already_in_model=false;
                    set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                        
                        tuple<int,int,int,int> initial_tuple=it3->first;
                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                        int initial_type=it3->second;
                        prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                    }
                    set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                    construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                    if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                        if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                            already_in_model=true;
                        }
                    }
                    //end of enforce
                    
                    if (already_in_model==false) {
                        spin_no_periodic_struct to_compare_spin_struct;
                        bool inside_some=false;
                        if (z_range>1) {
                            for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                                
                                spin_no_periodic_struct temp_spin_struct=*it2;
                                to_compare_spin_struct=temp_spin_struct.remove_z_neg();
                                
                                if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                    inside_some=true;
                                    break;
                                }
                                
                            }
                        }
                        
                        
                        //                        cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                        //                        printblock(truncated_initial_spin_struct_now.spin);
                        //                        cout<<"\nwhat is to_compare_spin_struct";
                        //                        printblock(to_compare_spin_struct.spin);
                        //                        cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                        
                        
                        if (inside_some==false) {
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            
                            break;
                        }
                    }
                    
                    
                }
                if (z_range>1)
                {
                    spin_no_periodic_struct truncated_initial_spin_struct_now=spin_struct_now.remove_z_neg();
                    
                    //firstly_enforce_we_do_not_loop_back
                    bool already_in_model=false;
                    set<tuple<int,int,int,int,int> > prototype_set_here_judge_whether_have_already;
                    for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                        
                        tuple<int,int,int,int> initial_tuple=it3->first;
                        int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                        int initial_type=it3->second;
                        prototype_set_here_judge_whether_have_already.insert(make_tuple(x,y,z,p,initial_type));
                    }
                    set<set<tuple<int,int,int,int,int> > > set_of_equivalent_test;
                    construct_set_of_equivalent_from_prototype_set(set_of_equivalent_test, prototype_set_here_judge_whether_have_already, x_range, y_range, z_range);
                    if (revert_index_to_equivalent_sets.count(set_of_equivalent_test)>0.1) {
                        if (significant_index.count( revert_index_to_equivalent_sets[set_of_equivalent_test])>0.1) {
                            already_in_model=true;
                        }
                    }
                    //end of enforce
                    
                    if (already_in_model==false) {
                        spin_no_periodic_struct to_compare_spin_struct;
                        bool inside_some=false;
                        if (z_range>1) {
                            for (set<spin_no_periodic_struct>::iterator it2=lowest_energy_spin_struct_set.begin(); it2!=lowest_energy_spin_struct_set.end(); it2++) {
                                
                                spin_no_periodic_struct temp_spin_struct=*it2;
                                to_compare_spin_struct=temp_spin_struct.remove_z_pos();
                                
                                if (truncated_initial_spin_struct_now==to_compare_spin_struct) {
                                    inside_some=true;
                                    break;
                                }
                                
                            }
                        }
                        
                        
                        //                        cout<<"\n debug here 84612381 what is truncated_initial_spin_struct_now "<<endl;
                        //                        printblock(truncated_initial_spin_struct_now.spin);
                        //                        cout<<"\nwhat is to_compare_spin_struct";
                        //                        printblock(to_compare_spin_struct.spin);
                        //                        cout<<"\nlet's check whether inside_some is true or false"<<inside_some;
                        
                        
                        if (inside_some==false) {
                            for (map<tuple<int,int,int,int>,int>::iterator it3=truncated_initial_spin_struct_now.spin.begin(); it3!=truncated_initial_spin_struct_now.spin.end(); it3++) {
                                
                                tuple<int,int,int,int> initial_tuple=it3->first;
                                int x=initial_tuple.get<0>(),y=initial_tuple.get<1>(),z=initial_tuple.get<2>(),p=initial_tuple.get<3>();
                                int initial_type=it3->second;
                                prototype_set_here.insert(make_tuple(x,y,z,p,initial_type));
                            }
                            
                            break;
                        }
                    }
                    
                    
                }
            }
            
        }
    }
    
}






void dedicated1D_alg(map<set<tuple<int,int,int,int,int> >, double> J, int x_range,map<int, int> component ,spin_periodic_struct &spin_struct_min,string id, solver_variable global_parameters, double &lowerbound_from_compat,bool obscenely_verbose)
{
    GRBEnv env_graph = GRBEnv();
    if (!obscenely_verbose) {
        env_graph.set(GRB_IntParam_OutputFlag, 0);
    }
    GRBModel model_graph  = GRBModel(env_graph);
    model_graph.set(GRB_IntAttr_ModelSense, 1);
    
    vector<spin_no_periodic_struct> vector_of_all_local_spin_struct;
    
    long long total_number=1;
    for (map<int, int>::iterator it1=component.begin(); it1!=component.end(); it1++) {
        int p=it1->first;
        int comp_num=it1->second;
        total_number*=myPow(comp_num, x_range);
    }
    
    for (long long spin_struct_number=0; spin_struct_number<total_number; spin_struct_number++) {
        map<tuple<int,int,int,int>,int> spin_temp_now;
        for (int i=1; i<=x_range; i++) {
            for (map<int, int>::iterator it1=component.begin(); it1!=component.end(); it1++) {
                
                int p=it1->first;
                int modulor=it1->second;

                long long divider=1;
                for (map<int, int>::iterator it2=component.begin(); it2!=component.end(); it2++) {
                    
                    int p_temp=it2->first;
                    int comp_num=it2->second;
                    if (p<p_temp)
                        divider*=myPow(comp_num, i);
                    else
                        divider*=myPow(comp_num, i-1);
                    
                    
                }
                spin_temp_now[make_tuple(i,1,1,p)]= (spin_struct_number/divider)%modulor;
                
            }
        }
        vector_of_all_local_spin_struct.push_back(spin_temp_now);
        
    }
    
//    cout<<"\n debug to check whether all spin_struct is added";
//
    map<spin_no_periodic_struct,GRBVar> spin_struct_variable;
    map<spin_no_periodic_struct, set<spin_no_periodic_struct> > left_truncated_to_original;
    map<spin_no_periodic_struct, set<spin_no_periodic_struct> > right_truncated_to_original;
    
    
    for (long long i=0; i<vector_of_all_local_spin_struct.size();i++ ) {
//        printblock(vector_of_all_local_spin_struct[i].spin);
        double energy_temp=0;
        map<set<tuple<int, int, int, int, int> >, double> cluster_type_temp;
        spin_no_periodic_struct spin_struct_temp=vector_of_all_local_spin_struct[i];
        map<tuple<int, int, int, int>, int>spin_now=vector_of_all_local_spin_struct[i].spin;
        
        calculate_cluster_type_and_energy_no_periodic(J, energy_temp, cluster_type_temp,  spin_now);
        
        spin_struct_variable[spin_now]=model_graph.addVar(0, 1, energy_temp, GRB_CONTINUOUS);
        
        spin_no_periodic_struct left_truncated=spin_struct_temp.remove_x_neg();
        spin_no_periodic_struct right_truncated=spin_struct_temp.remove_x_pos();
        
        left_truncated_to_original[left_truncated].insert(spin_struct_temp);
        right_truncated_to_original[right_truncated].insert(spin_struct_temp);
    }
    
    model_graph.update();
    GRBLinExpr total=0;
    for (long long i=0; i<vector_of_all_local_spin_struct.size();i++ ) {
        total+=spin_struct_variable[vector_of_all_local_spin_struct[i]];
    }
    model_graph.addConstr(total==1);
    
    for (map<spin_no_periodic_struct, set<spin_no_periodic_struct> >::iterator it1=left_truncated_to_original.begin(); it1!=left_truncated_to_original.end(); it1++) {
        spin_no_periodic_struct truncated_spin=it1->first;
        set<spin_no_periodic_struct> left_extension=it1->second;
        set<spin_no_periodic_struct> right_extension=right_truncated_to_original[truncated_spin];
        GRBLinExpr lhs=0,rhs=0;
        for (set<spin_no_periodic_struct>::iterator it1=left_extension.begin(); it1!=left_extension.end(); it1++) {
            lhs+=spin_struct_variable[*it1];
        }
        for (set<spin_no_periodic_struct>::iterator it1=right_extension.begin(); it1!=right_extension.end() ;it1++ ) {
            rhs+=spin_struct_variable[*it1];
        }
        model_graph.addConstr(lhs==rhs);
    }
    
    model_graph.optimize();
    
    lowerbound_from_compat=model_graph.get(GRB_DoubleAttr_ObjVal);
    
    set<spin_no_periodic_struct> positive_spin_struct;
    for (map<spin_no_periodic_struct, GRBVar>::iterator it1=spin_struct_variable.begin(); it1!=spin_struct_variable.end(); it1++) {
        spin_no_periodic_struct spin_temp=it1->first;
        GRBVar variable_now=it1->second;
        if (variable_now.get(GRB_DoubleAttr_X)>1e-5) {
            positive_spin_struct.insert(spin_temp);
        }
    }
    
//    for (set<spin_no_periodic_struct>::iterator it1=positive_spin_struct.begin(); it1!=positive_spin_struct.end(); it1++) {
//        spin_no_periodic_struct spin_struct_now=*it1;
//        printblock(spin_struct_now.spin);
//    }

    int overlapping_pointer=0;
    vector<spin_no_periodic_struct> final_spin_struct_representation;
    final_spin_struct_representation.push_back(*positive_spin_struct.begin());
    bool stop_signal=false;
    while (stop_signal==false) {
        
        spin_no_periodic_struct last_spin=final_spin_struct_representation[final_spin_struct_representation.size()-1];
        spin_no_periodic_struct left_truncated_last_spin=last_spin.remove_x_neg();
        int initialsize=final_spin_struct_representation.size();
        int final_size=initialsize;
        
        for (set<spin_no_periodic_struct>::iterator it1=positive_spin_struct.begin(); it1!=positive_spin_struct.end(); it1++) {
            
            spin_no_periodic_struct spin_to_examine=*it1;
            spin_no_periodic_struct truncated_right_spin_to_examine=spin_to_examine.remove_x_pos();
            if (left_truncated_last_spin==truncated_right_spin_to_examine) {
                final_spin_struct_representation.push_back(spin_to_examine);
                final_size=final_spin_struct_representation.size();
                //check over lappings
                for (int i=0; i<final_spin_struct_representation.size()-1; i++) {
                    if (final_spin_struct_representation[i]==spin_to_examine) {
                        overlapping_pointer=i;
                        stop_signal=true;
                    }
                }
                break;
            }
        }
        assert(final_size==initialsize+1);
        
        
        
    }
    
    
    if (overlapping_pointer>0) {
        final_spin_struct_representation.erase(final_spin_struct_representation.begin(),final_spin_struct_representation.begin()+overlapping_pointer);
    }
    
    
//    cout<<"\n debug here check overlapping final spin representation ";
//    for (int i=0; i<final_spin_struct_representation.size(); i++) {
//        printblock(final_spin_struct_representation[i].spin);
//    }
    
    tuple<int,int,int,int,int,int> periodicity=make_tuple(final_spin_struct_representation.size()-1,0,1,0,0,1);
    
    map<tuple<int,int,int,int>, int> spin_output;
    
    
    for (int i=0; i<final_spin_struct_representation.size()-1; i++) {
        map<tuple<int,int,int,int>, int> spin_struct_at_i=final_spin_struct_representation[i].spin;
        for (int x=i+1; x<i+1+x_range; x++) {
            for (map<int, int>::iterator it2=component.begin(); it2!=component.end(); it2++) {
                int p=it2->first;
                spin_output[make_tuple(x,1,1,p)]=spin_struct_at_i[make_tuple(x-i,1,1,p)];
            }
        }
    }
    
//    printblock(spin_output);
    
    spin_struct_min.spin=spin_output;
    spin_struct_min.periodicity=periodicity;
}

