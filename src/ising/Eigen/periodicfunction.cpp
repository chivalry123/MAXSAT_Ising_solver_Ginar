#include "periodicfunction.h"

using namespace std;
using namespace ::boost::tuples;
using namespace ::boost;

void periodic(int a0, int a1, int a2, int a3, int a4, int a5,
              map<set<tuple<int,int,int,int,int> >, double> &J,
              int x_range, int y_range, int z_range,
              map<int,int> componentnumber,
              map<tuple<int,int,int,int>,int> &spin,
              double &averageenergy,
              map< tuple<int,int,int,int,int,int>, map<set<tuple<int,int,int,int,int> >, double> > &clustertypeperiodic,
              double min_so_far,
              std::string id,
              bool stop_at_first_iteration,
              bool stop_when_upperbound_smaller_than_predicted_lower_bound,
              bool stop_till_exact_result,
              bool obscenely_verbose,
              solver_variable global_parameters){
    
    int size_of_super_cell=a0*a2*a5;

    if (obscenely_verbose) {
        cout << "Calculating periodic UB with a0 = " << a0;
        cout << ", a1 = " << a1;
        cout << ", a2 = " << a2;
        cout << ", a3 = " << a3;
        cout << ", a4 = " << a4;
        cout << ", a5 = " << a5 << endl;
    }
    
    
//    cout << "in periodic x_range, y_range, z_range "<<make_tuple(x_range,y_range,z_range);

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
        J_integral[it->first]=round(it->second);
    }

    for (int i = 1; i < a0 + 1; i++) {
        for (int j = 1; j < a2 + 1; j++) {
            for (int k = 1; k < a5 + 1; k++) {
                for (map<int, int>::iterator cit = componentnumber.begin(); cit!=componentnumber.end(); cit++) {
                    int position = cit->first;
                    int element_count = cit->second;
                    for (int z = 0; z < element_count; z++) {
                        if (z > 0.1) {
                            softer_variable_encode[make_tuple(i, j, k, position, z)] = softer_count_variable;
                            softer_variable_decode[softer_count_variable] = make_tuple(i, j, k, position, z);
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
                for (map<int, int>::iterator cit=componentnumber.begin(); cit!=componentnumber.end(); cit++) {
                    int position = cit->first;
                    int element_count = cit->second;
                    for (int mu = 0; mu < element_count; mu++) {
                        if (mu > 0.1) {
                            int i0 = positive_modulo(((i-1)-floor_int_division(j-(k-1)/a5*a4,a2)*a1-(k-1)/a5*a3), a0) + 1;
                            int j0 = positive_modulo(((j-1)-(k-1)/a5*a4), a2) + 1;
                            int k0 = (k-1) % a5 + 1;
                            softer_variable_encode[make_tuple(i, j, k, position, mu)] =
                                softer_variable_encode[make_tuple(i0, j0, k0, position, mu)];
                        }
                    }
                }
            }
        }
    }

    for (int i = 1; i < a0 + 1; i++) {
        for (int j = 1; j < a2 + 1; j++) {
            for (int k = 1; k < a5 + 1; k++) {
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

    for (int i = 1; i < a0 + 1;i++) {
        for (int j = 1; j < a2 + 1; j++) {
            for (int k = 1; k < a5 + 1; k++) {
                for(set< set<tuple<int,int,int,int,int> > >::iterator itv = nonzerolist.begin();
                                                                    itv != nonzerolist.end(); ++itv) {
                    tuple<int,int,int,int,int> first, second;
                    set<tuple<int,int,int,int,int> >::iterator itera1 = (*itv).begin();
                    set<tuple<int,int,int,int,int> > tuple_set = *itv;
                    if (J_integral[tuple_set] > 0.1 ) {
                        set<long long> set_of_terms;
                        set<long long> softer_set_of_terms;
                        for(set<tuple<int,int,int,int,int> >::iterator itv1 = tuple_set.begin();
                                                                        itv1 != tuple_set.end(); itv1++) {
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
                                                                    itv1 != tuple_set.end(); itv1++) {
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
                            softer_maxsat_model[softer_set_of_terms] -= J_integral[tuple_set];
                        }
                        overall_constant -= J_integral[tuple_set];
                    }
                }
            }
        }
    }

    ofstream softer_input;
    string softer_filename = id + "_" + to_string_is(a0)+ "." + to_string_is(a1)  + "." + to_string_is(a2)+ "." + to_string_is(a3) + "." + to_string_is(a4)  + "." + to_string_is(a5);
    softer_filename += "_softer_periodic.wcnf";
    softer_input.open(softer_filename.c_str());
    softer_input << "p wcnf " << softer_variable_decode.size() <<
                    " " << softer_maxsat_model.size() + softer_maxsat_model_hard.size() <<
                    " " << softer_max_element;

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

    long long criteria = overall_constant + round(min_so_far*a0*a2*a5) + 100;
    if (stop_at_first_iteration) {
        criteria = 0;
    }else if (stop_when_upperbound_smaller_than_predicted_lower_bound) {
        criteria = (overall_constant + round(min_so_far*a0*a2*a5) ) *1.11;
    }else if (stop_till_exact_result) {
        criteria = 1e17;
    }

    criteria = 1e17;

    
    string criteria_string = lexical_cast<string>(criteria);
    srand (time(NULL));
    long long seed_long=rand() % 1000 + 1;
    string seed_str=lexical_cast<string>(seed_long);
    
    double cut_off_time=1e6;
    
    if (softer_variable_decode.size()>=30) {
        cut_off_time=global_parameters.incomplete_time_cap;
    }
    else{
        cut_off_time=pow(2.0,softer_variable_decode.size())*0.05;
    }
    
    if (cut_off_time>global_parameters.incomplete_time_cap) {
        cut_off_time=global_parameters.incomplete_time_cap;
    }
    
    string cut_off_time_str=lexical_cast<string>(cut_off_time);
    
    
    string result_now;
    
    if (global_parameters.use_incomplete_solver==false||(global_parameters.use_incomplete_solver==true&&size_of_super_cell<global_parameters.N_to_start_incomplete)) {
         result_now = exec("./CCLS_to_akmaxsat " + softer_filename + " " + criteria_string);
    }
    else if (global_parameters.use_incomplete_solver==true&&size_of_super_cell>=global_parameters.N_to_start_incomplete){
        if (global_parameters.which_incomplete_solver=="CCEHC") {
             result_now = exec("./CCEHC -inst " + softer_filename + " -seed " + seed_str + " -t "+cut_off_time_str);
        }
        else if(global_parameters.which_incomplete_solver=="CCLS2015")
        {
            string string_tmp;
            string_tmp="./"+global_parameters.which_incomplete_solver+" -inst " + softer_filename + " -seed " + seed_str + " -t "+cut_off_time_str;
            result_now = exec(string_tmp);
        }
        else if(global_parameters.which_incomplete_solver=="Dist1"||global_parameters.which_incomplete_solver=="Dist2")
        {
            string string_tmp;
            string_tmp="./"+global_parameters.which_incomplete_solver+ " " + softer_filename + " " + cut_off_time_str + seed_str;
            cout<<string_tmp;
            result_now = exec(string_tmp);
        }
        
        
    }
    
    
    std::remove(softer_filename.c_str());
    
    cout<<result_now ;
    
    vector<string> result_pieces;
    split(result_now, '\n', result_pieces);
    string s_line;
    vector<string> v_line;

    long long o_value;
    for (vector<string>::iterator it = result_pieces.begin(); it != result_pieces.end(); it++) {
        string temp_string = *it;
        vector<string> temp_segment;
        split(temp_string, ' ', temp_segment);
        string first_segment = *temp_segment.begin();
        if (first_segment == "s") {
            s_line = temp_string;
        }
        if (first_segment=="o") {
            vector<string>::iterator it1=temp_segment.begin();
            o_value=boost::lexical_cast<long long>(*(++it1));
        }
        if (first_segment == "v") {
            v_line.push_back(temp_string);
            for (vector<string>::iterator it1 = temp_segment.begin(); it1 != temp_segment.end(); it1++) {
                if (it1 != temp_segment.begin()) {
                    long long number_now = boost::lexical_cast<long long>(*it1);
                    if (number_now < -0.1) {
                        soft_result[-number_now] = 0;
                    }else if (number_now > 0.1) {
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

    for (int i = 1; i < a0 + x_range; i++) {
        for (int j = 1; j < a2 + y_range; j++) {
            for (int k = 1; k < a5 + z_range; k++) {
                for (map<int, int>::iterator cit = componentnumber.begin();
                                                            cit != componentnumber.end(); cit++) {
                    int position = cit->first;
                    int element_count = cit->second;
                    for (int mu = 0; mu < element_count; mu++) {
                        int i0 = positive_modulo(((i-1)-floor_int_division(j-(k-1)/a5*a4,a2)*a1-(k-1)/a5*a3), a0) + 1;
                        int j0 = positive_modulo(((j-1)-(k-1)/a5*a4), a2) + 1;
                        int k0 = ((k-1) % a5) + 1;
                        if (s_result[make_tuple(i0, j0, k0, position, mu)] == 1){
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
        cout << "\nUB block is: ";
        printblock(spin);
        cout << endl;
    }

    #pragma omp critical(dataupdate)
    {
        for (set< set<tuple<int,int,int,int,int> > >::iterator it = nonzerolist.begin();
             it != nonzerolist.end(); it++){
            set<tuple<int,int,int,int,int> > thisset = *it;
            clustertypeperiodic[make_tuple(a0, a1, a2, a3, a4, a5)][*it] = 0;
            for (int i = 1; i < a0 + 1; i++) {
                for (int j = 1; j < a2 + 1;j++) {
                    for (int k = 1; k < a5 + 1; k++) {
                        int all_matched_indicator = 0;
                        for (set<tuple<int,int,int,int,int> >::iterator it1 = thisset.begin();
                             it1 != thisset.end(); it1++){
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
                            clustertypeperiodic[make_tuple(a0, a1, a2, a3, a4, a5)][*it] += 1.0/(a0*a2*a5);
                        }
                    }
                }
            }
        }
    }
    



    if (obscenely_verbose){
        cout << "UB configuration: ";
//        printmapfromsets(clustertypeperiodic[make_tuple(a0,a1,a2,a3,a4,a5)]);
        cout << endl;
    }
    
    map<set<tuple<int,int,int,int,int> >, double> cluster_type_here;
    #pragma omp critical(dataupdate)
    {
        cluster_type_here = clustertypeperiodic[make_tuple(a0,a1,a2,a3,a4,a5)];

    }
    averageenergy=0;
    for (map<set<tuple<int,int,int,int,int> >, double>::iterator it = cluster_type_here.begin();
                                                                it != cluster_type_here.end(); it++) {
        averageenergy += J[it->first]*it->second;
    }

    if (soft_result.empty()) {
        averageenergy = 1e10;
    }

    if (obscenely_verbose) {
        cout << "\nPeriodic UB average energy: " << averageenergy;
        cout << endl;
    }
    
//    cout << "\n periodicity is :"<<make_tuple(a0,a1,a2,a3,a4,a5);
//    cout << "\n Debug Periodic UB based on E average energy: " << averageenergy;
//    double UB_based_on_o_is = (o_value - overall_constant)/(a0*a2*a5);
//    cout << "\n UB based on o_value is: " << UB_based_on_o_is;
//    if (fabs(UB_based_on_o_is-averageenergy)>1e5*(abs(UB_based_on_o_is)+abs(averageenergy))) {
//       cout<<"\n alert! UB_based_on_o_is is significantly different from UB";
//    }
    
    
    
    
    softer_input.close();
}

void periodic_slave()
{
    mpi::communicator mpi_world;

//    cout<<"\n hello I am slave:"<<mpi_world.rank()<<endl;
    
//    initialize my slave
    
    std::string id="IS"+to_string(mpi_world.rank());
    int max_sites = 50;
    int num_loops = 4;
    double constant=0;
    double prec = 0.00001;
    bool translation_algorithm = false;
    bool basic_exact_mode = false;
    bool pseudo_mode = true;
    bool pseudo_mode_with_proof = false;
    bool verbose = true;
    bool very_verbose = false;
    bool obscenely_verbose = false;
    bool input_PRIM_output_PRIMOUT_mode=false;
    double limit_dimension=100000;
    bool mu_translation_algorithm=false;
    double mu_constant=0;
    bool work_with_mu=false;
    bool scan_chemical_potential=false;
    bool new_cluster_algorithm=false;
    int use_new_pair_terms=0;
    int use_new_triplet_terms=0;
    bool output_more_states=false;
    bool output_states_below_hull=false;
    double how_much_lower=5e-3;
    double output_states_how_sparse=4e-3;
    bool use_level_method=true;
    bool use_weighted_dual_average=false;
    solver_variable global_parameters;
    global_parameters.use_gradual_introduction_new_variables=false;
    global_parameters.do_monte_carlo=false;
    global_parameters.dedicated1D=false;
    global_parameters.ternary_alg=false;
    global_parameters.model_numbers=50;
    global_parameters.ternary_output_states_above_hull=false;
    global_parameters.ternary_x_min=0;
    global_parameters.ternary_x_max=1;
    global_parameters.ternary_y_min=0;
    global_parameters.ternary_y_max=1;
    global_parameters.ternary_z_min=0;
    global_parameters.ternary_z_max=1;
    global_parameters.ternary_output_states_above_hull_gap=0.001;
    global_parameters.ternary_output_states_above_hull_ceil=0.03;
    global_parameters.ternary_debug=false;
    global_parameters.mu1=0;
    global_parameters.mu2=0;
    
    broadcast(mpi_world,max_sites,0);
    map<set<tuple<int,int,int,int,int> >, double> J,J_in;
    broadcast(mpi_world,J,0);
    broadcast(mpi_world,prec,0);
    broadcast(mpi_world,num_loops,0);
    broadcast(mpi_world,basic_exact_mode,0);
    broadcast(mpi_world,pseudo_mode,0);
    broadcast(mpi_world,pseudo_mode_with_proof,0);
    broadcast(mpi_world,verbose,0);
    broadcast(mpi_world,very_verbose,0);
    broadcast(mpi_world,obscenely_verbose,0);
    broadcast(mpi_world,input_PRIM_output_PRIMOUT_mode,0);
    broadcast(mpi_world,limit_dimension,0);
    broadcast(mpi_world,constant,0);
    
    map< set<tuple<int,int,int,int,int> >, double>  mu;

    broadcast(mpi_world,mu,0);
    broadcast(mpi_world,mu_translation_algorithm,0);
    broadcast(mpi_world,mu_constant,0);
    broadcast(mpi_world,work_with_mu,0);
    broadcast(mpi_world,scan_chemical_potential,0);
    broadcast(mpi_world,new_cluster_algorithm,0);
    broadcast(mpi_world,use_new_pair_terms,0);
    broadcast(mpi_world,use_new_triplet_terms,0);
    broadcast(mpi_world,output_more_states,0);
    broadcast(mpi_world,output_states_below_hull,0);
    broadcast(mpi_world,how_much_lower,0);
    broadcast(mpi_world,output_states_how_sparse,0);
    broadcast(mpi_world,use_level_method,0);
    broadcast(mpi_world,use_weighted_dual_average,0);
    broadcast(mpi_world,global_parameters,0);
    
    
    
//    cout<<"\n hello I am slave: "<<mpi_world.rank()<<" received the global parameters"<<endl;

    
    bool kill_all_slave=false;
    
    broadcast(mpi_world,kill_all_slave,0);
    
    while (kill_all_slave==false) {
        
        
        vector<tuple<int,int,int,int,int,int> > periodicity_vector;
        
//        std::cout << "[SLAVE: " << mpi_world.rank()<< "] I am waiting to receive periodicity vector"<<endl;
        
        periodicity_vector.clear();
        
        broadcast(mpi_world,periodicity_vector,0);
        
        vector<double> mu_vector_numeric;

        
        broadcast(mpi_world,mu_vector_numeric,0);


//        std::cout << "[SLAVE: " << mpi_world.rank()<< "] I am waiting to J "<<endl;
        
//        J.clear();
        
        map<set<tuple<int,int,int,int,int> >, double> J;

        broadcast(mpi_world,J,0);
        
        map<set<tuple<int,int,int,int,int> >, double> J_fixed_part;
        double constant,mu_constant;
        
        broadcast(mpi_world,J_fixed_part,0);
        broadcast(mpi_world,constant,0);
        broadcast(mpi_world,mu_constant,0);
        
        
        
        int x_range,y_range,z_range;
        double min_bound;
        map<int,int> component;
        
//        cout<< "[SLAVE: " << mpi_world.rank()<<"] hello I am slave: "<<mpi_world.rank()<<" received 3rd set of the global parameters"<<endl;
        
        broadcast(mpi_world,x_range,0);
        broadcast(mpi_world,y_range,0);
        
//        cout<< "[SLAVE: " << mpi_world.rank()<<"] hello I am slave: "<<mpi_world.rank()<<" received 4th set of the global parameters"<<endl;
        broadcast(mpi_world,z_range,0);
        broadcast(mpi_world,component,0);
//        cout<< "[SLAVE: " << mpi_world.rank()<<"] hello I am slave: "<<mpi_world.rank()<<" received 5th set of the global parameters"<<endl;
        
        
        broadcast(mpi_world,min_bound,0);
        
//        cout<< "[SLAVE: " << mpi_world.rank()<<"] hello I am slave: "<<mpi_world.rank()<<" received 6th set of the global parameters"<<endl;
        
        bool stop = false;
        string status="initial";
        mpi_world.send(0, 0, status);
//        
//        cout<< "[SLAVE: " << mpi_world.rank()<<"] hello I am slave: "<<mpi_world.rank()<<" I am going to send to root status"<<endl;
        
        
        mpi_world.recv(0, 0, stop);
        
//        cout<< "[SLAVE: " << mpi_world.rank()<<"] hello I am slave: "<<mpi_world.rank()<<" I am going to receive stop signal from root"<<endl;
        
        
        while(!stop) {
            
            map< tuple<int,int,int,int,int,int>, map<set<tuple<int,int,int,int,int> >, double> > clustertype_periodic;
            
            // Wait for new job
            
            unsigned int job_id = 0;
            
//            cout<< "[SLAVE: " << mpi_world.rank()<<"] hello I am slave: "<<mpi_world.rank()<<" I am going to receive job id from root"<<endl;
            
            mpi_world.recv(0, 0, job_id);
            
//            cout<< "[SLAVE: " << mpi_world.rank()<<"] hello I am slave: "<<mpi_world.rank()<<" I am going to receive stop min bound from root"<<endl;
            
            
            mpi_world.recv(0, 1, min_bound);
            
            
//            std::cout << "[SLAVE: " << mpi_world.rank()<< "] Received job " << job_id << " from MASTER.\n"<<endl;
            
            // Perform "job"
            
            map<tuple<int,int,int,int>,int>  spin_tmp;
            double energy_tmp=0;
            
            tuple<int,int,int,int,int,int>  periodicity_now;
            double mu_now;
            {
                int i=job_id;
                if  (global_parameters.activate_large_cell_algo_binary)
                {
                    periodicity_now=periodicity_vector[0];
                    mu_now=mu_vector_numeric[i];
                    J.clear();
                    J_in.clear();
                    
                    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
                        mu[it1->first]=mu_now;
                    }
                    
                    
                    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_fixed_part.begin(); it1!=J_fixed_part.end(); it1++) {
                        J_in[it1->first]=it1->second;
                        if (mu.count(it1->first)==1) {
                            J_in[it1->first]+=mu[it1->first];
                        }
                    }
                    
                    
                    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_in.begin(); it1!=J_in.end(); it1++) {
                        J[it1->first]=J_in[it1->first]/global_parameters.prec;
                    }
                    
//                    do something on J;
//                std::cout << "[SLAVE: " << mpi_world.rank()<< "] my mu is "<< mu_now<<endl;
//                std::cout << "[SLAVE: " << mpi_world.rank()<< "] my J is is "<<endl;
//                    printmapfromsets(J);


                    
                    
                }
                else if (!global_parameters.activate_large_cell_algo_binary)
                {
                    periodicity_now=periodicity_vector[i];
                }
                
//                cout<<"\n periodicity_now is "<<periodicity_now;
                int a0=periodicity_now.get<0>();
                int a1=periodicity_now.get<1>();
                int a2=periodicity_now.get<2>();
                int a3=periodicity_now.get<3>();
                int a4=periodicity_now.get<4>();
                int a5=periodicity_now.get<5>();
                
//                cout << "in periodic x_range, y_range, z_range "<<make_tuple(x_range,y_range,z_range);

                
                {
                    if (pseudo_mode) {
                        
                        
                        
                        periodic(a0, a1, a2, a3, a4, a5,
                                 J, x_range, y_range, z_range, component,
                                 spin_tmp,
                                 energy_tmp,
                                 clustertype_periodic, min_bound, id,
                                 true, false, false,
                                 obscenely_verbose,global_parameters);
                        
                        //                    {
                        //                        spin_periodic[make_tuple(a0,a1,a2,a3,a4,a5)]=spin_tmp;
                        //                        energy_periodic[make_tuple(a0,a1,a2,a3,a4,a5)]=energy_tmp;
                        //                    }
                    }
                    else if (basic_exact_mode){
                        
                        
                        
                        
                        periodic(a0, a1, a2, a3, a4, a5,
                                 J, x_range, y_range, z_range, component,
                                 spin_tmp,
                                 energy_tmp,
                                 clustertype_periodic, min_bound, id,
                                 false, true, false,
                                 obscenely_verbose,global_parameters);
                        
                        //                    {
                        //                        spin_periodic[make_tuple(a0,a1,a2,a3,a4,a5)]=spin_tmp;
                        //                        energy_periodic[make_tuple(a0,a1,a2,a3,a4,a5)]=energy_tmp;
                        //                    }
                        
                        
                    }
                    
                }
            }
            
            
            
            // Notify master that the job is done
            
//            std::cout << "[SLAVE: " << mpi_world.rank()<< "] Done with job " << job_id << ". Notifying MASTER.\n"<<endl;
            
            status="finish";
            
            mpi_world.send(0, 0,status);
//            std::cout << "[SLAVE: " << mpi_world.rank()<< "] I have sent status back to master"<<endl;
            
            //note here, I deliberated exclude cluster_type_periodic which may be the culprit of a lot of bad computing performance for realistic system
            
            //send periodicity, spin_tmp and energy_tmp
            
            //        and remember to send more detail;
            
            
            mpi_world.send(0, 1,spin_tmp);
//            std::cout << "[SLAVE: " << mpi_world.rank()<< "] I have sent spin_tmp back to master"<<endl;
            mpi_world.send(0, 2,energy_tmp);
//            std::cout << "[SLAVE: " << mpi_world.rank()<< "] I have sent energy_tmp back to master"<<endl;
            mpi_world.send(0, 3,periodicity_now);
//            std::cout << "[SLAVE: " << mpi_world.rank()<< "] I have sent periodicity_now back to master"<<endl;
            
            
            if (global_parameters.activate_large_cell_algo_binary) {
                mpi_world.send(0, 6,mu_now);
            }
            
            
            if (!global_parameters.ternary_alg) {
                double formation_energy_out, concentration_out;
                get_concentration_formation_energy_from_spin_and_J_in(periodicity_now, mu, J_fixed_part, spin_tmp, constant, mu_constant, formation_energy_out, concentration_out);
                
                mpi_world.send(0, 4,formation_energy_out);
//                std::cout << "[SLAVE: " << mpi_world.rank()<< "] I have sent formation_energy_out back to master"<<endl;
                
                mpi_world.send(0, 5,concentration_out);
//                std::cout << "[SLAVE: " << mpi_world.rank()<< "] I have sent concentration_out back to master"<<endl;
            }
            
            

            
            
            // Check if a new job is coming
            
            mpi_world.recv(0, 0, stop);
//            std::cout << "[SLAVE: " << mpi_world.rank()<< "] I have receive stop from master and stop is "<<stop<<endl;
            
            
        }
        
        
//        std::cout << "[SLAVE: " << mpi_world.rank()<< "] I am waiting for signal whether kill_all_slave"<<endl;

        
        broadcast(mpi_world,kill_all_slave,0);
        
//        std::cout << "[SLAVE: " << mpi_world.rank()<< "] I receive kill_all_slave signal as "<<kill_all_slave<<endl;


    }

    
    
//    std::cout << "~~~~~~~~ Rank " << mpi_world.rank() << " is exiting ~~~~~~~~~~~\n"<<endl;
    

    
    
    
    
}





