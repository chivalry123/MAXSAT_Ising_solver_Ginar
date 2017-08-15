/* This module defines overall interface to the generalized Ising solver.
 * See solver.h for details
 *
 * Author: Wenxuan Huang
 * Maintainer: Wenxuan Huang, Daniil Kitchaev
 * Date: 15 March, 2015
 *
 * Copyright: Wenxuan Huang (C) 2014, All rights reserved.
 */



#include "common.h"

//boost::mpi::environment mpi_env;
//boost::mpi::communicator mpi_world;

using namespace std;
using namespace ::boost::tuples;
using namespace ::boost;
using Eigen::MatrixXd;
using namespace Eigen;
using namespace ::boost::mpi;


void read_from_file(std::string &id,
                    int &max_sites,
                    map< set<tuple<int,int,int,int,int> >, double> &J,
                    double &prec,
                    int &num_loops,
                    bool &translation_algorithm,
                    bool &basic_exact_mode,
                    bool &pseudo_mode,
                    bool &pseudo_mode_with_proof,
                    bool &verbose,
                    bool &very_verbose,
                    bool &obscenely_verbose,
                    bool &input_PRIM_output_PRIMOUT_mode,
                    double &limit_dimension,
                    double &constant,
                    map< set<tuple<int,int,int,int,int> >, double> &mu,
                    bool &mu_translation_algorithm,
                    double &mu_constant,
                    bool &work_with_mu,
                    bool &scan_chemical_potential,
                    bool &new_cluster_algorithm,
                    int &use_new_pair_terms,
                    int &use_new_triplet_terms,
                    bool &output_more_states,
                    bool &output_states_below_hull,
                    double &how_much_lower,
                    double &output_states_how_sparse,
                    bool &use_level_method,
                    bool &use_weighted_dual_average,
                    solver_variable &global_parameters);




int main (void){
    
    mpi::environment mpi_env;
    mpi::communicator mpi_world;
    
    
    char buffer[8192];
    
    setvbuf(stdout, buffer, _IOFBF, sizeof(buffer));
    
    cout.precision(15);
    cout<<std::scientific;
    std::string id="IS"+to_string(mpi_world.rank());
//    std::string id="IS0";
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
    global_parameters.use_incomplete_solver=false;
    global_parameters.which_incomplete_solver="CCEHC";
    global_parameters.incomplete_time_cap=600;
    global_parameters.N_to_start_incomplete=100;
    global_parameters.Major_periodicity=false;
    global_parameters.activate_large_cell_algo_binary=false;

    double lower_bound, upper_bound, exact_lowerbound,formation_energy_UB,formation_energy_LB;
    map< set<tuple<int,int,int,int,int> >, double> J, lowerboundclustertype, upperboundclustertype, J_for_proof,mu;
    map< tuple<int,int,int,int>,int> unitcell;
    tuple<int,int,int,int,int,int> periodicity;
    map< tuple<int,int,int,int>, int> cellrepresentation;
    map<tuple<int,int,int,int,int,int>, map<tuple<int,int,int,int>,int> > map_periodicity_to_spin;
    

    
    
    if (mpi_world.rank()==0)
    {
        read_from_file(id, max_sites, J, prec, num_loops, translation_algorithm,
                       basic_exact_mode, pseudo_mode, pseudo_mode_with_proof,
                       verbose, very_verbose, obscenely_verbose,input_PRIM_output_PRIMOUT_mode,limit_dimension,constant,mu,mu_translation_algorithm,mu_constant,work_with_mu,scan_chemical_potential,new_cluster_algorithm,use_new_pair_terms,use_new_triplet_terms, output_more_states,output_states_below_hull,how_much_lower,output_states_how_sparse,use_level_method,use_weighted_dual_average,global_parameters);
        
    }
    else
    {
        bool restart_signal=true;
        
        periodic_slave();
        return 0;
    }
    
    

    broadcast(mpi_world,max_sites,0);
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
    
    
    
    
    
    
    
//    cout<<"\npositive_modulo(-4, 10)"<<positive_modulo(-4, 10);
    if (!scan_chemical_potential&&!global_parameters.do_monte_carlo&&!global_parameters.ternary_alg&&!global_parameters.ternary_debug&&!global_parameters.activate_large_cell_algo_binary) {
        run_solver(max_sites, J,
                   lowerboundclustertype, upperboundclustertype, cellrepresentation,
                   lower_bound, upper_bound, exact_lowerbound, unitcell, periodicity,
                   J_for_proof,
                   id,mu,mu_constant,work_with_mu,formation_energy_UB,formation_energy_LB,new_cluster_algorithm,use_new_pair_terms,use_new_triplet_terms,map_periodicity_to_spin,use_level_method,use_weighted_dual_average,global_parameters,
                   prec, num_loops,
                   basic_exact_mode, pseudo_mode, pseudo_mode_with_proof,
                   verbose, very_verbose, obscenely_verbose,input_PRIM_output_PRIMOUT_mode,limit_dimension,constant);
        
        cout<<"\n tell me what is upper bound: "<<upper_bound<<" what is lower bound: "<<lower_bound<<" what is exact lower bound "<< exact_lowerbound;

        
    }
    else if(scan_chemical_potential)
    {
        map<int,int> components;
        int x_range, y_range, z_range;
        calculate_range_from_J(J, x_range, y_range, z_range, components);
        vector<map< set<tuple<int,int,int,int,int> >, double> > cluster_type_vector;
        vector<double> formation_energy_vector;
        map< set<tuple<int,int,int,int,int> >, double> cluster_type_temp;
        vector<vector<int> > to_scan_line;
        double formation_energy_temp=0;
        vector<double> concentration_vector;
        map< set<tuple<int,int,int,int,int> >, double> J_fixed_part,J_in;
        vector<map< tuple<int,int,int,int>,int> >unit_cell_vector;
        vector<tuple<int,int,int,int,int,int> >periodicity_vector;
        
        vector<map< tuple<int,int,int,int>,int> > output_more_states_unit_cell_vector;
        vector<tuple<int,int,int,int,int,int> > output_more_states_periodicity_vector;
        vector<double> output_more_states_concentration_vector;
        vector<int> output_more_states_super_cell_size_vector;
        vector<double> output_more_states_formation_energy_vector;
        
        
        

        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J.begin(); it1!=J.end(); it1++) {
            if (mu.count(it1->first)==0) {
                J_fixed_part[it1->first]=it1->second;
            }
            else{
                J_fixed_part[it1->first]=J[it1->first]-mu[it1->first];
            }
        }
        
        formation_energy_temp=0;
        cluster_type_temp.clear();
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J.begin(); it1!=J.end(); it1++) {
            cluster_type_temp[it1->first]=0;
        }
        calculate_formation_energy(J, mu, cluster_type_temp, constant, mu_constant, formation_energy_temp);
        formation_energy_vector.push_back(formation_energy_temp);
        cluster_type_vector.push_back(cluster_type_temp);
        concentration_vector.push_back(0);
        map< tuple<int,int,int,int>,int> unit_cell_temp;
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
            set<tuple<int,int,int,int,int> > set_temp=it1->first;
            for (set<tuple<int,int,int,int,int> >::iterator it2=set_temp.begin(); it2!=set_temp.end(); it2++) {
                tuple<int,int,int,int,int> tuple_temp=*it2;
                unit_cell_temp[make_tuple(1,1,1,tuple_temp.get<3>())]=0;
            }
        }
        unit_cell_vector.push_back(unit_cell_temp);
        
        unit_cell_temp.clear();
        periodicity_vector.push_back(make_tuple(1,0,1,0,0,1));
        
        formation_energy_temp=0;
        cluster_type_temp.clear();
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J.begin(); it1!=J.end(); it1++) {
            cluster_type_temp[it1->first]=1;
        }
        calculate_formation_energy(J, mu, cluster_type_temp, constant, mu_constant, formation_energy_temp);
        formation_energy_vector.push_back(formation_energy_temp);
        cluster_type_vector.push_back(cluster_type_temp);
        cluster_type_temp.clear();
        
        
        formation_energy_temp=0;
        concentration_vector.push_back(1);
        unit_cell_temp.clear();
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
            set<tuple<int,int,int,int,int> > set_temp=it1->first;
            for (set<tuple<int,int,int,int,int> >::iterator it2=set_temp.begin(); it2!=set_temp.end(); it2++) {
                tuple<int,int,int,int,int> tuple_temp=*it2;
                unit_cell_temp[make_tuple(1,1,1,tuple_temp.get<3>())]=tuple_temp.get<4>();
            }
        }
        unit_cell_vector.push_back(unit_cell_temp);
        unit_cell_temp.clear();
        periodicity_vector.push_back(make_tuple(1,0,1,0,0,1));
        
        
        vector<int> first_line_to_check;
        first_line_to_check.push_back(0);
        first_line_to_check.push_back(1);
        to_scan_line.push_back(first_line_to_check);
        
        std::system("rm -r GS_solutions");
        std::system("mkdir GS_solutions");
        
        std::system("mkdir GS_solutions/x0");
        print_poscar_out(periodicity_vector[0], unit_cell_vector[0], components);
        std::system("cp POSCAR_OUT GS_solutions/x0/");
        print_poscar_no_vacancy(periodicity_vector[0], unit_cell_vector[0], components);
        std::system("cp POSCAR GS_solutions/x0/");
        
        
        std::system("mkdir GS_solutions/x1");
        print_poscar_out(periodicity_vector[1], unit_cell_vector[1], components);
        std::system("cp POSCAR_OUT GS_solutions/x1/");
        print_poscar_no_vacancy(periodicity_vector[1], unit_cell_vector[1], components);
        std::system("cp POSCAR GS_solutions/x1/");


        
        while (!to_scan_line.empty()) {
            J_in.clear();
            vector<int> the_line_to_check=to_scan_line[0];
            cout<<"\nthe line to chekc is";
            printvector(the_line_to_check);
            cout<<"  corresponding to concentration:"<<concentration_vector[the_line_to_check[0]]<<" "<<concentration_vector[the_line_to_check[1]]<<endl;
            map< set<tuple<int,int,int,int,int> >, double>correlation_difference,correlation_1,correlation_2;
            correlation_1=cluster_type_vector[the_line_to_check[0]];
            correlation_2=cluster_type_vector[the_line_to_check[1]];
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=correlation_1.begin(); it1!=correlation_1.end(); it1++) {
                correlation_difference[it1->first]=correlation_1[it1->first]-correlation_2[it1->first];
            }
            
            double value_temp=0,value_temp_2=0,value_temp_3=0;
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_fixed_part.begin(); it1!=J_fixed_part.end(); it1++) {
                value_temp-=correlation_difference[it1->first]*it1->second;
            }
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
                value_temp_2+=correlation_difference[it1->first];
            }
            
            value_temp_3=value_temp/value_temp_2;
            double chemical_potential_value=value_temp_3;
            
//            cout<<"\nvalue temp 3 here is"<<value_temp_3<<endl<<"\nvalue_temp_2 is"<<value_temp_2<<endl<<"value_temp_3 is"<<value_temp_3<<endl;
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
                mu[it1->first]=value_temp_3;
            }
            
            cout<<"\ntell me what is mu"<<endl;
            printmapfromsets(mu);
            
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_fixed_part.begin(); it1!=J_fixed_part.end(); it1++) {
                J_in[it1->first]=it1->second;
                if (mu.count(it1->first)==1) {
                    J_in[it1->first]+=mu[it1->first];
                }
            }
            
            
            
            run_solver(max_sites, J_in,
                       lowerboundclustertype, upperboundclustertype, cellrepresentation,
                       lower_bound, upper_bound, exact_lowerbound, unitcell, periodicity,
                       J_for_proof,
                       id,mu,mu_constant,work_with_mu,formation_energy_UB,formation_energy_LB,new_cluster_algorithm,use_new_pair_terms,use_new_triplet_terms,map_periodicity_to_spin,use_level_method,use_weighted_dual_average,global_parameters,
                       prec, num_loops,
                       basic_exact_mode, pseudo_mode, pseudo_mode_with_proof,
                       verbose, very_verbose, obscenely_verbose,input_PRIM_output_PRIMOUT_mode,limit_dimension,constant,J_fixed_part);
            

            
            
            if (output_more_states||output_states_below_hull) {
                for (map<tuple<int,int,int,int,int,int>, map<tuple<int,int,int,int>,int> >::iterator it1=map_periodicity_to_spin.begin(); it1!=map_periodicity_to_spin.end(); it1++) {
                    tuple<int,int,int,int,int,int> periodicity_now=it1->first;
                    map<tuple<int,int,int,int>,int>  spin_now=it1->second;
                    
                    double formation_energy_out;
                    double concentration_out;
                    
                    concentration_out=global_parameters.formation_periodic_to_concentration_transfer[periodicity_now];
                    formation_energy_out=global_parameters.formation_periodic_to_energy_transfer[periodicity_now];
//                    get_concentration_formation_energy_from_spin_and_J_in(periodicity_now, mu, J_fixed_part, spin_now, constant, mu_constant, formation_energy_out, concentration_out);
                    
                    bool to_add_this_structure=true;

                    for (int ct1=0;ct1<output_more_states_concentration_vector.size();ct1++)
                    {
                        double concentration_to_compare=output_more_states_concentration_vector[ct1];
                        double concentration_closeness=0;
                        concentration_closeness=fabs(concentration_out-concentration_to_compare);
                        if (concentration_closeness<= 1e-6)
                        {
                            double formation_energy_to_compare=output_more_states_formation_energy_vector[ct1];
                            int super_cell_size_to_compare=output_more_states_super_cell_size_vector[ct1];
                            
                            
                            //check to replace structure edits 109
                            if (formation_energy_out-formation_energy_to_compare>=-1e-7 && formation_energy_out-formation_energy_to_compare<= 1e-7 && periodicity_now.get<0>()*periodicity_now.get<2>()*periodicity_now.get<5>()< super_cell_size_to_compare)
                            {
                                output_more_states_concentration_vector[ct1]=(concentration_out);
                                output_more_states_periodicity_vector[ct1]=(periodicity_now);
                                output_more_states_unit_cell_vector[ct1]=(spin_now);
                                output_more_states_super_cell_size_vector[ct1]=(periodicity_now.get<0>()*periodicity_now.get<2>()*periodicity_now.get<5>());
                                output_more_states_formation_energy_vector[ct1]=(formation_energy_out);
                            }
                            
                            
                            if (formation_energy_out-formation_energy_to_compare>=-1e-5 && formation_energy_out-formation_energy_to_compare<= output_states_how_sparse)
                            {
                                
                                to_add_this_structure=false;
                                break;
//                                if (periodicity_now.get<0>()*periodicity_now.get<2>()*periodicity_now.get<5>()>= super_cell_size_to_compare)
//                                {
//                                to_add_this_structure=false;
//                                break;
//                                }
                            }
                            
                        }
                        
                    }
                    
//                    cout<<"\n formation energy is "<<formation_energy_out<<" concentration out is "<<concentration_out;
                    if (to_add_this_structure)
                    {
                        
                        output_more_states_concentration_vector.push_back(concentration_out);
                        output_more_states_periodicity_vector.push_back(periodicity_now);
                        output_more_states_unit_cell_vector.push_back(spin_now);
                        output_more_states_super_cell_size_vector.push_back(periodicity_now.get<0>()*periodicity_now.get<2>()*periodicity_now.get<5>());
                        output_more_states_formation_energy_vector.push_back(formation_energy_out);
                    }
                    
                    
                }
            }

            
            
            
            
            double concentration_temp=0;
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
                concentration_temp+=upperboundclustertype[it1->first];
            }
            concentration_temp=concentration_temp/mu.size();
            
            cout<<"let me check what is formation energy_UB: "<<formation_energy_UB;
            cout<<".  concentration temp is"<<concentration_temp;
            bool inside_indicator=false;
            for (int i=0; i<concentration_vector.size();i++)
            {
                if (fabs(concentration_vector[i]-concentration_temp)<1e-5) {
                    inside_indicator=true;
                }
            }
            
            if (inside_indicator==false) {
                concentration_vector.push_back(concentration_temp);
                cluster_type_vector.push_back(upperboundclustertype);
                formation_energy_vector.push_back(formation_energy_UB);
                int current_index=concentration_vector.size()-1;
                double lower_concentration=-1,upper_concentration=2;
                int lower_index=-1000,upper_index=1000;
                for (int i=0;i<concentration_vector.size()-1;i++) {
                    if (concentration_vector[i]>lower_concentration&&concentration_vector[i]<concentration_temp) {
                        lower_concentration=concentration_vector[i];
                        lower_index=i;
                    }
                    if (concentration_vector[i]<upper_concentration&&concentration_vector[i]>concentration_temp) {
                        upper_concentration=concentration_vector[i];
                        upper_index=i;
                    }
                }
                vector<int> temp_vector;
                temp_vector.clear();
                temp_vector.push_back(lower_index);
                temp_vector.push_back(current_index);
                to_scan_line.push_back(temp_vector);
                temp_vector.clear();
                temp_vector.push_back(current_index);
                temp_vector.push_back(upper_index);
                to_scan_line.push_back(temp_vector);
                
                unit_cell_vector.push_back(unitcell);
                periodicity_vector.push_back(periodicity);
                
                std::ostringstream ss;
                ss << std::fixed << std::setprecision(6);
                ss << concentration_vector[current_index];
                
                string strings_temp="mkdir GS_solutions/x"+ss.str();
                const char* string_command=strings_temp.c_str();
                std::system(string_command);
                print_poscar_out(periodicity_vector[current_index], unit_cell_vector[current_index], components);
                strings_temp="cp POSCAR_OUT GS_solutions/x"+ss.str()+"/";
                string_command=strings_temp.c_str();
                std::system(string_command);
                
                print_poscar_no_vacancy(periodicity_vector[current_index], unit_cell_vector[current_index], components);
                strings_temp="cp POSCAR GS_solutions/x"+ss.str()+"/";
                string_command=strings_temp.c_str();
                std::system(string_command);
                
            }
            to_scan_line.erase(to_scan_line.begin());
        }
        
        
        map<double,double> concentration_formation_energy;
        map<double,int> concentration_index;
        map<double,int> concentration_super_cell_size;
        map<double,tuple<int,int,int,int,int,int> > concentration_periodicity;

        cout<<"\n after everything showing the results:";
        for (int i=0; i<concentration_vector.size(); i++) {
            concentration_formation_energy[concentration_vector[i]]=formation_energy_vector[i];
            concentration_index[concentration_vector[i]]=i;
            tuple<int,int,int,int,int,int> periodicity_temp;
            periodicity_temp=periodicity_vector[i];
            concentration_super_cell_size[concentration_vector[i]]=periodicity_temp.get<0>()*periodicity_temp.get<2>()*periodicity_temp.get<5>();
            concentration_periodicity[concentration_vector[i]]=periodicity_temp;
        }
        
        stringstream to_print_file;
        to_print_file<<setw(10)<<"x"<<setw(20)<<"formation_energy"<<setw(20)<<"super_cell_size"<<endl;
        for (map<double, double>::iterator it1=concentration_formation_energy.begin(); it1!=concentration_formation_energy.end(); it1++) {
            
            to_print_file<<setw(10)<< std::fixed << std::setprecision(6)<<it1->first<<setw(20)<< std::fixed << std::setprecision(8)<<it1->second<<setw(20)<< std::fixed <<concentration_super_cell_size[it1->first]<<endl;
            
        }
        cout<<"\n"<<to_print_file.str();
        
        const char *path="./GS_solutions/hull.txt";
        std::ofstream file(path);
        file<<to_print_file.str();
        
        
        
        stringstream to_print_file_debug;
        to_print_file_debug<<setw(10)<<"x"<<setw(20)<<"formation_energy"<<setw(20)<<"super_cell_size"<<setw(20)<<"periodicity"<<endl;
        for (map<double, double>::iterator it1=concentration_formation_energy.begin(); it1!=concentration_formation_energy.end(); it1++) {
            
            to_print_file_debug<<setw(10)<< std::fixed << std::setprecision(6)<<it1->first<<setw(20)<< std::fixed << std::setprecision(8)<<it1->second<<setw(20)<< std::fixed <<concentration_super_cell_size[it1->first]<<setw(20)<< std::fixed <<concentration_periodicity[it1->first]<<endl;
            
        }
        cout<<"\n"<<to_print_file_debug.str();
        
        const char *path1="./GS_solutions/hull_debug.txt";
        std::ofstream file1(path1);
        file1<<to_print_file_debug.str();
        
        
        
        if (output_more_states||output_states_below_hull)
        {
            map<double,vector<double> > output_more_states_concentration_formation_energy;
            map<double, vector<int> > output_more_states_concentration_index;
            map<double,vector<int> > output_more_states_concentration_super_cell_size;
            map<double,vector<tuple<int,int,int,int,int,int> > > output_more_states_concentration_periodicity;

            
            {
                vector<map< tuple<int,int,int,int>,int> > output_more_states_unit_cell_vector;
                vector<tuple<int,int,int,int,int,int> > output_more_states_periodicity_vector;
                vector<double> output_more_states_concentration_vector;
                vector<int> output_more_states_super_cell_size_vector;
            }
            
            
            bool still_searching=true;
            
            int number=0;
            
            
            while (still_searching) {
//                still_searching=false;
//                cout<<"\ninside the long loop "<<number;
//                number++;
                
                int index_to_be_added=-1;
                
                
                double maximal_concentration;
                double formation_energy_to_compare_over;
                
//                cout<<"\n output_more_states_concentration_formation_energy.size()>1 ? "<<output_more_states_concentration_formation_energy.size();
                if (output_more_states_concentration_formation_energy.size()>0.1) {
                     maximal_concentration=(output_more_states_concentration_formation_energy.rbegin()->first);
                    vector<double> vector_formation_energy_temp=(output_more_states_concentration_formation_energy.rbegin()->second);
                     formation_energy_to_compare_over=*vector_formation_energy_temp.rbegin();
                }
                else{
                    maximal_concentration=0;
                    formation_energy_to_compare_over=-1e10;
                }
                
//                cout<<"\nwhat is maximal concentration "<<maximal_concentration<<" and formation_energy_to_compare_over "<<formation_energy_to_compare_over;

                bool still_possible_to_add_index=false;
                
                for (int i=0; i<output_more_states_unit_cell_vector.size(); i++) {
                    if (output_more_states_concentration_vector[i]>maximal_concentration-1e-8) {
                        if ((output_more_states_concentration_vector[i]>maximal_concentration+1e-8)||output_more_states_formation_energy_vector[i]>formation_energy_to_compare_over+output_states_how_sparse+1e-8) {
                            still_possible_to_add_index=true;
                            index_to_be_added=i;
                            break;
                        }

                    }
                }
                
                if (!still_possible_to_add_index) {
                    still_searching=false;
                }
                
                if (still_possible_to_add_index) {
                    for (int i=0; i<output_more_states_unit_cell_vector.size(); i++) {
                        if (output_more_states_concentration_vector[i]>maximal_concentration-1e-8) {
                            if ((output_more_states_concentration_vector[i]>maximal_concentration+1e-8)||output_more_states_formation_energy_vector[i]>formation_energy_to_compare_over+output_states_how_sparse+1e-8) {
                                if ((output_more_states_concentration_vector[i]<output_more_states_concentration_vector[index_to_be_added]-1e-8)||(fabs(output_more_states_concentration_vector[i]-output_more_states_concentration_vector[index_to_be_added])<1e-8&&output_more_states_formation_energy_vector[i]<output_more_states_formation_energy_vector[index_to_be_added]-1e-8)||(fabs(output_more_states_concentration_vector[i]-output_more_states_concentration_vector[index_to_be_added])<1e-8&&fabs(output_more_states_formation_energy_vector[i]-output_more_states_formation_energy_vector[index_to_be_added])<1e-8&&output_more_states_super_cell_size_vector[i]<output_more_states_super_cell_size_vector[index_to_be_added])) {
                                    index_to_be_added=i;
                                }
                            }
                        }
                    }
                    
                    
                    int temp_index_to_be_added=index_to_be_added;
                    int min_super_cell_size=1e7;
                    
                    for (int i=0; i<output_more_states_unit_cell_vector.size(); i++) {
                        if ((fabs(output_more_states_concentration_vector[i]-output_more_states_concentration_vector[temp_index_to_be_added])<1e-8&&fabs(output_more_states_formation_energy_vector[i]-output_more_states_formation_energy_vector[temp_index_to_be_added])<output_states_how_sparse/2&&output_more_states_super_cell_size_vector[i]<min_super_cell_size)) {
                            index_to_be_added=i;
                            min_super_cell_size=output_more_states_super_cell_size_vector[i];
                        }
                    }
                    
                    
                    double concentration_to_be_added=output_more_states_concentration_vector[index_to_be_added];
                    int concentration_index_to_be_added=index_to_be_added;
                    int concentration_super_cell_size_to_be_added=output_more_states_super_cell_size_vector[index_to_be_added];
                    double formation_energy_to_be_added=output_more_states_formation_energy_vector[index_to_be_added];
                    tuple<int,int,int,int,int,int> concentration_periodicity_to_be_added=output_more_states_periodicity_vector[index_to_be_added];
                    
                    
                    if (concentration_to_be_added>maximal_concentration+1e-8) {
                        output_more_states_concentration_formation_energy[concentration_to_be_added].push_back(formation_energy_to_be_added);
                        output_more_states_concentration_index[concentration_to_be_added].push_back(index_to_be_added);
                        output_more_states_concentration_super_cell_size[concentration_to_be_added].push_back(concentration_super_cell_size_to_be_added);
                        output_more_states_concentration_periodicity[concentration_to_be_added].push_back(concentration_periodicity_to_be_added);
                        
                    }
                    else {
                        output_more_states_concentration_formation_energy[maximal_concentration].push_back(formation_energy_to_be_added);
                        output_more_states_concentration_index[maximal_concentration].push_back(index_to_be_added);
                        output_more_states_concentration_super_cell_size[maximal_concentration].push_back(concentration_super_cell_size_to_be_added);
                        output_more_states_concentration_periodicity[maximal_concentration].push_back(concentration_periodicity_to_be_added);
                    }
                }

            }
        
//            print out the output_more_states_concentration_super_cell_size
            
            if (output_more_states)
            {
                stringstream output_more_states_string;
                output_more_states_string<<setw(10)<<"x"<<setw(20)<<"formation_energy"<<setw(20)<<"super_cell_size"<<setw(20)<<"corresponding_index"<<endl;
                
                for (map<double,vector<double> >::iterator it1=output_more_states_concentration_formation_energy.begin(); it1!=output_more_states_concentration_formation_energy.end(); it1++) {
                    
                    double concentration_now=it1->first;
                    vector<double> energy_vector=it1->second;
                    int size_to_enumerate_over=energy_vector.size();
                    
                    for (int i=0; i<size_to_enumerate_over; i++) {
                        
                        output_more_states_string<<setw(10)<< std::fixed << std::setprecision(6)<<concentration_now<<setw(20)<< std::fixed << std::setprecision(8)<<energy_vector[i]<<setw(20)<< std::fixed <<output_more_states_concentration_super_cell_size[concentration_now][i]<<setw(20)<<output_more_states_concentration_index[concentration_now][i]<<endl;
                    }
                    
                }
                
                cout<<output_more_states_string.str();
                
                
                std::system("rm -r more_low_E_states");
                std::system("mkdir more_low_E_states");
                
                const char *path="./more_low_E_states/states_info.txt";
                std::ofstream file(path);
                file<<output_more_states_string.str();
                
                for (map<double,vector<double> >::iterator it1=output_more_states_concentration_formation_energy.begin(); it1!=output_more_states_concentration_formation_energy.end(); it1++) {
                    
                    double concentration_now=it1->first;
                    vector<double> energy_vector=it1->second;
                    int size_to_enumerate_over=energy_vector.size();
                    
                    for (int i=0; i<size_to_enumerate_over; i++) {
                        
                        string command_string="mkdir more_low_E_states/index"+lexical_cast<string>(output_more_states_concentration_index[concentration_now][i]);
                        std::system(command_string.c_str());
                        int index=output_more_states_concentration_index[concentration_now][i];
                        print_poscar_out(output_more_states_periodicity_vector[index], output_more_states_unit_cell_vector[index], components);
                        command_string="cp POSCAR_OUT more_low_E_states/index"+lexical_cast<string>(output_more_states_concentration_index[concentration_now][i]);
                        std::system(command_string.c_str());
                        
                    }
                    
                }

            }
        
            
            if (output_states_below_hull) {
//                firstly read the hull;
                map<double,double> hull_map;
                std::ifstream t1("input_hull.in");
                std::stringstream buffer1;
                buffer1 << t1.rdbuf();
                string input_hull_file=buffer1.str();
                
                
                vector<string> input_hull_line;
                split(input_hull_file, '\n', input_hull_line);
                
                for (int i=1; i<input_hull_line.size(); i++) {
                    string this_line=input_hull_line[i];
                    vector<string> this_line_vector;
                    
                    split(this_line, ' ', this_line_vector);
                    
                    double concentration_now= lexical_cast<double>(this_line_vector[0]) ;
                    double formation_enery_now=lexical_cast<double>(this_line_vector[1]);
                    hull_map[concentration_now]=formation_enery_now;
                    
                }
                
//                cout<<"\nwhat is hull_map";
//                printmap(hull_map);
                

                
                map<double,vector<double> > output_states_below_hull_concentration_formation_energy_map;
                map<double, vector<int> > output_states_below_hull_concentration_index_map;
                map<double,vector<int> > output_states_below_hull_concentration_super_cell_size_map;
                map<double,vector<tuple<int,int,int,int,int,int> > > output_states_below_hull_concentration_periodicity_map;
                
                for (map<double,vector<double> >::iterator it1=output_more_states_concentration_formation_energy.begin(); it1!=output_more_states_concentration_formation_energy.end(); it1++) {
                    
                    double concentration_now=it1->first;
                    vector<double> energy_vector=it1->second;
                    int size_to_enumerate_over=energy_vector.size();
                    double energy_at_hull=1e80;
                    evaluate_energy_at_hull(concentration_now,hull_map,energy_at_hull);
                    
                    for (int i=0; i<size_to_enumerate_over; i++) {
                        
                        if (energy_vector[i]<energy_at_hull-how_much_lower) {
                            output_states_below_hull_concentration_formation_energy_map[concentration_now].push_back(energy_vector[i]);
                            output_states_below_hull_concentration_index_map[concentration_now].push_back(output_more_states_concentration_index[concentration_now][i]);
                            output_states_below_hull_concentration_super_cell_size_map[concentration_now].push_back(output_more_states_concentration_super_cell_size[concentration_now][i]);
                            output_states_below_hull_concentration_periodicity_map[concentration_now].push_back(output_more_states_concentration_periodicity[concentration_now][i]);
                        }
                        
                    }
                }
                
                
                
                stringstream output_states_below_hull_string;
                output_states_below_hull_string<<setw(10)<<"x"<<setw(20)<<"formation_energy"<<setw(20)<<"super_cell_size"<<setw(20)<<"corresponding_index"<<endl;
                
                for (map<double,vector<double> >::iterator it1=output_states_below_hull_concentration_formation_energy_map.begin(); it1!=output_states_below_hull_concentration_formation_energy_map.end(); it1++) {
                    
                    double concentration_now=it1->first;
                    vector<double> energy_vector=it1->second;
                    int size_to_enumerate_over=energy_vector.size();
                    
                    for (int i=0; i<size_to_enumerate_over; i++) {
                        
                        output_states_below_hull_string<<setw(10)<< std::fixed << std::setprecision(6)<<concentration_now<<setw(20)<< std::fixed << std::setprecision(8)<<energy_vector[i]<<setw(20)<< std::fixed <<output_states_below_hull_concentration_super_cell_size_map[concentration_now][i]<<setw(20)<<output_states_below_hull_concentration_index_map[concentration_now][i]<<endl;
                    }
                    
                }
                
                cout<<output_states_below_hull_string.str();
                
                
                std::system("rm -r below_the_hull_states");
                std::system("mkdir below_the_hull_states");
                
                const char *path="./below_the_hull_states/below_hull_states_info.txt";
                std::ofstream file(path);
                file<<output_states_below_hull_string.str();
                
                
                
                
                
                //debug mode:
                
                stringstream output_states_below_hull_string_debug;
                output_states_below_hull_string_debug<<setw(10)<<"x"<<setw(20)<<"formation_energy"<<setw(20)<<"super_cell_size"<<setw(20)<<"corresponding_index"<<setw(20)<<"periodicity"<<endl;
                
                for (map<double,vector<double> >::iterator it1=output_states_below_hull_concentration_formation_energy_map.begin(); it1!=output_states_below_hull_concentration_formation_energy_map.end(); it1++) {
                    
                    double concentration_now=it1->first;
                    vector<double> energy_vector=it1->second;
                    int size_to_enumerate_over=energy_vector.size();
                    
                    for (int i=0; i<size_to_enumerate_over; i++) {
                        
                        output_states_below_hull_string_debug<<setw(10)<< std::fixed << std::setprecision(6)<<concentration_now<<setw(20)<< std::fixed << std::setprecision(8)<<energy_vector[i]<<setw(20)<< std::fixed <<output_states_below_hull_concentration_super_cell_size_map[concentration_now][i]<<setw(20)<<output_states_below_hull_concentration_index_map[concentration_now][i]<<setw(20)<<output_states_below_hull_concentration_periodicity_map[concentration_now][i]<<endl;
                    }
                    
                }
                
                cout<<output_states_below_hull_string_debug.str();
                
                
                
                const char *path_debug="./below_the_hull_states/below_hull_states_info_debug.txt";
                std::ofstream file_debug(path_debug);
                file_debug<<output_states_below_hull_string_debug.str();
                //end of debug mode
                
                
                
                
                
                
                
                
                for (map<double,vector<double> >::iterator it1=output_states_below_hull_concentration_formation_energy_map.begin(); it1!=output_states_below_hull_concentration_formation_energy_map.end(); it1++) {
                    
                    double concentration_now=it1->first;
                    vector<double> energy_vector=it1->second;
                    int size_to_enumerate_over=energy_vector.size();
                    
                    for (int i=0; i<size_to_enumerate_over; i++) {
                        
                        string command_string="mkdir below_the_hull_states/index"+lexical_cast<string>(output_states_below_hull_concentration_index_map[concentration_now][i]);
                        std::system(command_string.c_str());
                        int index=output_states_below_hull_concentration_index_map[concentration_now][i];
                        print_poscar_out(output_more_states_periodicity_vector[index], output_more_states_unit_cell_vector[index], components);
                        command_string="cp POSCAR_OUT below_the_hull_states/index"+lexical_cast<string>(output_more_states_concentration_index[concentration_now][i]);
                        std::system(command_string.c_str());
                        
                        print_poscar_no_vacancy(output_more_states_periodicity_vector[index], output_more_states_unit_cell_vector[index], components);
                        command_string="cp POSCAR below_the_hull_states/index"+lexical_cast<string>(output_more_states_concentration_index[concentration_now][i]);
                        std::system(command_string.c_str());
                        
                    }
                    
                }
                

                
            }
            
        }
    }
    else if(global_parameters.activate_large_cell_algo_binary)
    {
        map<int,int> components;
        int x_range, y_range, z_range;
        calculate_range_from_J(J, x_range, y_range, z_range, components);
        vector<map< set<tuple<int,int,int,int,int> >, double> > cluster_type_vector;
        vector<double> formation_energy_vector;
        map< set<tuple<int,int,int,int,int> >, double> cluster_type_temp;
        vector<vector<int> > to_scan_line;
        double formation_energy_temp=0;
        vector<double> concentration_vector;
        map< set<tuple<int,int,int,int,int> >, double> J_fixed_part,J_in;
        vector<map< tuple<int,int,int,int>,int> >unit_cell_vector;
        vector<tuple<int,int,int,int,int,int> >periodicity_vector;
        
        vector<map< tuple<int,int,int,int>,int> > output_more_states_unit_cell_vector;
        vector<tuple<int,int,int,int,int,int> > output_more_states_periodicity_vector;
        vector<double> output_more_states_concentration_vector;
        vector<int> output_more_states_super_cell_size_vector;
        vector<double> output_more_states_formation_energy_vector;
        

        
        
        bool abandon=true;
        
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J.begin(); it1!=J.end(); it1++) {
            if (mu.count(it1->first)==0) {
                J_fixed_part[it1->first]=it1->second;
            }
            else{
                J_fixed_part[it1->first]=J[it1->first]-mu[it1->first];
            }
        }
        
//        formation_energy_temp=0;
//        cluster_type_temp.clear();
//        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J.begin(); it1!=J.end(); it1++) {
//            cluster_type_temp[it1->first]=0;
//        }
//        calculate_formation_energy(J, mu, cluster_type_temp, constant, mu_constant, formation_energy_temp);
//        formation_energy_vector.push_back(formation_energy_temp);
//        cluster_type_vector.push_back(cluster_type_temp);
//        concentration_vector.push_back(0);
//        map< tuple<int,int,int,int>,int> unit_cell_temp;
//        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
//            set<tuple<int,int,int,int,int> > set_temp=it1->first;
//            for (set<tuple<int,int,int,int,int> >::iterator it2=set_temp.begin(); it2!=set_temp.end(); it2++) {
//                tuple<int,int,int,int,int> tuple_temp=*it2;
//                unit_cell_temp[make_tuple(1,1,1,tuple_temp.get<3>())]=0;
//            }
//        }
//        unit_cell_vector.push_back(unit_cell_temp);
//        
//        unit_cell_temp.clear();
//        periodicity_vector.push_back(make_tuple(1,0,1,0,0,1));
//        
//        formation_energy_temp=0;
//        cluster_type_temp.clear();
//        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J.begin(); it1!=J.end(); it1++) {
//            cluster_type_temp[it1->first]=1;
//        }
//        calculate_formation_energy(J, mu, cluster_type_temp, constant, mu_constant, formation_energy_temp);
//        formation_energy_vector.push_back(formation_energy_temp);
//        cluster_type_vector.push_back(cluster_type_temp);
//        cluster_type_temp.clear();
//        
//        
//        formation_energy_temp=0;
//        concentration_vector.push_back(1);
//        unit_cell_temp.clear();
//        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
//            set<tuple<int,int,int,int,int> > set_temp=it1->first;
//            for (set<tuple<int,int,int,int,int> >::iterator it2=set_temp.begin(); it2!=set_temp.end(); it2++) {
//                tuple<int,int,int,int,int> tuple_temp=*it2;
//                unit_cell_temp[make_tuple(1,1,1,tuple_temp.get<3>())]=tuple_temp.get<4>();
//            }
//        }
//        unit_cell_vector.push_back(unit_cell_temp);
//        unit_cell_temp.clear();
//        periodicity_vector.push_back(make_tuple(1,0,1,0,0,1));
//        
//        
//        vector<int> first_line_to_check;
//        first_line_to_check.push_back(0);
//        first_line_to_check.push_back(1);
//        to_scan_line.push_back(first_line_to_check);
//        
        std::system("rm -r GS_solutions");
        std::system("mkdir GS_solutions");
//
//        std::system("mkdir GS_solutions/x0");
//        print_poscar_out(periodicity_vector[0], unit_cell_vector[0], components);
//        std::system("cp POSCAR_OUT GS_solutions/x0/");
//        print_poscar_no_vacancy(periodicity_vector[0], unit_cell_vector[0], components);
//        std::system("cp POSCAR GS_solutions/x0/");
//        
//        
//        std::system("mkdir GS_solutions/x1");
//        print_poscar_out(periodicity_vector[1], unit_cell_vector[1], components);
//        std::system("cp POSCAR_OUT GS_solutions/x1/");
//        print_poscar_no_vacancy(periodicity_vector[1], unit_cell_vector[1], components);
//        std::system("cp POSCAR GS_solutions/x1/");
        
        
        vector<double> mu_vector_numeric;
        
        global_parameters.large_cell_binary_mu_final_todo=global_parameters.large_cell_binary_mu_final;
        global_parameters.large_cell_binary_mu_init_todo=global_parameters.large_cell_binary_mu_init;
        global_parameters.large_cell_binary_mu_step_todo=global_parameters.large_cell_binary_mu_step;
        
        for (int N=0 ; N< (global_parameters.large_cell_binary_mu_final_todo-global_parameters.large_cell_binary_mu_init_todo)/global_parameters.large_cell_binary_mu_step_todo+0.1; N++) {
            double mu_now=global_parameters.large_cell_binary_mu_init_todo+global_parameters.large_cell_binary_mu_step_todo*N;
            mu_vector_numeric.push_back(mu_now);
        }

        bool should_i_stop=false;
        
        while (!should_i_stop) {
            J_in.clear();
//            vector<int> the_line_to_check=to_scan_line[0];
//            cout<<"\nthe line to chekc is";
//            printvector(the_line_to_check);
//            cout<<"  corresponding to concentration:"<<concentration_vector[the_line_to_check[0]]<<" "<<concentration_vector[the_line_to_check[1]]<<endl;
//            map< set<tuple<int,int,int,int,int> >, double>correlation_difference,correlation_1,correlation_2;
//            correlation_1=cluster_type_vector[the_line_to_check[0]];
//            correlation_2=cluster_type_vector[the_line_to_check[1]];
//            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=correlation_1.begin(); it1!=correlation_1.end(); it1++) {
//                correlation_difference[it1->first]=correlation_1[it1->first]-correlation_2[it1->first];
//            }
//            
//            double value_temp=0,value_temp_2=0,value_temp_3=0;
//            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_fixed_part.begin(); it1!=J_fixed_part.end(); it1++) {
//                value_temp-=correlation_difference[it1->first]*it1->second;
//            }
//            
//            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
//                value_temp_2+=correlation_difference[it1->first];
//            }
//            
//            value_temp_3=value_temp/value_temp_2;
//            double chemical_potential_value=value_temp_3;
//            
//            //            cout<<"\nvalue temp 3 here is"<<value_temp_3<<endl<<"\nvalue_temp_2 is"<<value_temp_2<<endl<<"value_temp_3 is"<<value_temp_3<<endl;
            
            
//            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
//                mu[it1->first]=value_temp_3;
//            }
//            
//            cout<<"\ntell me what is mu"<<endl;
//            printmapfromsets(mu);
//            
//            
//            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_fixed_part.begin(); it1!=J_fixed_part.end(); it1++) {
//                J_in[it1->first]=it1->second;
//                if (mu.count(it1->first)==1) {
//                    J_in[it1->first]+=mu[it1->first];
//                }
//            }
            J_in=J ;//this is only used for getting x_range y_range in the next step...;
            
            
            run_solver(max_sites, J_in,
                       lowerboundclustertype, upperboundclustertype, cellrepresentation,
                       lower_bound, upper_bound, exact_lowerbound, unitcell, periodicity,
                       J_for_proof,
                       id,mu,mu_constant,work_with_mu,formation_energy_UB,formation_energy_LB,new_cluster_algorithm,use_new_pair_terms,use_new_triplet_terms,map_periodicity_to_spin,use_level_method,use_weighted_dual_average,global_parameters,
                       prec, num_loops,
                       basic_exact_mode, pseudo_mode, pseudo_mode_with_proof,
                       verbose, very_verbose, obscenely_verbose,input_PRIM_output_PRIMOUT_mode,limit_dimension,constant,J_fixed_part);
            
            
//            
//            
//            if (output_more_states||output_states_below_hull) {
//                for (map<tuple<int,int,int,int,int,int>, map<tuple<int,int,int,int>,int> >::iterator it1=map_periodicity_to_spin.begin(); it1!=map_periodicity_to_spin.end(); it1++) {
//                    tuple<int,int,int,int,int,int> periodicity_now=it1->first;
//                    map<tuple<int,int,int,int>,int>  spin_now=it1->second;
//                    
//                    double formation_energy_out;
//                    double concentration_out;
//                    
//                    concentration_out=global_parameters.formation_periodic_to_concentration_transfer[periodicity_now];
//                    formation_energy_out=global_parameters.formation_periodic_to_energy_transfer[periodicity_now];
//                    //                    get_concentration_formation_energy_from_spin_and_J_in(periodicity_now, mu, J_fixed_part, spin_now, constant, mu_constant, formation_energy_out, concentration_out);
//                    
//                    bool to_add_this_structure=true;
//                    
//                    for (int ct1=0;ct1<output_more_states_concentration_vector.size();ct1++)
//                    {
//                        double concentration_to_compare=output_more_states_concentration_vector[ct1];
//                        double concentration_closeness=0;
//                        concentration_closeness=fabs(concentration_out-concentration_to_compare);
//                        if (concentration_closeness<= 1e-6)
//                        {
//                            double formation_energy_to_compare=output_more_states_formation_energy_vector[ct1];
//                            int super_cell_size_to_compare=output_more_states_super_cell_size_vector[ct1];
//                            
//                            
//                            //check to replace structure edits 109
//                            if (formation_energy_out-formation_energy_to_compare>=-1e-7 && formation_energy_out-formation_energy_to_compare<= 1e-7 && periodicity_now.get<0>()*periodicity_now.get<2>()*periodicity_now.get<5>()< super_cell_size_to_compare)
//                            {
//                                output_more_states_concentration_vector[ct1]=(concentration_out);
//                                output_more_states_periodicity_vector[ct1]=(periodicity_now);
//                                output_more_states_unit_cell_vector[ct1]=(spin_now);
//                                output_more_states_super_cell_size_vector[ct1]=(periodicity_now.get<0>()*periodicity_now.get<2>()*periodicity_now.get<5>());
//                                output_more_states_formation_energy_vector[ct1]=(formation_energy_out);
//                            }
//                            
//                            
//                            if (formation_energy_out-formation_energy_to_compare>=-1e-5 && formation_energy_out-formation_energy_to_compare<= output_states_how_sparse)
//                            {
//                                
//                                to_add_this_structure=false;
//                                break;
//                                //                                if (periodicity_now.get<0>()*periodicity_now.get<2>()*periodicity_now.get<5>()>= super_cell_size_to_compare)
//                                //                                {
//                                //                                to_add_this_structure=false;
//                                //                                break;
//                                //                                }
//                            }
//                            
//                        }
//                        
//                    }
//                    
//                    //                    cout<<"\n formation energy is "<<formation_energy_out<<" concentration out is "<<concentration_out;
//                    if (to_add_this_structure)
//                    {
//                        
//                        output_more_states_concentration_vector.push_back(concentration_out);
//                        output_more_states_periodicity_vector.push_back(periodicity_now);
//                        output_more_states_unit_cell_vector.push_back(spin_now);
//                        output_more_states_super_cell_size_vector.push_back(periodicity_now.get<0>()*periodicity_now.get<2>()*periodicity_now.get<5>());
//                        output_more_states_formation_energy_vector.push_back(formation_energy_out);
//                    }
//                    
//                    
//                }
//            }
//            
            
            
            
            
//            double concentration_temp=0;
//            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
//                concentration_temp+=upperboundclustertype[it1->first];
//            }
//            concentration_temp=concentration_temp/mu.size();
//            
//            cout<<"let me check what is formation energy_UB: "<<formation_energy_UB;
//            cout<<".  concentration temp is"<<concentration_temp;
//            bool inside_indicator=false;
//            for (int i=0; i<concentration_vector.size();i++)
//            {
//                if (fabs(concentration_vector[i]-concentration_temp)<1e-5) {
//                    inside_indicator=true;
//                }
//            }
//            
//            if (inside_indicator==false) {
//                concentration_vector.push_back(concentration_temp);
//                cluster_type_vector.push_back(upperboundclustertype);
//                formation_energy_vector.push_back(formation_energy_UB);
//                int current_index=concentration_vector.size()-1;
//                double lower_concentration=-1,upper_concentration=2;
//                int lower_index=-1000,upper_index=1000;
//                for (int i=0;i<concentration_vector.size()-1;i++) {
//                    if (concentration_vector[i]>lower_concentration&&concentration_vector[i]<concentration_temp) {
//                        lower_concentration=concentration_vector[i];
//                        lower_index=i;
//                    }
//                    if (concentration_vector[i]<upper_concentration&&concentration_vector[i]>concentration_temp) {
//                        upper_concentration=concentration_vector[i];
//                        upper_index=i;
//                    }
//                }
//                vector<int> temp_vector;
//                temp_vector.clear();
//                temp_vector.push_back(lower_index);
//                temp_vector.push_back(current_index);
//                to_scan_line.push_back(temp_vector);
//                temp_vector.clear();
//                temp_vector.push_back(current_index);
//                temp_vector.push_back(upper_index);
//                to_scan_line.push_back(temp_vector);
//                
//                unit_cell_vector.push_back(unitcell);
//                periodicity_vector.push_back(periodicity);
//                
//                std::ostringstream ss;
//                ss << std::fixed << std::setprecision(6);
//                ss << concentration_vector[current_index];
//                
//                string strings_temp="mkdir GS_solutions/x"+ss.str();
//                const char* string_command=strings_temp.c_str();
//                std::system(string_command);
//                print_poscar_out(periodicity_vector[current_index], unit_cell_vector[current_index], components);
//                strings_temp="cp POSCAR_OUT GS_solutions/x"+ss.str()+"/";
//                string_command=strings_temp.c_str();
//                std::system(string_command);
//                
//                print_poscar_no_vacancy(periodicity_vector[current_index], unit_cell_vector[current_index], components);
//                strings_temp="cp POSCAR GS_solutions/x"+ss.str()+"/";
//                string_command=strings_temp.c_str();
//                std::system(string_command);
//                
//            }
//            to_scan_line.erase(to_scan_line.begin());
            
            should_i_stop=true;
        }
        
        
        global_parameters.formation_mu_to_concentration_transfer;
        global_parameters.formation_mu_to_energy_transfer;
        global_parameters.map_mu_to_spin_transfer;
//        cout<<"\n global_parameters.formation_mu_to_concentration_transfer is"<<endl;
//        printmap(global_parameters.formation_mu_to_concentration_transfer);
//        cout<<"\n global_parameters.formation_mu_to_energy_transfer is"<<endl;
//        printmap(global_parameters.formation_mu_to_energy_transfer);
//        cout<<"\n global_parameters.map_mu_to_spin_transfer is"<<endl;
//        printmapofmap(global_parameters.map_mu_to_spin_transfer);
        tuple<int,int,int,int,int,int> periodicity_always=make_tuple(global_parameters.large_cell_binary_cubic_periodicity,0,global_parameters.large_cell_binary_cubic_periodicity,0,0,global_parameters.large_cell_binary_cubic_periodicity);
        map<double, string > mu_to_string;
        
        for (map<double, double >::iterator it1=global_parameters.formation_mu_to_energy_transfer.begin(); it1!=global_parameters.formation_mu_to_energy_transfer.end(); it1++) {
            std::ostringstream ss_temp;
            ss_temp << std::fixed << std::setprecision(6);
            ss_temp << it1->first ;
            mu_to_string[it1->first]=ss_temp.str();
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6);
            ss << global_parameters.formation_mu_to_concentration_transfer[it1->first];
            string strings_temp="mkdir GS_solutions/x"+ss.str()+"mu"+mu_to_string[it1->first];
            const char* string_command=strings_temp.c_str();
            std::system(string_command);
            print_poscar_out(periodicity_always, global_parameters.map_mu_to_spin_transfer[it1->first], components);
            strings_temp="cp POSCAR_OUT GS_solutions/x"+ss.str()+"mu"+mu_to_string[it1->first]+"/";
            string_command=strings_temp.c_str();
            std::system(string_command);

            print_poscar_no_vacancy(periodicity_always, global_parameters.map_mu_to_spin_transfer[it1->first], components);
            strings_temp="cp POSCAR GS_solutions/x"+ss.str()+"mu"+mu_to_string[it1->first]+"/";
            string_command=strings_temp.c_str();
            std::system(string_command);
            
        }
        
//        map<double,double> concentration_formation_energy;
//        map<double,int> concentration_index;
//        map<double,int> concentration_super_cell_size;
//        map<double,tuple<int,int,int,int,int,int> > concentration_periodicity;
//        
//        cout<<"\n after everything showing the results:";
//        for (int i=0; i<concentration_vector.size(); i++) {
//            concentration_formation_energy[concentration_vector[i]]=formation_energy_vector[i];
//            concentration_index[concentration_vector[i]]=i;
//            tuple<int,int,int,int,int,int> periodicity_temp;
//            periodicity_temp=periodicity_vector[i];
//            concentration_super_cell_size[concentration_vector[i]]=periodicity_temp.get<0>()*periodicity_temp.get<2>()*periodicity_temp.get<5>();
//            concentration_periodicity[concentration_vector[i]]=periodicity_temp;
//        }
        
        stringstream to_print_file;
        to_print_file<<setw(10)<<"x"<<setw(20)<<"formation_energy"<<setw(20)<<"super_cell_size"<<setw(20)<<"mu"<<endl;
        for (map<double, double >::iterator it1=global_parameters.formation_mu_to_energy_transfer.begin(); it1!=global_parameters.formation_mu_to_energy_transfer.end(); it1++) {
            
            to_print_file<<setw(10)<< std::fixed << std::setprecision(6)<<global_parameters.formation_mu_to_concentration_transfer[it1->first]<<setw(20)<< std::fixed << std::setprecision(8)<<it1->second<<setw(20)<< std::fixed <<myPow(global_parameters.large_cell_binary_cubic_periodicity, 3) <<setw(20)<< std::fixed<<mu_to_string[it1->first]<< endl;
            
        }
        cout<<"\n"<<to_print_file.str();
        
        const char *path="./GS_solutions/hull.txt";
        std::ofstream file(path);
        file<<to_print_file.str();
        
        
//        stringstream to_print_file_debug;
//        to_print_file_debug<<setw(10)<<"x"<<setw(20)<<"formation_energy"<<setw(20)<<"super_cell_size"<<setw(20)<<"periodicity"<<endl;
//        for (map<double, double>::iterator it1=concentration_formation_energy.begin(); it1!=concentration_formation_energy.end(); it1++) {
//            
//            to_print_file_debug<<setw(10)<< std::fixed << std::setprecision(6)<<it1->first<<setw(20)<< std::fixed << std::setprecision(8)<<it1->second<<setw(20)<< std::fixed <<concentration_super_cell_size[it1->first]<<setw(20)<< std::fixed <<concentration_periodicity[it1->first]<<endl;
//            
//        }
//        cout<<"\n"<<to_print_file_debug.str();
//        
//        const char *path1="./GS_solutions/hull_debug.txt";
//        std::ofstream file1(path1);
//        file1<<to_print_file_debug.str();
        
    }
    else if (global_parameters.do_monte_carlo)
    {
        spin_periodic_struct monte_spin_struct;
        double monte_energy;
        simulated_annealing(J,global_parameters,monte_energy,monte_spin_struct,constant);
        
        
    }
    else if (global_parameters.ternary_alg){
        cout<<"\n in ternary alg";
//        cout<<"\n constant is "<<constant;

        mu.clear();
        mu.insert(global_parameters.mu_tern_first_group.begin(),global_parameters.mu_tern_first_group.end());
        mu.insert(global_parameters.mu_tern_second_group.begin(),global_parameters.mu_tern_second_group.end());

        
        map<set< tuple<int,int,int,int,int> >,double > J_tern_in;
        convert_J_input_tern_unaltered_to_J_tern_in(global_parameters.J_input_tern_unaltered, J_tern_in);

        
//        merging a few maps
        J_tern_in.insert(global_parameters.mu_tern_first_group.begin(),global_parameters.mu_tern_first_group.end());
        J_tern_in.insert(global_parameters.mu_tern_second_group.begin(),global_parameters.mu_tern_second_group.end());
        
//        cout<<"\n debug again check J_tern_in";
//        printmapfromsets(J_tern_in);
        
        
        map<int,int> components;
        int x_range, y_range, z_range;
        calculate_range_from_J(J_tern_in, x_range, y_range, z_range, components);
        
        
        vector<map< set<tuple<int,int,int,int,int> >, double> > cluster_type_vector;
        vector<double> formation_energy_vector;
        vector<double> mu1_vector, mu2_vector;
        map< set<tuple<int,int,int,int,int> >, double> cluster_type_temp;
        vector<vector<int> > to_scan_facet;
        set<set<int> > already_scanned_facet;
        double formation_energy_temp=0;
        
        vector<tuple<double,double> > concentration_vector;
        map< set<tuple<int,int,int,int,int> >, double> J_fixed_part,J_in;
        vector<map< tuple<int,int,int,int>,int> >unit_cell_vector;
        vector<tuple<int,int,int,int,int,int> >periodicity_vector;

        
        
        vector<map< tuple<int,int,int,int>,int> > output_more_states_unit_cell_vector;
        vector<tuple<int,int,int,int,int,int> > output_more_states_periodicity_vector;
        vector<tuple<double,double> > output_more_states_concentration_vector;
        vector<int> output_more_states_super_cell_size_vector;
        vector<double> output_more_states_formation_energy_vector;
        
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_tern_in.begin(); it1!=J_tern_in.end(); it1++) {
            if (mu.count(it1->first)==0) {
                J_fixed_part[it1->first]=it1->second;
            }
            else{
                J_fixed_part[it1->first]=J_tern_in[it1->first]-mu[it1->first];
            }
        }
        
        
        //the 0 point
        formation_energy_temp=0;
        cluster_type_temp.clear();
        
        
        map< tuple<int,int,int,int>,int> unit_cell_temp;
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
            set<tuple<int,int,int,int,int> > set_temp=it1->first;
            for (set<tuple<int,int,int,int,int> >::iterator it2=set_temp.begin(); it2!=set_temp.end(); it2++) {
                tuple<int,int,int,int,int> tuple_temp=*it2;
                unit_cell_temp[make_tuple(1,1,1,tuple_temp.get<3>())]=0;
            }
        }
        
        unit_cell_temp=extend_periodic_spin(make_tuple(1,0,1,0,0,1), unit_cell_temp, J_tern_in);
        
        
        
        tuple<double,double> concentration_temp;
        
        get_concentration_formation_energy_from_spin_and_J_in_ternary(make_tuple(1,0,1,0,0,1), global_parameters.mu_tern_first_group, global_parameters.mu_tern_second_group, J_fixed_part, unit_cell_temp, constant, mu_constant, formation_energy_temp, concentration_temp,cluster_type_temp);
        
        
        formation_energy_vector.push_back(formation_energy_temp);
        cluster_type_vector.push_back(cluster_type_temp);
        concentration_vector.push_back(concentration_temp);
        
        
        
        unit_cell_vector.push_back(unit_cell_temp);
        unit_cell_temp.clear();
        periodicity_vector.push_back(make_tuple(1,0,1,0,0,1));
        
        mu1_vector.push_back(10000);
        mu2_vector.push_back(10000);
        
        
        
        //get the first mu point
        formation_energy_temp=0;
        cluster_type_temp.clear();
        
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=global_parameters.mu_tern_first_group.begin(); it1!=global_parameters.mu_tern_first_group.end(); it1++) {
            set<tuple<int,int,int,int,int> > set_temp=it1->first;
            for (set<tuple<int,int,int,int,int> >::iterator it2=set_temp.begin(); it2!=set_temp.end(); it2++) {
                tuple<int,int,int,int,int> tuple_temp=*it2;
                unit_cell_temp[make_tuple(1,1,1,tuple_temp.get<3>())]=tuple_temp.get<4>();
            }
        }
        
        unit_cell_temp=extend_periodic_spin(make_tuple(1,0,1,0,0,1), unit_cell_temp, J_tern_in);
        

        
        get_concentration_formation_energy_from_spin_and_J_in_ternary(make_tuple(1,0,1,0,0,1), global_parameters.mu_tern_first_group, global_parameters.mu_tern_second_group, J_fixed_part, unit_cell_temp, constant, mu_constant, formation_energy_temp, concentration_temp,cluster_type_temp);
        
        
        formation_energy_vector.push_back(formation_energy_temp);
        cluster_type_vector.push_back(cluster_type_temp);
        concentration_vector.push_back(concentration_temp);
        
        
        
        
        unit_cell_vector.push_back(unit_cell_temp);
        unit_cell_temp.clear();
        periodicity_vector.push_back(make_tuple(1,0,1,0,0,1));
        
        
        mu1_vector.push_back(-10000);
        mu2_vector.push_back(0);
        
        //get second mu point
        
        
        formation_energy_temp=0;
        cluster_type_temp.clear();
        
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=global_parameters.mu_tern_second_group.begin(); it1!=global_parameters.mu_tern_second_group.end(); it1++) {
            set<tuple<int,int,int,int,int> > set_temp=it1->first;
            for (set<tuple<int,int,int,int,int> >::iterator it2=set_temp.begin(); it2!=set_temp.end(); it2++) {
                tuple<int,int,int,int,int> tuple_temp=*it2;
                unit_cell_temp[make_tuple(1,1,1,tuple_temp.get<3>())]=tuple_temp.get<4>();
            }
        }
        
        unit_cell_temp=extend_periodic_spin(make_tuple(1,0,1,0,0,1), unit_cell_temp, J_tern_in);
        
        
        get_concentration_formation_energy_from_spin_and_J_in_ternary(make_tuple(1,0,1,0,0,1), global_parameters.mu_tern_first_group, global_parameters.mu_tern_second_group, J_fixed_part, unit_cell_temp, constant, mu_constant, formation_energy_temp, concentration_temp,cluster_type_temp);
        
        formation_energy_vector.push_back(formation_energy_temp);
        cluster_type_vector.push_back(cluster_type_temp);
        concentration_vector.push_back(concentration_temp);
        
        
//        cout<<"\nafter formation energy is "<<formation_energy_temp;
//        cout<<"\n cluster_type_temp is";
//        printmapfromsets(cluster_type_temp);
//        cout<<"\n unit cell temp is ";
//        printmap(unit_cell_temp);
//        cout<<"\n concentration_temp is: "<<concentration_temp;
        
        
        
        unit_cell_vector.push_back(unit_cell_temp);
        unit_cell_temp.clear();
        periodicity_vector.push_back(make_tuple(1,0,1,0,0,1));
        
        mu1_vector.push_back(0);
        mu2_vector.push_back(-10000);
        
        vector<int> first_facet_to_check;
        first_facet_to_check.push_back(0);
        first_facet_to_check.push_back(1);
        first_facet_to_check.push_back(2);
        to_scan_facet.push_back(first_facet_to_check);
    
        std::system("rm -r GS_solutions");
        std::system("mkdir GS_solutions");
        
        bool tern_alg=true;
        std::system("mkdir GS_solutions/x0.000000y0.000000");
        print_poscar_out(periodicity_vector[0], unit_cell_vector[0], components,tern_alg);
        std::system("cp POSCAR_OUT GS_solutions/x0.000000y0.000000");
        print_poscar_no_vacancy(periodicity_vector[0], unit_cell_vector[0], components,tern_alg);
        std::system("cp POSCAR GS_solutions/x0.000000y0.000000/");
        
        
        std::system("mkdir GS_solutions/x1.000000y0.000000");
        print_poscar_out(periodicity_vector[1], unit_cell_vector[1], components,tern_alg);
        std::system("cp POSCAR_OUT GS_solutions/x1.000000y0.000000/");
        print_poscar_no_vacancy(periodicity_vector[1], unit_cell_vector[1], components,tern_alg);
        std::system("cp POSCAR GS_solutions/x1.000000y0.000000/");
        
        
        std::system("mkdir GS_solutions/x0.000000y1.000000");
        print_poscar_out(periodicity_vector[2], unit_cell_vector[2], components,tern_alg);
        std::system("cp POSCAR_OUT GS_solutions/x0.000000y1.000000/");
        print_poscar_no_vacancy(periodicity_vector[2], unit_cell_vector[2], components,tern_alg);
        std::system("cp POSCAR GS_solutions/x0.000000y1.000000/");
        
        
        while (!to_scan_facet.empty()) {
        
            double first_chemical_potential=0;
            double second_chemical_potential=0;
            
            J_in.clear();
            vector<int> the_facet_to_check=to_scan_facet[0];
            cout<<"\nthe facet to chekc is";
            printvector(the_facet_to_check);
            cout<<"  corresponding to concentration:"<<concentration_vector[the_facet_to_check[0]]<<" "<<concentration_vector[the_facet_to_check[1]]<<" "<<concentration_vector[the_facet_to_check[2]]<<endl;
            map< set<tuple<int,int,int,int,int> >, double>correlation_difference_1_2,correlation_difference_1_3,correlation_1,correlation_2,correlation_3;
            correlation_1=cluster_type_vector[the_facet_to_check[0]];
            correlation_2=cluster_type_vector[the_facet_to_check[1]];
            correlation_3=cluster_type_vector[the_facet_to_check[2]];
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=correlation_1.begin(); it1!=correlation_1.end(); it1++) {
                correlation_difference_1_2[it1->first]=correlation_1[it1->first]-correlation_2[it1->first];
                correlation_difference_1_3[it1->first]=correlation_1[it1->first]-correlation_3[it1->first];
            }
            
            double neg_corr_1_2_dot_J_fix=0,neg_corr_1_3_dot_J_fix=0,corr_1_2_dot_MU_1=0,corr_1_2_dot_MU_2=0,corr_1_3_dot_MU_1=0,corr_1_3_dot_MU_2=0;
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_fixed_part.begin(); it1!=J_fixed_part.end(); it1++) {
//                value_temp-=correlation_difference[it1->first]*it1->second;
                neg_corr_1_2_dot_J_fix-=correlation_difference_1_2[it1->first]*it1->second;
                neg_corr_1_3_dot_J_fix-=correlation_difference_1_3[it1->first]*it1->second;
            }
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=global_parameters.mu_tern_first_group.begin(); it1!=global_parameters.mu_tern_first_group.end(); it1++) {
                
                corr_1_2_dot_MU_1+=correlation_difference_1_2[it1->first];
                corr_1_3_dot_MU_1+=correlation_difference_1_3[it1->first];
            }
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=global_parameters.mu_tern_second_group.begin(); it1!=global_parameters.mu_tern_second_group.end(); it1++) {
                
                corr_1_2_dot_MU_2+=correlation_difference_1_2[it1->first];
                corr_1_3_dot_MU_2+=correlation_difference_1_3[it1->first];
            }
            
            MatrixXd M1(2,2);
            M1(0,0)=corr_1_2_dot_MU_1;
            M1(0,1)=corr_1_2_dot_MU_2;
            M1(1,0)=corr_1_3_dot_MU_1;
            M1(1,1)=corr_1_3_dot_MU_2;
            VectorXd V1(2);
            V1(0)=neg_corr_1_2_dot_J_fix;
            V1(1)=neg_corr_1_3_dot_J_fix;
            VectorXd x(2);
            x=M1.inverse()*V1;
            
            
//            cout<<"\nlet me check M1.inverse() :\n"<<M1.inverse();
//            cout<<"\nlet me check x vector:\n"<<x;
//            cout<<"\ncheck if x[0]>1e10:"<<(x[0]>1e10);
//            cout<<"\n M1.determinant():"<<M1.determinant();
            
            if (fabs(M1.determinant())<1e-10) {
                cout<<"\nthis facet is colinnear, skip";
                to_scan_facet.erase(to_scan_facet.begin());
                continue;
            }
            
            
            first_chemical_potential=x(0);
            second_chemical_potential=x(1);
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=global_parameters.mu_tern_first_group.begin(); it1!=global_parameters.mu_tern_first_group.end(); it1++) {
                mu[it1->first]=first_chemical_potential;

            }
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=global_parameters.mu_tern_second_group.begin(); it1!=global_parameters.mu_tern_second_group.end(); it1++) {
                mu[it1->first]=second_chemical_potential;
            }
            
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_fixed_part.begin(); it1!=J_fixed_part.end(); it1++) {
                J_in[it1->first]=it1->second;
                if (mu.count(it1->first)==1) {
                    J_in[it1->first]+=mu[it1->first];
                }
            }
            

            run_solver(max_sites, J_in,
                       lowerboundclustertype, upperboundclustertype, cellrepresentation,
                       lower_bound, upper_bound, exact_lowerbound, unitcell, periodicity,
                       J_for_proof,
                       id,mu,mu_constant,work_with_mu,formation_energy_UB,formation_energy_LB,new_cluster_algorithm,use_new_pair_terms,use_new_triplet_terms,map_periodicity_to_spin,use_level_method,use_weighted_dual_average,global_parameters,
                       prec, num_loops,
                       basic_exact_mode, pseudo_mode, pseudo_mode_with_proof,
                       verbose, very_verbose, obscenely_verbose,input_PRIM_output_PRIMOUT_mode,limit_dimension,constant);
            
            
            
            if (global_parameters.ternary_output_states_above_hull)
            {
                
                for (map<tuple<int,int,int,int,int,int>, map<tuple<int,int,int,int>,int> >::iterator it1=map_periodicity_to_spin.begin(); it1!=map_periodicity_to_spin.end(); it1++) {
                    tuple<int,int,int,int,int,int> periodicity_now=it1->first;
                    map<tuple<int,int,int,int>,int>  spin_now=it1->second;
                    map< set<tuple<int,int,int,int,int> >, double> cluster_type_temp_1;
                    double formation_energy_out;
                    tuple<double,double> concentration_out;
                    
                    get_concentration_formation_energy_from_spin_and_J_in_ternary(periodicity_now, global_parameters.mu_tern_first_group, global_parameters.mu_tern_second_group, J_fixed_part, spin_now, constant, mu_constant, formation_energy_out, concentration_out,cluster_type_temp_1);
                    
                    bool to_add_this_structure=true;
                    
                    for (int ct1=0;ct1<output_more_states_concentration_vector.size();ct1++)
                    {
                        tuple<double,double> concentration_to_compare=output_more_states_concentration_vector[ct1];
                        double concentration_closeness=0;
                        concentration_closeness=fabs(concentration_out.get<0>()-concentration_to_compare.get<0>())+fabs(concentration_out.get<1>()-concentration_to_compare.get<1>());
                        if (concentration_closeness<= 1e-6)
                        {
                            double formation_energy_to_compare=output_more_states_formation_energy_vector[ct1];
                            int super_cell_size_to_compare=output_more_states_super_cell_size_vector[ct1];
                            if (formation_energy_out-formation_energy_to_compare>=-1e-5 && formation_energy_out-formation_energy_to_compare<= global_parameters.ternary_output_states_above_hull_gap)
                            {
                                if (periodicity_now.get<0>()*periodicity_now.get<2>()*periodicity_now.get<5>()>= super_cell_size_to_compare)
                                {
                                    to_add_this_structure=false;
                                }
                            }
                            
                        }
                        
                        
                        
                    }
                    
                    if (to_add_this_structure)
                    {
                    
                        output_more_states_concentration_vector.push_back(concentration_out);
                        output_more_states_periodicity_vector.push_back(periodicity_now);
                        output_more_states_unit_cell_vector.push_back(spin_now);
                        output_more_states_super_cell_size_vector.push_back(periodicity_now.get<0>()*periodicity_now.get<2>()*periodicity_now.get<5>());
                        output_more_states_formation_energy_vector.push_back(formation_energy_out);
                    
                    }
                    
                }
                1;
//                details to be added later
            }
            
            
//            cout<<"\nlet me check what is formation energy_UB: "<<formation_energy_UB;
            
            tuple<double,double> concentration_temp;
            double formation_energy_temp;
            map<set<tuple<int, int, int, int, int> >, double> clustertype_out_temp;
            get_concentration_formation_energy_from_spin_and_J_in_ternary(periodicity, global_parameters.mu_tern_first_group, global_parameters.mu_tern_second_group, J_fixed_part, unitcell, constant, mu_constant, formation_energy_temp, concentration_temp, clustertype_out_temp);
            
            cout<<"\n the other way to compute formation energy is "<<formation_energy_temp;
            cout<<"\n the concentration is "<<concentration_temp;
            cout<<"\n the chemical potential used is "<<make_tuple(first_chemical_potential,second_chemical_potential);
            
            
            bool inside_indicator=false;

            
            
            
            for (int i=0; i<concentration_vector.size();i++)
            {
                
                if (fabs(concentration_vector[i].get<0>()-concentration_temp.get<0>())+fabs(concentration_vector[i].get<1>()-concentration_temp.get<1>())<1e-7) {
                    inside_indicator=true;
                }
            }
            
            
            if (inside_indicator==false) {
                concentration_vector.push_back(concentration_temp);
                cluster_type_vector.push_back(clustertype_out_temp);
                formation_energy_vector.push_back(formation_energy_temp);
                unit_cell_vector.push_back(unitcell);
                periodicity_vector.push_back(periodicity);
                
                
                //ok, but what is happening now?
                long long current_index=concentration_vector.size()-1;
                
                vector<int> temp_vector;
                temp_vector.clear();
                temp_vector.push_back(the_facet_to_check[0]);
                temp_vector.push_back(the_facet_to_check[1]);
                temp_vector.push_back(current_index);
                to_scan_facet.push_back(temp_vector);
                
                temp_vector.clear();
                temp_vector.push_back(the_facet_to_check[1]);
                temp_vector.push_back(the_facet_to_check[2]);
                temp_vector.push_back(current_index);
                to_scan_facet.push_back(temp_vector);
                
                temp_vector.clear();
                temp_vector.push_back(the_facet_to_check[0]);
                temp_vector.push_back(the_facet_to_check[2]);
                temp_vector.push_back(current_index);
                to_scan_facet.push_back(temp_vector);
                
                std::ostringstream ss1,ss2;
                ss1 << std::fixed << std::setprecision(6);
                ss1 << concentration_vector[current_index].get<0>();
                ss2 << std::fixed << std::setprecision(6);
                ss2 << concentration_vector[current_index].get<1>();
                
                string strings_temp="mkdir GS_solutions/x"+ss1.str()+"y"+ss2.str();
                const char* string_command=strings_temp.c_str();
                std::system(string_command);
                
//                print_poscar_out(periodicity_vector[current_index], unit_cell_vector[current_index], components);
                bool tern_alg=true;
                print_poscar_out(periodicity_vector[current_index], unit_cell_vector[current_index], components,tern_alg);
                strings_temp="cp POSCAR_OUT GS_solutions/x"+ss1.str()+"y"+ss2.str()+"/";
                string_command=strings_temp.c_str();
                std::system(string_command);
                
                print_poscar_no_vacancy(periodicity_vector[current_index], unit_cell_vector[current_index], components,tern_alg);
                strings_temp="cp POSCAR GS_solutions/x"+ss1.str()+"y"+ss2.str()+"/";
                string_command=strings_temp.c_str();
                std::system(string_command);
                
                mu1_vector.push_back(first_chemical_potential);
                mu2_vector.push_back(second_chemical_potential);
                
                
            }
            
            
            
            
            to_scan_facet.erase(to_scan_facet.begin());
            set <int> scanned_facet_this_time(the_facet_to_check.begin(),the_facet_to_check.end());
            already_scanned_facet.insert(scanned_facet_this_time);
            
            if (inside_indicator==false) {
                vector<vector<double> > corre_dot_ex_corre_dot_ey_f_fix_dot_corre_vector;
                for(int li=0;li<unit_cell_vector.size();li++){
                    
                    map< set<tuple<int,int,int,int,int> >, double> cluster_type_for_li=cluster_type_vector[li];
                    
                    double corre_dot_ex=0,corre_dot_ey=0,f_fix_dot_corre=0;
                    
                    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=global_parameters.mu_tern_first_group.begin(); it1!=global_parameters.mu_tern_first_group.end(); it1++) {
                        corre_dot_ex+=cluster_type_for_li[it1->first];
                    }
                    
                    
                    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=global_parameters.mu_tern_second_group.begin(); it1!=global_parameters.mu_tern_second_group.end(); it1++) {
                        corre_dot_ey+=cluster_type_for_li[it1->first];
                    }
                    
                    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_fixed_part.begin(); it1!=J_fixed_part.end(); it1++) {
                        //                value_temp-=correlation_difference[it1->first]*it1->second;
                        f_fix_dot_corre+=cluster_type_for_li[it1->first]*it1->second;
                    }
                    
                    vector <double> v1_temp;
                    v1_temp.push_back(corre_dot_ex);
                    v1_temp.push_back(corre_dot_ey);
                    v1_temp.push_back(f_fix_dot_corre);
                    corre_dot_ex_corre_dot_ey_f_fix_dot_corre_vector.push_back(v1_temp);
                    
                }
                
                ofstream qhull_input;
                string qhull_file_name = id + "_" + "qhull_input.in";
                qhull_input.open(qhull_file_name.c_str());
                
                stringstream qhull_stringstream;
                qhull_stringstream<<"3\n"<<unit_cell_vector.size()<<"\n";
                for (int li=0; li<corre_dot_ex_corre_dot_ey_f_fix_dot_corre_vector.size(); li++)
                {
                    qhull_stringstream<<corre_dot_ex_corre_dot_ey_f_fix_dot_corre_vector[li][0]<<" "<<corre_dot_ex_corre_dot_ey_f_fix_dot_corre_vector[li][1]<<" "<<corre_dot_ex_corre_dot_ey_f_fix_dot_corre_vector[li][2]<<"\n";
                }
                qhull_input<<qhull_stringstream.str();
                qhull_input.close();
                
                string result_now = exec("qhull Fv <" + qhull_file_name);
                cout<<"\n show me the qhull result:\n";
                cout<<result_now;
                
                vector<string> result_pieces;
                split(result_now, '\n', result_pieces);
                
                to_scan_facet.clear();
                
                for (int li=1;li<result_pieces.size();li++){
                    string temp_string=result_pieces[li];
                    vector<string> temp_segment;
                    split(temp_string, ' ', temp_segment);
                    string first_segment = *temp_segment.begin();
                    assert(lexical_cast<int>(first_segment) >=3);
                    assert(temp_segment.size()>=4);
                    
                    vector<int> this_one_facet;
                    set<int> this_one_facet_set;
                    for (int l2=1; l2<temp_segment.size(); l2++) {
                        int num_here=lexical_cast<int>(temp_segment[l2]);
                        if (l2<=3) {
                            this_one_facet.push_back(num_here);
                        }
                        this_one_facet_set.insert(num_here);
                    }
                    
                    
                    
                    bool this_one_facet_is_new=true;
                    for (set<set<int> >::iterator it1=already_scanned_facet.begin(); it1!=already_scanned_facet.end(); it1++) {
                        set<int> one_scanned_facet=*it1;
                        if (includes(this_one_facet_set.begin(), this_one_facet_set.end(), one_scanned_facet.begin(), one_scanned_facet.end())) {
                            this_one_facet_is_new=false;
                            break;
                        }
                    }
                    
                    if (this_one_facet_is_new) {
                        to_scan_facet.push_back(this_one_facet);
                    }
                    
                    
                };
                
                
                
                
            }
            
            
            
            
            
            
            
        };
        
        
        map<tuple<double,double>,double> concentration_formation_energy;
        map<tuple<double,double>,int> concentration_index;
        map<tuple<double,double>,int> concentration_super_cell_size;
        
        map<tuple<double,double>,tuple<double,double> > concentration_mu1_mu2;
        
        cout<<"\n after everything showing the results:";
        for (int i=0; i<concentration_vector.size(); i++) {
            
            concentration_formation_energy[concentration_vector[i]]=formation_energy_vector[i];
            concentration_index[concentration_vector[i]]=i;
            tuple<int,int,int,int,int,int> periodicity_temp;
            periodicity_temp=periodicity_vector[i];
            concentration_super_cell_size[concentration_vector[i]]=periodicity_temp.get<0>()*periodicity_temp.get<2>()*periodicity_temp.get<5>();
            concentration_mu1_mu2[concentration_vector[i]]=make_tuple(mu1_vector[i],mu2_vector[i]);
            
        }
        
        stringstream to_print_file;
    to_print_file<<setw(10)<<"x"<<setw(10)<<"y"<<setw(20)<<"formation_energy"<<setw(20)<<"super_cell_size"<<setw(25)<<"mu_x"<<setw(25)<<"mu_y"<<endl;
        
        for (map<tuple<double,double>, double>::iterator it1=concentration_formation_energy.begin(); it1!=concentration_formation_energy.end(); it1++) {
            
            to_print_file<<setw(10)<< std::fixed << std::setprecision(6)<<it1->first.get<0>()<<setw(10)<< std::fixed << std::setprecision(6)<<it1->first.get<1>()<<setw(20)<< std::fixed << std::setprecision(8)<<it1->second<<setw(20)<< std::fixed <<concentration_super_cell_size[it1->first]<<setw(25)<< std::fixed << std::setprecision(6)<<-concentration_mu1_mu2[it1->first].get<0>()<<setw(25)<< std::fixed << std::setprecision(6)<<-concentration_mu1_mu2[it1->first].get<1>()<<endl;
            
            
        }
        cout<<"\n"<<to_print_file.str();
        
        const char *path="./GS_solutions/hull.txt";
        std::ofstream file(path);
        file<<to_print_file.str();
        
        
        if (global_parameters.ternary_output_states_above_hull)
        {
            
            map< tuple<double,double>, vector<double> > output_more_states_concentration_formation_energy;
            map<tuple<double,double>, vector<int> > output_more_states_concentration_index;
            map<tuple<double,double>,vector<int> > output_more_states_concentration_super_cell_size;
            map<tuple<double,double>,vector<tuple<int,int,int,int,int,int> > > output_more_states_concentration_periodicity;
            
            
            {
                vector<double> output_more_states_formation_energy_vector;
                vector<map< tuple<int,int,int,int>,int> > output_more_states_unit_cell_vector;
                vector<tuple<int,int,int,int,int,int> > output_more_states_periodicity_vector;
                vector<double> output_more_states_concentration_vector;
                vector<int> output_more_states_super_cell_size_vector;
            }
            
        
            map<int, vector<int> > supercell_size_to_index;
            
            for(int size=1; size<=max_sites;size++)
            {
                for (int i=0; i<output_more_states_super_cell_size_vector.size();i++)
                {
                    if(output_more_states_super_cell_size_vector[i]==size)
                    {
                       supercell_size_to_index[size].push_back(i);
                    }
                }
            }
            
            for(map<int, vector<int> >::iterator it1=supercell_size_to_index.begin();it1!=supercell_size_to_index.end();it1++)
            {
                int super_cell_size=it1->first;
                vector<int> index_list=it1->second;
                for (int index_ct=0;index_ct<index_list.size();index_ct++)
                {
                    map< tuple<int,int,int,int>,int> unit_cell_now=output_more_states_unit_cell_vector[index_list[index_ct]];
                    tuple<int,int,int,int,int,int> periodicity_now=output_more_states_periodicity_vector[index_list[index_ct]];
                    tuple<double,double> concentration_now=output_more_states_concentration_vector[index_list[index_ct]];
                    double formation_energy_now=output_more_states_formation_energy_vector[index_list[index_ct]];
                    int super_cell_size_now=output_more_states_super_cell_size_vector[index_list[index_ct]];
                    int index_now=index_list[index_ct];
                    bool add_this_into_map=true;
                    
                    tuple<double,double> concentration_to_compare_standard=concentration_now;
                    for(map< tuple<double,double>, vector<double> >::iterator it2=output_more_states_concentration_formation_energy.begin();it2!=output_more_states_concentration_formation_energy.end();it2++)
                    {
                        tuple<double,double> concentration_to_compare=it2->first;
                        vector<double> formation_energy_to_compare_vector=it2->second;
                                                                                       
                        bool concentration_matched=false;
                        double concentration_closeness=fabs(concentration_now.get<0>()-concentration_to_compare.get<0>())+fabs(concentration_now.get<1>()-concentration_to_compare.get<1>());



                        if  (concentration_closeness<1e-7)
                        {
                            concentration_matched=true;
                            concentration_to_compare_standard=concentration_to_compare;
                                                                                      
                                                                                           
                            for(int to_compare_ct=0;to_compare_ct<formation_energy_to_compare_vector.size();to_compare_ct++)
                            {
                                double formation_energy_to_compare=formation_energy_to_compare_vector[to_compare_ct];
                                int super_cell_size_to_compare=output_more_states_concentration_super_cell_size[concentration_to_compare][to_compare_ct];
                                if (formation_energy_now-formation_energy_to_compare>=-1e-6 && formation_energy_now-formation_energy_to_compare<=global_parameters.ternary_output_states_above_hull_gap && super_cell_size_to_compare<=super_cell_size_now)
                                {
                                    add_this_into_map=false;
                                }
                                
                                if(add_this_into_map==false)
                                {
                                    break;
                                }

                            }
                            
                            if(add_this_into_map==false)
                            {
                                break;
                            }

                        }
                        
                    }
                    
                    
                    double x_obj=concentration_now.get<0>();
                    double y_obj=concentration_now.get<1>();
                    double z_obj=1-x_obj-y_obj;
                    
                    if (add_this_into_map==true&&!(x_obj<global_parameters.ternary_x_min-1e-7 || x_obj>global_parameters.ternary_x_max+1e-7 || y_obj<global_parameters.ternary_y_min-1e-7|| y_obj>global_parameters.ternary_y_max+1e-7 || z_obj<global_parameters.ternary_z_min-1e-7||z_obj>global_parameters.ternary_z_max+1e-7))
                    {
                        output_more_states_concentration_formation_energy[concentration_to_compare_standard].push_back(formation_energy_now);
                        output_more_states_concentration_index[concentration_to_compare_standard].push_back(index_now);
                        output_more_states_concentration_super_cell_size[concentration_to_compare_standard].push_back(super_cell_size_now);
                        output_more_states_concentration_periodicity[concentration_to_compare_standard].push_back(periodicity_now);
                    }
                }
            }
            
            
            for(map< tuple<double,double>, vector<double> >::iterator it1=output_more_states_concentration_formation_energy.begin();it1!=output_more_states_concentration_formation_energy.end();it1++)
            {
                tuple<double,double> concentration_now=it1->first;
            
                GRBEnv env = GRBEnv();
                GRBModel m = GRBModel(env);
                
                double x_obj=concentration_now.get<0>();
                double y_obj=concentration_now.get<1>();
                double z_obj=1-x_obj-y_obj;
                

                
                vector<tuple<double,double,double> > concentration_formation_E;
                for (map<tuple<double,double>,double>::iterator it1=concentration_formation_energy.begin();it1!=concentration_formation_energy.end();it1++)
                {
                    tuple<double,double> concentration_here=it1->first;
                    double formation_energy_here=it1->second;
                    concentration_formation_E.push_back(make_tuple(concentration_here.get<0>(),concentration_here.get<1>(),formation_energy_here));
                    
                }
                
                
                vector<GRBVar> weight_vector;
                for (int i=0; i<concentration_formation_E.size(); i++) {
                    weight_vector.push_back(m.addVar(0, 1, concentration_formation_E[i].get<2>(), GRB_CONTINUOUS));
                }
                
                m.update();
                
                GRBLinExpr x_constraint=0,y_constraint=0,sum_constraint=0;
                for (int i=0; i<concentration_formation_E.size(); i++) {
                    x_constraint+=weight_vector[i]*concentration_formation_E[i].get<0>();
                    y_constraint+=weight_vector[i]*concentration_formation_E[i].get<1>();
                    sum_constraint+=weight_vector[i];
                }
                
                m.addConstr(x_constraint==x_obj);
                m.addConstr(y_constraint==y_obj);
                m.addConstr(sum_constraint==1);
                m.optimize();
                double hull_value=m.get(GRB_DoubleAttr_ObjVal);
                
                for (int ct=0;ct < output_more_states_concentration_formation_energy[concentration_now].size() ; ct++) {
                    
                    if (output_more_states_concentration_formation_energy[concentration_now][ct]-hull_value>global_parameters.ternary_output_states_above_hull_ceil) {
                    
                        output_more_states_concentration_formation_energy[concentration_now].erase(output_more_states_concentration_formation_energy[concentration_now].begin()+ct);
                        
                        output_more_states_concentration_index[concentration_now].erase(output_more_states_concentration_index[concentration_now].begin()+ct);
                        
                        output_more_states_concentration_super_cell_size[concentration_now].erase(output_more_states_concentration_super_cell_size[concentration_now].begin()+ct);
                        
                       output_more_states_concentration_periodicity[concentration_now].erase(output_more_states_concentration_periodicity[concentration_now].begin()+ct);
                        
                        ct--;
                    }

                    
                    
                }
                
            }
            
            
            stringstream output_more_states_string;
            output_more_states_string<<setw(10)<<"x"<<setw(10)<<"y"<<setw(20)<<"formation_energy"<<setw(20)<<"super_cell_size"<<setw(20)<<"corresponding_index"<<endl;
            
            for (map<tuple<double,double>,vector<double> >::iterator it1=output_more_states_concentration_formation_energy.begin(); it1!=output_more_states_concentration_formation_energy.end(); it1++) {
                
                tuple<double,double> concentration_now=it1->first;
                vector<double> energy_vector=it1->second;
                int size_to_enumerate_over=energy_vector.size();
                
                for (int i=0; i<size_to_enumerate_over; i++) {
                    
                    output_more_states_string<<setw(10)<< std::fixed << std::setprecision(6)<<concentration_now.get<0>()<<setw(10)<< std::fixed << std::setprecision(6)<<concentration_now.get<1>()<<setw(20)<< std::fixed << std::setprecision(8)<<energy_vector[i]<<setw(20)<< std::fixed <<output_more_states_concentration_super_cell_size[concentration_now][i]<<setw(20)<<output_more_states_concentration_index[concentration_now][i]<<endl;
                }
                
            }
            
            cout<<output_more_states_string.str();
            
            std::system("rm -r more_low_E_states");
            std::system("mkdir more_low_E_states");
            
            const char *path="./more_low_E_states/states_info.txt";
            std::ofstream file(path);
            file<<output_more_states_string.str();
            
            for (map<tuple<double,double>,vector<double> >::iterator it1=output_more_states_concentration_formation_energy.begin(); it1!=output_more_states_concentration_formation_energy.end(); it1++) {
                
                tuple<double,double> concentration_now=it1->first;
                vector<double> energy_vector=it1->second;
                int size_to_enumerate_over=energy_vector.size();
                
                for (int i=0; i<size_to_enumerate_over; i++) {
                    
                    string command_string="mkdir more_low_E_states/index"+lexical_cast<string>(output_more_states_concentration_index[concentration_now][i]);
                    std::system(command_string.c_str());
                    int index=output_more_states_concentration_index[concentration_now][i];
                    bool tern_alg=true;
                    print_poscar_out(output_more_states_periodicity_vector[index], output_more_states_unit_cell_vector[index], components,tern_alg);
                    
                    command_string="cp POSCAR_OUT more_low_E_states/index"+lexical_cast<string>(output_more_states_concentration_index[concentration_now][i]);
                    std::system(command_string.c_str());
                    
                    print_poscar_no_vacancy(output_more_states_periodicity_vector[index], output_more_states_unit_cell_vector[index], components,tern_alg);
                    
                    command_string="cp POSCAR more_low_E_states/index"+lexical_cast<string>(output_more_states_concentration_index[concentration_now][i]);
                    std::system(command_string.c_str());
                    
                }
                
            }
            
            
            //.... waits to be programmed
            1;
        }
        
        
        
        
    }
    else if (global_parameters.ternary_debug){
                cout<<"\n in ternary_debug";
//        cout<<"\n constant is "<<constant;

        mu.clear();
        mu.insert(global_parameters.mu_tern_first_group.begin(),global_parameters.mu_tern_first_group.end());
        mu.insert(global_parameters.mu_tern_second_group.begin(),global_parameters.mu_tern_second_group.end());

        
        map<set< tuple<int,int,int,int,int> >,double > J_tern_in;
        convert_J_input_tern_unaltered_to_J_tern_in(global_parameters.J_input_tern_unaltered, J_tern_in);

        
//        merging a few maps
        J_tern_in.insert(global_parameters.mu_tern_first_group.begin(),global_parameters.mu_tern_first_group.end());
        J_tern_in.insert(global_parameters.mu_tern_second_group.begin(),global_parameters.mu_tern_second_group.end());
        
//        cout<<"\n debug again check J_tern_in";
//        printmapfromsets(J_tern_in);
        
        
        map<int,int> components;
        int x_range, y_range, z_range;
        calculate_range_from_J(J_tern_in, x_range, y_range, z_range, components);
        
        
        vector<map< set<tuple<int,int,int,int,int> >, double> > cluster_type_vector;
        vector<double> formation_energy_vector;
        vector<double> mu1_vector, mu2_vector;
        map< set<tuple<int,int,int,int,int> >, double> cluster_type_temp;
        double formation_energy_temp=0;
        
        vector<tuple<double,double> > concentration_vector;
        map< set<tuple<int,int,int,int,int> >, double> J_fixed_part,J_in;
        vector<map< tuple<int,int,int,int>,int> >unit_cell_vector;
        vector<tuple<int,int,int,int,int,int> >periodicity_vector;
        
        vector<map< tuple<int,int,int,int>,int> > output_more_states_unit_cell_vector;
        vector<tuple<int,int,int,int,int,int> > output_more_states_periodicity_vector;
        vector<tuple<double,double> > output_more_states_concentration_vector;
        vector<int> output_more_states_super_cell_size_vector;
        vector<double> output_more_states_formation_energy_vector;
        
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_tern_in.begin(); it1!=J_tern_in.end(); it1++) {
            if (mu.count(it1->first)==0) {
                J_fixed_part[it1->first]=it1->second;
            }
            else{
                J_fixed_part[it1->first]=J_tern_in[it1->first]-mu[it1->first];
            }
        }
        
        
        //the 0 point
        formation_energy_temp=0;
        cluster_type_temp.clear();
        
        
        bool tern_alg=true;

        
        int temp_int=1;
        while (temp_int==1) {
            temp_int++;

            
            
            
            double first_chemical_potential=global_parameters.mu1;
            double second_chemical_potential=global_parameters.mu2;
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=global_parameters.mu_tern_first_group.begin(); it1!=global_parameters.mu_tern_first_group.end(); it1++) {
                mu[it1->first]=first_chemical_potential;

            }
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=global_parameters.mu_tern_second_group.begin(); it1!=global_parameters.mu_tern_second_group.end(); it1++) {
                mu[it1->first]=second_chemical_potential;
            }
            
            
            for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J_fixed_part.begin(); it1!=J_fixed_part.end(); it1++) {
                J_in[it1->first]=it1->second;
                if (mu.count(it1->first)==1) {
                    J_in[it1->first]+=mu[it1->first];
                }
            }

            run_solver(max_sites, J_in,
                       lowerboundclustertype, upperboundclustertype, cellrepresentation,
                       lower_bound, upper_bound, exact_lowerbound, unitcell, periodicity,
                       J_for_proof,
                       id,mu,mu_constant,work_with_mu,formation_energy_UB,formation_energy_LB,new_cluster_algorithm,use_new_pair_terms,use_new_triplet_terms,map_periodicity_to_spin,use_level_method,use_weighted_dual_average,global_parameters,
                       prec, num_loops,
                       basic_exact_mode, pseudo_mode, pseudo_mode_with_proof,
                       verbose, very_verbose, obscenely_verbose,input_PRIM_output_PRIMOUT_mode,limit_dimension,constant);
            
//            cout<<"\nlet me check what is formation energy_UB: "<<formation_energy_UB;
            
            tuple<double,double> concentration_temp;
            double formation_energy_temp;
            map<set<tuple<int, int, int, int, int> >, double> clustertype_out_temp;
            get_concentration_formation_energy_from_spin_and_J_in_ternary(periodicity, global_parameters.mu_tern_first_group, global_parameters.mu_tern_second_group, J_fixed_part, unitcell, constant, mu_constant, formation_energy_temp, concentration_temp, clustertype_out_temp);
            
            cout<<"\n periodicity is "<<periodicity;
            cout<<"\n energy including mu is "<<upper_bound;
            cout<<"\n the other way to compute formation energy is "<<formation_energy_temp;
            cout<<"\n the concentration is "<<concentration_temp;
            
            
            cout<<"\n the way i think about upper bound is : ";
            
            cout<<"\n";
            bool inside_indicator=false;
            
        };
    
    
    }
    
    bool kill_all_slave=true;
    broadcast(mpi_world,kill_all_slave,0);
    
    
    return 0;
}

void read_from_file(std::string &id,
                    int &max_sites,
                    map< set<tuple<int,int,int,int,int> >, double> &J,
                    double &prec,
                    int &num_loops,
                    bool &translation_algorithm,
                    bool &basic_exact_mode,
                    bool &pseudo_mode,
                    bool &pseudo_mode_with_proof,
                    bool &verbose,
                    bool &very_verbose,
                    bool &obscenely_verbose,
                    bool &input_PRIM_output_PRIMOUT_mode,
                    double &limit_dimension,
                    double &constant,
                    map< set<tuple<int,int,int,int,int> >, double> &mu,
                    bool &mu_translation_algorithm,
                    double &mu_constant,
                    bool &work_with_mu,
                    bool &scan_chemical_potential,
                    bool &new_cluster_algorithm,
                    int &use_new_pair_terms,
                    int &use_new_triplet_terms,
                    bool &output_more_states,
                    bool &output_states_below_hull,
                    double &how_much_lower,
                    double &output_states_how_sparse,
                    bool &use_level_method,
                    bool &use_weighted_dual_average,
                    solver_variable &global_parameters){

        map< set<tuple<int,int,int,int,int> >, double> J_negative1_positive1,mu_negative1_positive1;
        set< tuple<int,int,int,int,int> >setoftupletemp;
        double constant_term = 0;
        int solver_mode = 1;
        int verbosity = 1;

        std::ifstream t2("J_config.in");
        std::stringstream buffer2;
        buffer2 << t2.rdbuf();
        string J_config_file=buffer2.str();

        vector<string> J_config_line;
        split(J_config_file, '\n', J_config_line);
        for (vector<string>::iterator it=J_config_line.begin(); it!=J_config_line.end(); it++) {
            vector<string> J_config_line_segment;
            string temp=*it;
            split(temp, ' ', J_config_line_segment);
            vector<string>::iterator segment_iterator=J_config_line_segment.begin();
            string first_segment=*segment_iterator;
            vector<string> equal_left_right;
            split(first_segment, '=', equal_left_right);
            vector<string>::iterator equal_iterator=equal_left_right.begin();
            string equal_left=*equal_iterator;
            equal_iterator++;
            string equal_right=*equal_iterator;

            if (equal_left=="NSITES") {
                max_sites=lexical_cast<int>(equal_right) ;
            }else if (equal_left=="NLOOPS") {
                num_loops=lexical_cast<int>(equal_right) ;
            }else if (equal_left=="LABEL") {
                id=equal_right;
            }else if (equal_left=="PREC") {
                prec=lexical_cast<double>(equal_right) ;
                global_parameters.prec=lexical_cast<double>(equal_right);
            }else if (equal_left=="MODE_JPLUSMINUS") {
                translation_algorithm=lexical_cast<bool>(equal_right) ;
            }else if (equal_left=="MODE_SOLVER") {
                solver_mode=lexical_cast<int>(equal_right) ;
            }else if (equal_left=="MODE_VERBOSITY") {
                verbosity = lexical_cast<int>(equal_right) ;
            }
            else if (equal_left=="input_PRIM_output_PRIMOUT_mode") {
                input_PRIM_output_PRIMOUT_mode=lexical_cast<int>(equal_right) ;
            }
            else if (equal_left=="LIMIT_DIMENSION") {
                limit_dimension=lexical_cast<int>(equal_right) ;
            }
            else if (equal_left=="work_with_mu") {
                work_with_mu=lexical_cast<bool>(equal_right) ;
            }
            else if (equal_left=="scan_chemical_potential") {
                scan_chemical_potential=lexical_cast<bool>(equal_right) ;
            }
            else if (equal_left=="new_cluster_algorithm"){
                new_cluster_algorithm=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="use_new_pair_terms"){
                use_new_pair_terms=lexical_cast<int>(equal_right);
            }
            else if (equal_left=="use_new_triplet_terms"){
                use_new_triplet_terms=lexical_cast<int>(equal_right);
            }
            else if (equal_left=="output_more_states"){
                output_more_states=lexical_cast<int>(equal_right);
            }
            else if (equal_left=="output_states_below_hull"){
                output_states_below_hull=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="how_much_lower"){
                how_much_lower=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="output_states_how_sparse"){
                output_states_how_sparse=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="use_level_method"){
                use_level_method=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="use_weighted_dual_average"){
                use_weighted_dual_average=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="use_gradual_introduction_new_variables"){
                global_parameters.use_gradual_introduction_new_variables=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="do_monte_carlo"){
                global_parameters.do_monte_carlo=lexical_cast<bool>(equal_right);
            }            else if (equal_left=="monte_carlo_periodicity"){
                string monte_periodicity_string=equal_right;
                vector<string> periodicity_segment;
                split(monte_periodicity_string, ',', periodicity_segment);
                int monte_a0=lexical_cast<int>(periodicity_segment[0]) ;
                int monte_a1=lexical_cast<int>(periodicity_segment[1]) ;
                int monte_a2=lexical_cast<int>(periodicity_segment[2]) ;
                int monte_a3=lexical_cast<int>(periodicity_segment[3]) ;
                int monte_a4=lexical_cast<int>(periodicity_segment[4]) ;
                int monte_a5=lexical_cast<int>(periodicity_segment[5]) ;
                global_parameters.monte_carlo_periodicity=make_tuple(monte_a0,monte_a1,monte_a2,monte_a3,monte_a4,monte_a5);
            }
            else if (equal_left=="initial_T"){
                global_parameters.initial_T=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="final_T"){
                global_parameters.final_T=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="T_decrement"){
                global_parameters.T_decrement=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="each_T_how_many_runs"){
                global_parameters.each_T_how_many_runs=lexical_cast<long long>(equal_right);
            }
            else if (equal_left=="dedicated1D"){
                global_parameters.dedicated1D=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="ternary_alg"){
                global_parameters.ternary_alg=lexical_cast<bool>(equal_right);
            }
            
            else if (equal_left=="ternary_output_states_above_hull"){
                global_parameters.ternary_output_states_above_hull=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="ternary_x_min"){
                global_parameters.ternary_x_min=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="ternary_x_max"){
                global_parameters.ternary_x_max=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="ternary_y_min"){
                global_parameters.ternary_y_min=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="ternary_y_max"){
                global_parameters.ternary_y_max=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="ternary_z_min"){
                global_parameters.ternary_z_min=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="ternary_z_max"){
                global_parameters.ternary_z_max=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="ternary_output_states_above_hull_gap"){
                global_parameters.ternary_output_states_above_hull_gap=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="ternary_output_states_above_hull_ceil"){
                global_parameters.ternary_output_states_above_hull_ceil=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="ternary_debug"){
                global_parameters.ternary_debug=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="mu1"){
                global_parameters.mu1=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="mu2"){
                global_parameters.mu2=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="use_incomplete_solver"){
                global_parameters.use_incomplete_solver=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="which_incomplete_solver"){
                global_parameters.which_incomplete_solver=equal_right;
            }
            else if (equal_left=="incomplete_time_cap"){
                global_parameters.incomplete_time_cap=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="N_to_start_incomplete"){
                global_parameters.N_to_start_incomplete=lexical_cast<int>(equal_right);
            }
            else if (equal_left=="Major_periodicity"){
                global_parameters.Major_periodicity=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="activate_large_cell_algo_binary"){
                global_parameters.activate_large_cell_algo_binary=lexical_cast<bool>(equal_right);
            }
            else if (equal_left=="large_cell_binary_cubic_periodicity"){
                global_parameters.large_cell_binary_cubic_periodicity=lexical_cast<long long>(equal_right);
            }
            else if (equal_left=="large_cell_binary_mu_init"){
                global_parameters.large_cell_binary_mu_init=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="large_cell_binary_mu_step"){
                global_parameters.large_cell_binary_mu_step=lexical_cast<double>(equal_right);
            }
            else if (equal_left=="large_cell_binary_mu_final"){
                global_parameters.large_cell_binary_mu_final=lexical_cast<double>(equal_right);
            }
        }
    
    

        if (solver_mode == 0) {
            basic_exact_mode = true;
            pseudo_mode = false;
            pseudo_mode = false;
        }else if (solver_mode == 1) {
            basic_exact_mode = false;
            pseudo_mode = true;
            pseudo_mode_with_proof = false;
        }else if (solver_mode == 2) {
            basic_exact_mode = false;
            pseudo_mode = true;
            pseudo_mode_with_proof = true;
        }else{
            cout << "Invalid solver mode given. Exiting." << endl;
            exit(1);
        }

        if (verbosity == 0) {
            verbose = false;
            very_verbose = false;
            obscenely_verbose = false;
        }else if (verbosity == 1) {
            verbose = true;
            very_verbose = false;
            obscenely_verbose = false;
        }else if (verbosity == 2) {
            verbose = true;
            very_verbose = true;
            obscenely_verbose = false;
        }else if (verbosity == 3) {
            verbose = true;
            very_verbose = true;
            obscenely_verbose = true;
        }else {
            cout << "Invalid verbosity mode given. Exiting." << endl;
            exit(1);
        }

        std::ifstream t1("J_in.in");
        std::stringstream buffer1;
        buffer1 << t1.rdbuf();
        string J_in_file=buffer1.str();

        vector<string> J_in_line;
        split(J_in_file, '\n', J_in_line);
        int line_format_indicator=0;
        for (vector<string>::iterator it=J_in_line.begin(); it!=J_in_line.end(); it++){
            string this_line=*it;
            if (it==J_in_line.begin()){
                vector<string> segment;
                split(this_line, ' ', segment);
                vector<string>::iterator segment_iterator=segment.begin();
                string temp_segment=*segment_iterator;
                if (temp_segment=="Constant"){
                    segment_iterator++;
                    constant_term=lexical_cast<double>(*segment_iterator);
                }
            }else{
                if (line_format_indicator==0) {
                    vector<string> segment;
                    split(this_line, ' ', segment);
                    vector<string>::iterator segment_iterator=segment.begin();
                    string temp_segment=*segment_iterator;
                    if (temp_segment=="Cluster") {
                        line_format_indicator=1;
                    }
                }else if (line_format_indicator==1){
                    setoftupletemp.clear();
                    vector<string> segment;
                    split(this_line, ' ', segment);

                    for (vector<string>::iterator it2=segment.begin(); it2!=segment.end(); it2++) {
                        string this_segment=*it2;
                        vector<string> number_vector;
                        split(this_segment, ',', number_vector);
                        int l1,l2,l3,l4,l5;
                        l1=lexical_cast<int>(number_vector[0]);
                        l2=lexical_cast<int>(number_vector[1]);
                        l3=lexical_cast<int>(number_vector[2]);
                        l4=lexical_cast<int>(number_vector[3]);
                        l5=lexical_cast<int>(number_vector[4]);
                        setoftupletemp.insert(make_tuple(l1,l2,l3,l4,l5));
                    }

                    line_format_indicator=2;
                }else if (line_format_indicator==2) {
                    vector<string> segment;
                    split(this_line, '=', segment);

                    if (translation_algorithm){
                        J_negative1_positive1[setoftupletemp]=lexical_cast<double>(segment[1]);
                    }else{
                        J[setoftupletemp]=lexical_cast<double>(segment[1]);
                    }

                    line_format_indicator=0;
                }
            }
        }
    
    
    
    

        if(translation_algorithm){
            convertJ_neg1_pos1toJ(J_negative1_positive1, J, constant_term);
        }
    
    
    constant=constant_term;
    
    if (work_with_mu) {
        std::ifstream t2("mu_config.in");
        std::stringstream buffer2;
        buffer2 << t2.rdbuf();
        string J_config_file=buffer2.str();
        
        vector<string> mu_config_line;
        split(J_config_file, '\n', mu_config_line);
        for (vector<string>::iterator it=mu_config_line.begin(); it!=mu_config_line.end(); it++) {
            vector<string> mu_config_line_segment;
            string temp=*it;
            split(temp, ' ', mu_config_line_segment);
            vector<string>::iterator segment_iterator=mu_config_line_segment.begin();
            string first_segment=*segment_iterator;
            vector<string> equal_left_right;
            split(first_segment, '=', equal_left_right);
            vector<string>::iterator equal_iterator=equal_left_right.begin();
            string equal_left=*equal_iterator;
            equal_iterator++;
            string equal_right=*equal_iterator;
            
            if (equal_left=="MODE_JPLUSMINUS") {
                mu_translation_algorithm =lexical_cast<bool>(equal_right) ;
            }
        }
        
        
        
        

        std::ifstream t1("mu_in.in");
        std::stringstream buffer1;
        buffer1 << t1.rdbuf();
        string mu_in_file=buffer1.str();
        
        vector<string> mu_in_line;
        split(mu_in_file, '\n', mu_in_line);
        int line_format_indicator=0;
        for (vector<string>::iterator it=mu_in_line.begin(); it!=mu_in_line.end(); it++){
            string this_line=*it;
            if (it==mu_in_line.begin()){
                vector<string> segment;
                split(this_line, ' ', segment);
                vector<string>::iterator segment_iterator=segment.begin();
                string temp_segment=*segment_iterator;
                if (temp_segment=="Constant"){
                    segment_iterator++;
                    mu_constant=lexical_cast<double>(*segment_iterator);
                }
            }else{
                if (line_format_indicator==0) {
                    vector<string> segment;
                    split(this_line, ' ', segment);
                    vector<string>::iterator segment_iterator=segment.begin();
                    string temp_segment=*segment_iterator;
                    if (temp_segment=="Cluster") {
                        line_format_indicator=1;
                    }
                }else if (line_format_indicator==1){
                    setoftupletemp.clear();
                    vector<string> segment;
                    split(this_line, ' ', segment);
                    
                    for (vector<string>::iterator it2=segment.begin(); it2!=segment.end(); it2++) {
                        string this_segment=*it2;
                        vector<string> number_vector;
                        split(this_segment, ',', number_vector);
                        int l1,l2,l3,l4,l5;
                        l1=lexical_cast<int>(number_vector[0]);
                        l2=lexical_cast<int>(number_vector[1]);
                        l3=lexical_cast<int>(number_vector[2]);
                        l4=lexical_cast<int>(number_vector[3]);
                        l5=lexical_cast<int>(number_vector[4]);
                        setoftupletemp.insert(make_tuple(l1,l2,l3,l4,l5));
                    }
                    
                    line_format_indicator=2;
                }else if (line_format_indicator==2) {
                    vector<string> segment;
                    split(this_line, '=', segment);
                    
                    if (translation_algorithm){
                        mu_negative1_positive1[setoftupletemp]=lexical_cast<double>(segment[1]);
                    }else{
                        mu[setoftupletemp]=lexical_cast<double>(segment[1]);
                    }
                    
                    line_format_indicator=0;
                }
            }
        }

        
        if(mu_translation_algorithm){
            convertJ_neg1_pos1toJ(mu_negative1_positive1, mu, mu_constant);
        }
        

        
        for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
            J[it1->first]+=it1->second;
        }
        constant+=mu_constant;
    }
    
    if (global_parameters.ternary_alg||global_parameters.ternary_debug) {
        constant=0;
        
        set< tuple<int,int,int,int,int,int> >setoftupletemp_tern;
        
        std::ifstream t1("J_in_tern_casm.in");
        std::stringstream buffer1;
        buffer1 << t1.rdbuf();
        string J_in_tern_casm_file=buffer1.str();
        
        vector<string> J_in_tern_casm_line;
        split(J_in_tern_casm_file, '\n', J_in_tern_casm_line);
        int line_format_indicator=0;
        for (vector<string>::iterator it=J_in_tern_casm_line.begin(); it!=J_in_tern_casm_line.end(); it++){
            string this_line=*it;
            if (it==J_in_tern_casm_line.begin()){
                vector<string> segment;
                split(this_line, ' ', segment);
                vector<string>::iterator segment_iterator=segment.begin();
                string temp_segment=*segment_iterator;
                
                if (temp_segment=="Constant"){
                    segment_iterator++;
                    constant=lexical_cast<double>(*segment_iterator);
                }
            }
            else{
                if (line_format_indicator==0) {
                    vector<string> segment;
                    split(this_line, ' ', segment);
                    vector<string>::iterator segment_iterator=segment.begin();
                    string temp_segment=*segment_iterator;
                    if (temp_segment=="Cluster") {
                        line_format_indicator=1;
                    }
                }else if (line_format_indicator==1){
                    setoftupletemp_tern.clear();
                    vector<string> segment;
                    split(this_line, ' ', segment);
                    
                    for (vector<string>::iterator it2=segment.begin(); it2!=segment.end(); it2++) {
                        string this_segment=*it2;
                        vector<string> number_vector;
                        split(this_segment, ',', number_vector);
                        int l1,l2,l3,l4,l5,l6;
                        l1=lexical_cast<int>(number_vector[0]);
                        l2=lexical_cast<int>(number_vector[1]);
                        l3=lexical_cast<int>(number_vector[2]);
                        l4=lexical_cast<int>(number_vector[3]);
                        l5=lexical_cast<int>(number_vector[4]);
                        l6=lexical_cast<int>(number_vector[5]);
                        setoftupletemp_tern.insert(make_tuple(l1,l2,l3,l4,l5,l6));
                    }
                    
                    line_format_indicator=2;
                }else if (line_format_indicator==2) {
                    vector<string> segment;
                    split(this_line, '=', segment);

                    global_parameters.J_input_tern_unaltered[setoftupletemp_tern]+=lexical_cast<double>(segment[1]);
                    
//                    if (translation_algorithm){
//                        J_negative1_positive1[setoftupletemp]=lexical_cast<double>(segment[1]);
//                    }else{
//                        J[setoftupletemp]=lexical_cast<double>(segment[1]);
//                    }
                    
                    line_format_indicator=0;
                }
            }
        }

        
        
        std::ifstream t3("mu_in_tern_casm.in");
        std::stringstream buffer3;
        buffer3 << t3.rdbuf();
        string mu_tern_in_file=buffer3.str();
        
        vector<string> mu_tern_in_line;
        split(mu_tern_in_file, '\n', mu_tern_in_line);
        line_format_indicator=0;
        int group_index=0;
        for (vector<string>::iterator it=mu_tern_in_line.begin(); it!=mu_tern_in_line.end(); it++){
            string this_line=*it;
            
            if (it==mu_tern_in_line.begin()){
                vector<string> segment;
                split(this_line, ' ', segment);
                vector<string>::iterator segment_iterator=segment.begin();
                string temp_segment=*segment_iterator;
                if (temp_segment=="Constant"){
                    segment_iterator++;
                    mu_constant=lexical_cast<double>(*segment_iterator);
                }
            }else{
                if (line_format_indicator==0) {
                    vector<string> segment;
                    split(this_line, ' ', segment);
                    vector<string>::iterator segment_iterator=segment.begin();
                    string temp_segment=*segment_iterator;
                    if (temp_segment=="Cluster") {
                        line_format_indicator=1;
                    }
                    
                }
                else if (line_format_indicator==1){
                    setoftupletemp.clear();
                    vector<string> segment;
                    split(this_line, ' ', segment);
                    assert(segment[0]=="Group");
                    group_index=lexical_cast<int>(segment[1]);
                    
                    line_format_indicator=2;
                }
                else if (line_format_indicator==2){
                    
                    setoftupletemp.clear();
                    vector<string> segment;
                    split(this_line, ' ', segment);
                    
                    for (vector<string>::iterator it2=segment.begin(); it2!=segment.end(); it2++) {
                        string this_segment=*it2;
                        vector<string> number_vector;
                        split(this_segment, ',', number_vector);
                        int l1,l2,l3,l4,l5;
                        l1=lexical_cast<int>(number_vector[0]);
                        l2=lexical_cast<int>(number_vector[1]);
                        l3=lexical_cast<int>(number_vector[2]);
                        l4=lexical_cast<int>(number_vector[3]);
                        l5=lexical_cast<int>(number_vector[4]);
                        setoftupletemp.insert(make_tuple(l1,l2,l3,l4,l5));
                    }
                    
                    line_format_indicator=3;
                }
                else if (line_format_indicator==3) {
                    assert(group_index==1||group_index==2);
                    vector<string> segment;
                    split(this_line, '=', segment);
                    
                    
                    if(group_index==1)
                        global_parameters.mu_tern_first_group[setoftupletemp]=lexical_cast<double>(segment[1]);
                    
                    if(group_index==2)
                        global_parameters.mu_tern_second_group[setoftupletemp]=lexical_cast<double>(segment[1]);
                    
                    line_format_indicator=0;
                }
            }
        }
        
        constant+=mu_constant;
        
//        cout<<"\ndebug here print global_parameters.J_input_tern_unaltered";
//        printmapfromsets(global_parameters.J_input_tern_unaltered);
//        cout<<"\n debug here to print mu_tern_first_group";
//        printmapfromsets(global_parameters.mu_tern_first_group);
//        cout<<"\n debug here to print mu_tern_second_group";
//        printmapfromsets(global_parameters.mu_tern_second_group);
        
    }
    
    
    
    
    
}

