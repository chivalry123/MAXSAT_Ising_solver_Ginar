#include "tools.h"

using namespace std;
using namespace ::boost::tuples;
using namespace ::boost;

void matchnumber(map<tuple<int,int,int,int,int,int>, double >& dictionary,double number,set<tuple<int,int,int,int,int,int> > &list){
    list.clear();
    for (map<tuple<int,int,int,int,int,int>, double >::iterator ait=dictionary.begin(); ait!=dictionary.end(); ait++)
    {
        if (fabs(ait->second-number)<fabs(number)*1e-5) {
            list.insert(ait->first);
        }
    }
}

long long myPow(long long x, long long p){
    if (p == 0) return 1;
    if (p == 1) return x;

    long tmp = myPow(x, p/2);
    if (p%2 == 0) return tmp * tmp;
    else return x * tmp * tmp;
}

std::string exec(std::string command) {
    //cout << "\n\nEvaluating command: " << command << "\n\n" << endl;
    const char* cmd = command.c_str();
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[256];
    std::string result = "";
    while(!feof(pipe)) {
    	if(fgets(buffer, 256, pipe) != NULL)
    		result += buffer;
    }
    pclose(pipe);
    return result;
}

std::vector<std::string> split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if(item.length() > 0)
            elems.push_back(item);
    }
    return elems;
}

void split_is(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if(item.length() > 0)
            elems.push_back(item);
    }
}

//void split_is(const std::string &s, string delim_regex, std::vector<std::string> &elems) {
//    // Terribly slow and inefficient, sorry. Should probably fix this at some point
//    boost::algorithm::split_regex(elems, s, regex(delim_regex));
////    regex(delim_regex);
//    std::vector<std::string> non_zero_elems;
//    for(int e = 0; e < elems.size(); e++){
//        if (elems[e].length() > 0){
//            non_zero_elems.push_back(elems[e]);
//        }
//    }
//    elems = non_zero_elems;
//}

std::string to_string_is(int n){
    std::string str = boost::lexical_cast<std::string>(n);
    return str;
}

std::string to_string_is(double d){
    std::string str = boost::lexical_cast<std::string>(d);
    return str;
}

int stoi_is(std::string s){
    return boost::lexical_cast<int>(s);
}

double stod_is(std::string s){
    return boost::lexical_cast<double>(s);
}

void convertspintostring(map<tuple<int,int,int,int>,int>& spin, string &spinstring) {
    set<tuple<int,int,int,int> > keylist;
    int imax = getimax<0>(spin);
    int jmax = getimax<1>(spin);
    int zmax = getimax<2>(spin);
    int location_max = getimax<3>(spin);
    int location_min = getimin<3>(spin);
    spinstring = "";
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
}

void printblock(map<tuple<int,int,int,int>,int> &thisblock){
    string tempstring3;
    convertspintostring(thisblock, tempstring3);
    cout << tempstring3 << "\n";
}

void printblocklist(vector<map<tuple<int,int,int,int>,int> > &thisblocklist){
    for(vector<map<tuple<int,int,int,int>,int> >::iterator it=thisblocklist.begin();it!=thisblocklist.end();it++ ){
        cout << "\n-----";
        printblock(*it);
    }
}

void calculate_formation_energy(map< set<tuple<int,int,int,int,int> >, double> J,map< set<tuple<int,int,int,int,int> >, double> mu, map< set<tuple<int,int,int,int,int> >, double>cluster_type_temp, double constant, double mu_constant ,double &formation_energy_temp)
{
    formation_energy_temp=0;
    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=J.begin(); it1!=J.end(); it1++) {
        formation_energy_temp+=cluster_type_temp[it1->first]*it1->second;
    }
    formation_energy_temp+=constant;
    
    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
        formation_energy_temp-=cluster_type_temp[it1->first]*it1->second;
    }
    formation_energy_temp-=mu_constant;
    
}

void print_poscar_out(tuple<int,int,int,int,int,int> periodicity,map< tuple<int,int,int,int>,int> unitcell,map<int,int> components,bool tern_alg)
{
    
    cout<<"\nprinting POSCAROUT\n";
    std::ifstream t3("PRIM");
    std::stringstream poscar_in;
    poscar_in << t3.rdbuf();
    string poscar_in_string=poscar_in.str();
    
    ofstream poscar_out;
    string poscar_out_filename="POSCAR_OUT";
    poscar_out.open (poscar_out_filename.c_str());
    
    vector<string> poscar_in_line;
    split(poscar_in_string, '\n', poscar_in_line);
    
    
    vector<string> primitive_x_in_line;
    vector<string> primitive_y_in_line;
    vector<string> primitive_z_in_line;
    vector<double> primitive_x_double;
    vector<double> primitive_y_double;
    vector<double> primitive_z_double;
    
    //        cout<<poscar_in_line[2]<<endl<<poscar_in_line[3]<<endl<<poscar_in_line[4]<<endl;
    split(poscar_in_line[2], ' ', primitive_x_in_line);
    split(poscar_in_line[3], ' ', primitive_y_in_line);
    split(poscar_in_line[4], ' ', primitive_z_in_line);
    primitive_x_double.push_back(lexical_cast<double >(primitive_x_in_line[0])) ;
    primitive_x_double.push_back(lexical_cast<double >(primitive_x_in_line[1])) ;
    primitive_x_double.push_back(lexical_cast<double >(primitive_x_in_line[2])) ;
    primitive_y_double.push_back(lexical_cast<double >(primitive_y_in_line[0])) ;
    primitive_y_double.push_back(lexical_cast<double >(primitive_y_in_line[1])) ;
    primitive_y_double.push_back(lexical_cast<double >(primitive_y_in_line[2])) ;
    primitive_z_double.push_back(lexical_cast<double >(primitive_z_in_line[0])) ;
    primitive_z_double.push_back(lexical_cast<double >(primitive_z_in_line[1])) ;
    primitive_z_double.push_back(lexical_cast<double >(primitive_z_in_line[2])) ;
    
    //        printvector(primitive_x_double);
    //        printvector(primitive_y_double);
    //        printvector(primitive_z_double);
    
    vector<double> gs_primitive_x_double;
    vector<double> gs_primitive_y_double;
    vector<double> gs_primitive_z_double;
    
    gs_primitive_x_double.push_back(primitive_x_double[0]*periodicity.get<0>());
    gs_primitive_x_double.push_back(primitive_x_double[1]*periodicity.get<0>());
    gs_primitive_x_double.push_back(primitive_x_double[2]*periodicity.get<0>());
    gs_primitive_y_double.push_back(primitive_y_double[0]*periodicity.get<2>()+primitive_x_double[0]*periodicity.get<1>());
    gs_primitive_y_double.push_back(primitive_y_double[1]*periodicity.get<2>()+primitive_x_double[1]*periodicity.get<1>());
    gs_primitive_y_double.push_back(primitive_y_double[2]*periodicity.get<2>()+primitive_x_double[2]*periodicity.get<1>());
    gs_primitive_z_double.push_back(primitive_z_double[0]*periodicity.get<5>()+primitive_y_double[0]*periodicity.get<4>()+primitive_x_double[0]*periodicity.get<3>());
    gs_primitive_z_double.push_back(primitive_z_double[1]*periodicity.get<5>()+primitive_y_double[1]*periodicity.get<4>()+primitive_x_double[1]*periodicity.get<3>());
    gs_primitive_z_double.push_back(primitive_z_double[2]*periodicity.get<5>()+primitive_y_double[2]*periodicity.get<4>()+primitive_x_double[2]*periodicity.get<3>());
    
    
    vector<string> looping_numbers;
    
    split(poscar_in_line[5], ' ', looping_numbers);
    
    //        cout<<endl<<poscar_in_line[5];
    //        cout<<"\nindicator2\n";
    //
    //        printvector(looping_numbers);
    //        cout<<looping_numbers.size();
    //
    //        vector<string> line_of_a_term;
    //        split(poscar_in_line[8], ' ', line_of_a_term);
    //        cout<<"\n line numbers for line of a term is "<<line_of_a_term.size()  ;
    
    int total_lines=0;
    for (int i=0; i<looping_numbers.size(); i++) {
        total_lines+=lexical_cast<int>(looping_numbers[i]);
    }
    //        cout<<"\n total lines are "<<total_lines;
    
    
    vector< vector<double> > active_coordinates;
    vector< vector<string> > active_components;
    vector< vector<double> > passive_coordinates;
    vector< vector<string> > passive_components;
    
    for (int line_number=7 ; line_number<7+total_lines; line_number++) {
        int index=line_number-6;
        vector<string> read_line;
        split(poscar_in_line[line_number], ' ' , read_line);
        
        if (read_line.size()>4) {
            vector<double> vector_temp;
            vector_temp.push_back(lexical_cast<double>(read_line[0]));
            vector_temp.push_back(lexical_cast<double>(read_line[1]));
            vector_temp.push_back(lexical_cast<double>(read_line[2]));
            vector<string> vector_component;
            for (int temp_index=3; temp_index<read_line.size(); temp_index++) {
                vector_component.push_back(read_line[temp_index]);
            }
            active_coordinates.push_back(vector_temp);
            active_components.push_back(vector_component);
        }
        else if(read_line.size()==4) {
            
            vector<double> vector_temp;
            vector_temp.push_back(lexical_cast<double>(read_line[0]));
            vector_temp.push_back(lexical_cast<double>(read_line[1]));
            vector_temp.push_back(lexical_cast<double>(read_line[2]));
            vector<string> vector_component;
            for (int temp_index=3; temp_index<read_line.size(); temp_index++) {
                vector_component.push_back(read_line[temp_index]);
            }
            passive_coordinates.push_back(vector_temp);
            passive_components.push_back(vector_component);
            
        }
        
    }
    
    poscar_out << fixed << setprecision(8);
    //        cout << fixed << setprecision(8);
    poscar_out<<poscar_in_line[0]<<endl<<poscar_in_line[1]<<endl<<gs_primitive_x_double[0]<<" "<<gs_primitive_x_double[1]<<" "<<gs_primitive_x_double[2]<<endl<<gs_primitive_y_double[0]<<" "<<gs_primitive_y_double[1]<<" "<<gs_primitive_y_double[2]<<endl<<gs_primitive_z_double[0]<<" "<<gs_primitive_z_double[1]<<" "<<gs_primitive_z_double[2]<<endl;
    
    int a0,a1,a2,a3,a4,a5;
    a0=periodicity.get<0>();
    a1=periodicity.get<1>();
    a2=periodicity.get<2>();
    a3=periodicity.get<3>();
    a4=periodicity.get<4>();
    a5=periodicity.get<5>();
    
    vector<vector<double> > coordinate_in_poscar_out_before_transformation_sorting;
    vector<string> species_in_poscar_out_before_transformation_sorting;
    
    for (int i=1; i<=a0; i++) {
        for (int j=1; j<=a2; j++) {
            for (int k=1; k<=a5; k++) {
                for (int p=1; p<=components.size(); p++) {
                    int component_index=unitcell[make_tuple(i,j,k,p)];
                    vector<double> position_before_transformation;
                    position_before_transformation.push_back(i-1+active_coordinates[p-1][0]);
                    position_before_transformation.push_back(j-1+active_coordinates[p-1][1]);
                    position_before_transformation.push_back(k-1+active_coordinates[p-1][2]);
                    coordinate_in_poscar_out_before_transformation_sorting.push_back(position_before_transformation);
                    if (!tern_alg) {
                        species_in_poscar_out_before_transformation_sorting.push_back(active_components[p-1][components[p]-component_index-1]);
                    }
                    else if(tern_alg){
                        if (component_index==0)
                            species_in_poscar_out_before_transformation_sorting.push_back(active_components[p-1][1]);
                        
                        if (component_index==1)
                            species_in_poscar_out_before_transformation_sorting.push_back(active_components[p-1][0]);
                        
                        if (component_index==2)
                            species_in_poscar_out_before_transformation_sorting.push_back(active_components[p-1][2]);
                        
                    }

                }
                
                for (int passive_p=0; passive_p<passive_coordinates.size(); passive_p++) {
                    vector<double> position_before_transformation;
                    position_before_transformation.push_back(i-1+passive_coordinates[passive_p][0]);
                    position_before_transformation.push_back(j-1+passive_coordinates[passive_p][1]);
                    position_before_transformation.push_back(k-1+passive_coordinates[passive_p][2]);
                    coordinate_in_poscar_out_before_transformation_sorting.push_back(position_before_transformation);
                    species_in_poscar_out_before_transformation_sorting.push_back(passive_components[passive_p][0]);
                }
            }
        }
    }
    
    //        cout<<"\n try printing the thing";
    //        for (int temp=0; temp<coordinate_in_poscar_out_before_transformation_sorting.size(); temp++) {
    //            cout<<endl;
    //            for (int i=0; i<3; i++) {
    //                cout<< coordinate_in_poscar_out_before_transformation_sorting[temp][i];
    //                cout<< " ";
    //            }
    //            cout<<species_in_poscar_out_before_transformation_sorting[temp];
    //        }
    
    vector<vector<double> > coordinate_in_poscar_out_after_transformation_before_sorting;
    vector<string> species_in_poscar_out_after_transformation_before_sorting=species_in_poscar_out_before_transformation_sorting;
    
    
    //        cout<<"\n try converting the thing";
    for (int temp=0; temp<coordinate_in_poscar_out_before_transformation_sorting.size(); temp++) {
        double x,y,z;
        x=coordinate_in_poscar_out_before_transformation_sorting[temp][0];
        y=coordinate_in_poscar_out_before_transformation_sorting[temp][1];
        z=coordinate_in_poscar_out_before_transformation_sorting[temp][2];
        vector<double> temp_vec;
        temp_vec.push_back(x/a0-y/a2*a1/a0+z*a1*a4/a0/a2/a5-z/a5*a3/a0);
        temp_vec.push_back(y/a2-z/a5*a4/a2);
        temp_vec.push_back(z/a5);
        coordinate_in_poscar_out_after_transformation_before_sorting.push_back(temp_vec);
    }
    //
    //        cout<<"\n try printing the thing after transformation";
    //        for (int temp=0; temp<coordinate_in_poscar_out_after_transformation_before_sorting.size(); temp++) {
    //            cout<<endl;
    //            for (int i=0; i<3; i++) {
    //                cout<< coordinate_in_poscar_out_after_transformation_before_sorting[temp][i];
    //                cout<< " ";
    //            }
    //            cout<<species_in_poscar_out_after_transformation_before_sorting[temp];
    //        }
    
    map<string, vector<vector<double> > > reverse_species_coordinate;
    for (int temp=0; temp<coordinate_in_poscar_out_after_transformation_before_sorting.size(); temp++) {
        
        reverse_species_coordinate[species_in_poscar_out_after_transformation_before_sorting[temp]].push_back(coordinate_in_poscar_out_after_transformation_before_sorting[temp]);
    }
    
    
    vector<vector<double> > coordinate_in_poscar_out_after_transformation_after_sorting;
    vector<string> species_in_poscar_out_after_transformation_after_sorting=species_in_poscar_out_before_transformation_sorting;
    for (map<string, vector<vector<double> > >::iterator it1=reverse_species_coordinate.begin(); it1!=reverse_species_coordinate.end(); it1++) {
        poscar_out<<it1->first<<" ";
    }
    poscar_out<<endl;
    for (map<string, vector<vector<double> > >::iterator it1=reverse_species_coordinate.begin(); it1!=reverse_species_coordinate.end(); it1++) {
        poscar_out<<(it1->second).size()<<" ";
    }
    poscar_out<<endl;
    poscar_out<<poscar_in_line[6];
    poscar_out<<endl;
    for (map<string, vector<vector<double> > >::iterator it1=reverse_species_coordinate.begin(); it1!=reverse_species_coordinate.end(); it1++) {
        vector<vector<double> > coordinates_here=it1->second;
        for (int temp=0; temp<coordinates_here.size(); temp++) {
            vector<double> vector_here=coordinates_here[temp];
            poscar_out<<vector_here[0]-floor(vector_here[0])<<" "<<vector_here[1]-floor(vector_here[1])<<" "<<vector_here[2]-floor(vector_here[2])<<" "<<it1->first<<endl;
        }
    }
}

void convert_to_starting_point_cluster(set<tuple<int,int,int,int,int> > prototype_set,set<tuple<int,int,int,int,int> > &returned_equivalent_set)
{

    
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
    int translation_z = 1-z_min;
    int translation_y = 1-y_min;
    int translation_x = 1-x_min;
    
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
    
   returned_equivalent_set=temp_equivalent_set;
    
}



void get_concentration_formation_energy_from_spin_and_J_in(tuple<int,int,int,int,int,int> periodicity_now,map< set<tuple<int,int,int,int,int> >, double > mu,map< set<tuple<int,int,int,int,int> >, double > J_fixed_part, map<tuple<int,int,int,int>,int>  spin_now, double constant, double mu_constant,double&formation_energy_out, double &concentration_out)
{
    
    map< set<tuple<int,int,int,int,int> >, double > clustertype;
    
    int a0=periodicity_now.get<0>();
    int a1=periodicity_now.get<1>();
    int a2=periodicity_now.get<2>();
    int a3=periodicity_now.get<3>();
    int a4=periodicity_now.get<4>();
    int a5=periodicity_now.get<5>();
    
    for (map< set<tuple<int,int,int,int,int> >, double >::iterator it = J_fixed_part.begin();
         it != J_fixed_part.end(); it++){
        set<tuple<int,int,int,int,int> > thisset = it->first;
        clustertype[thisset] = 0;
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
                        if (spin_now[make_tuple(i + x - 1, j + y - 1, k + z - 1, position)] != var) {
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
                        clustertype[thisset] += 1.0/(a0*a2*a5);
                    }
                }
            }
        }
    }
    
    
    double formation_energy=0;
    for (map< set<tuple<int,int,int,int,int> >, double >::iterator it = J_fixed_part.begin();it != J_fixed_part.end(); it++) {
        formation_energy+=it->second*clustertype[it->first];
    }
    formation_energy+=constant;
    formation_energy-=mu_constant;
    
//    cout<<"\ndebug here what is formation energy: "<<setw(20) << std::fixed << std::setprecision(8)<<formation_energy;
    
    double concentration_for_this_periodicity=0;
    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu.begin(); it1!=mu.end(); it1++) {
        if (clustertype.count(it1->first)==1) {
            concentration_for_this_periodicity+=clustertype[it1->first];
        }
        else{
            perror("error please check, error id 52919sdi");
            return ;
        }
    }
    concentration_for_this_periodicity=concentration_for_this_periodicity/mu.size();
    
    formation_energy_out=formation_energy;
    concentration_out=concentration_for_this_periodicity;
//    cout<<" concentration_for_this_periodicity is "<<concentration_for_this_periodicity;
}

void evaluate_energy_at_hull(double concentration_now, map<double,double> hull_map, double &energy_at_hull)
{
    
    for (map<double, double>::iterator it1=hull_map.begin(); it1!=hull_map.end(); it1++) {
        double left_concentration=it1->first;
        double left_formation_energy=it1->second;
        double right_concentration;
        double right_formation_energy;
        it1++;
        right_concentration=it1->first;
        right_formation_energy=it1->second;
        it1--;
        if (concentration_now>left_concentration-1e-8&&concentration_now<1e-8+right_concentration) {
            energy_at_hull=(right_concentration-concentration_now)/(right_concentration-left_concentration)*left_formation_energy+(concentration_now-left_concentration)/(right_concentration-left_concentration)*right_formation_energy;
            break;
        }
        
        
    }
    
}



void print_poscar_no_vacancy(tuple<int,int,int,int,int,int> periodicity,map< tuple<int,int,int,int>,int> unitcell,map<int,int> components,bool tern_alg)
{
    
    cout<<"\nprinting POSCAR\n";
    std::ifstream t3("PRIM");
    std::stringstream poscar_in;
    poscar_in << t3.rdbuf();
    string poscar_in_string=poscar_in.str();
    
    ofstream poscar_out;
    string poscar_out_filename="POSCAR";
    poscar_out.open (poscar_out_filename.c_str());
    
    vector<string> poscar_in_line;
    split(poscar_in_string, '\n', poscar_in_line);
    
    
    vector<string> primitive_x_in_line;
    vector<string> primitive_y_in_line;
    vector<string> primitive_z_in_line;
    vector<double> primitive_x_double;
    vector<double> primitive_y_double;
    vector<double> primitive_z_double;
    
    //        cout<<poscar_in_line[2]<<endl<<poscar_in_line[3]<<endl<<poscar_in_line[4]<<endl;
    split(poscar_in_line[2], ' ', primitive_x_in_line);
    split(poscar_in_line[3], ' ', primitive_y_in_line);
    split(poscar_in_line[4], ' ', primitive_z_in_line);
    primitive_x_double.push_back(lexical_cast<double >(primitive_x_in_line[0])) ;
    primitive_x_double.push_back(lexical_cast<double >(primitive_x_in_line[1])) ;
    primitive_x_double.push_back(lexical_cast<double >(primitive_x_in_line[2])) ;
    primitive_y_double.push_back(lexical_cast<double >(primitive_y_in_line[0])) ;
    primitive_y_double.push_back(lexical_cast<double >(primitive_y_in_line[1])) ;
    primitive_y_double.push_back(lexical_cast<double >(primitive_y_in_line[2])) ;
    primitive_z_double.push_back(lexical_cast<double >(primitive_z_in_line[0])) ;
    primitive_z_double.push_back(lexical_cast<double >(primitive_z_in_line[1])) ;
    primitive_z_double.push_back(lexical_cast<double >(primitive_z_in_line[2])) ;
    
    //        printvector(primitive_x_double);
    //        printvector(primitive_y_double);
    //        printvector(primitive_z_double);
    
    vector<double> gs_primitive_x_double;
    vector<double> gs_primitive_y_double;
    vector<double> gs_primitive_z_double;
    
    gs_primitive_x_double.push_back(primitive_x_double[0]*periodicity.get<0>());
    gs_primitive_x_double.push_back(primitive_x_double[1]*periodicity.get<0>());
    gs_primitive_x_double.push_back(primitive_x_double[2]*periodicity.get<0>());
    gs_primitive_y_double.push_back(primitive_y_double[0]*periodicity.get<2>()+primitive_x_double[0]*periodicity.get<1>());
    gs_primitive_y_double.push_back(primitive_y_double[1]*periodicity.get<2>()+primitive_x_double[1]*periodicity.get<1>());
    gs_primitive_y_double.push_back(primitive_y_double[2]*periodicity.get<2>()+primitive_x_double[2]*periodicity.get<1>());
    gs_primitive_z_double.push_back(primitive_z_double[0]*periodicity.get<5>()+primitive_y_double[0]*periodicity.get<4>()+primitive_x_double[0]*periodicity.get<3>());
    gs_primitive_z_double.push_back(primitive_z_double[1]*periodicity.get<5>()+primitive_y_double[1]*periodicity.get<4>()+primitive_x_double[1]*periodicity.get<3>());
    gs_primitive_z_double.push_back(primitive_z_double[2]*periodicity.get<5>()+primitive_y_double[2]*periodicity.get<4>()+primitive_x_double[2]*periodicity.get<3>());
    
    
    vector<string> looping_numbers;
    
    split(poscar_in_line[5], ' ', looping_numbers);
    
    //        cout<<endl<<poscar_in_line[5];
    //        cout<<"\nindicator2\n";
    //
    //        printvector(looping_numbers);
    //        cout<<looping_numbers.size();
    //
    //        vector<string> line_of_a_term;
    //        split(poscar_in_line[8], ' ', line_of_a_term);
    //        cout<<"\n line numbers for line of a term is "<<line_of_a_term.size()  ;
    
    int total_lines=0;
    for (int i=0; i<looping_numbers.size(); i++) {
        total_lines+=lexical_cast<int>(looping_numbers[i]);
    }
    //        cout<<"\n total lines are "<<total_lines;
    
    
    vector< vector<double> > active_coordinates;
    vector< vector<string> > active_components;
    vector< vector<double> > passive_coordinates;
    vector< vector<string> > passive_components;
    
    for (int line_number=7 ; line_number<7+total_lines; line_number++) {
        int index=line_number-6;
        vector<string> read_line;
        split(poscar_in_line[line_number], ' ' , read_line);
        
        if (read_line.size()>4) {
            vector<double> vector_temp;
            vector_temp.push_back(lexical_cast<double>(read_line[0]));
            vector_temp.push_back(lexical_cast<double>(read_line[1]));
            vector_temp.push_back(lexical_cast<double>(read_line[2]));
            vector<string> vector_component;
            for (int temp_index=3; temp_index<read_line.size(); temp_index++) {
                vector_component.push_back(read_line[temp_index]);
            }
            active_coordinates.push_back(vector_temp);
            active_components.push_back(vector_component);
        }
        else if(read_line.size()==4) {
            
            vector<double> vector_temp;
            vector_temp.push_back(lexical_cast<double>(read_line[0]));
            vector_temp.push_back(lexical_cast<double>(read_line[1]));
            vector_temp.push_back(lexical_cast<double>(read_line[2]));
            vector<string> vector_component;
            for (int temp_index=3; temp_index<read_line.size(); temp_index++) {
                vector_component.push_back(read_line[temp_index]);
            }
            passive_coordinates.push_back(vector_temp);
            passive_components.push_back(vector_component);
            
        }
        
    }
    
    poscar_out << fixed << setprecision(8);
    //        cout << fixed << setprecision(8);
    poscar_out<<poscar_in_line[0]<<endl<<poscar_in_line[1]<<endl<<gs_primitive_x_double[0]<<" "<<gs_primitive_x_double[1]<<" "<<gs_primitive_x_double[2]<<endl<<gs_primitive_y_double[0]<<" "<<gs_primitive_y_double[1]<<" "<<gs_primitive_y_double[2]<<endl<<gs_primitive_z_double[0]<<" "<<gs_primitive_z_double[1]<<" "<<gs_primitive_z_double[2]<<endl;
    
    int a0,a1,a2,a3,a4,a5;
    a0=periodicity.get<0>();
    a1=periodicity.get<1>();
    a2=periodicity.get<2>();
    a3=periodicity.get<3>();
    a4=periodicity.get<4>();
    a5=periodicity.get<5>();
    
    vector<vector<double> > coordinate_in_poscar_out_before_transformation_sorting;
    vector<string> species_in_poscar_out_before_transformation_sorting;
    
    for (int i=1; i<=a0; i++) {
        for (int j=1; j<=a2; j++) {
            for (int k=1; k<=a5; k++) {
                for (int p=1; p<=components.size(); p++) {
                    int component_index=unitcell[make_tuple(i,j,k,p)];
                    vector<double> position_before_transformation;
                    position_before_transformation.push_back(i-1+active_coordinates[p-1][0]);
                    position_before_transformation.push_back(j-1+active_coordinates[p-1][1]);
                    position_before_transformation.push_back(k-1+active_coordinates[p-1][2]);
                    coordinate_in_poscar_out_before_transformation_sorting.push_back(position_before_transformation);
                    if (!tern_alg) {
                        species_in_poscar_out_before_transformation_sorting.push_back(active_components[p-1][components[p]-component_index-1]);
                    }
                    else if(tern_alg){
                        if (component_index==0)
                            species_in_poscar_out_before_transformation_sorting.push_back(active_components[p-1][1]);
                        
                        if (component_index==1)
                            species_in_poscar_out_before_transformation_sorting.push_back(active_components[p-1][0]);
                        
                        if (component_index==2)
                            species_in_poscar_out_before_transformation_sorting.push_back(active_components[p-1][2]);
                        
                    }
                }
                
                for (int passive_p=0; passive_p<passive_coordinates.size(); passive_p++) {
                    vector<double> position_before_transformation;
                    position_before_transformation.push_back(i-1+passive_coordinates[passive_p][0]);
                    position_before_transformation.push_back(j-1+passive_coordinates[passive_p][1]);
                    position_before_transformation.push_back(k-1+passive_coordinates[passive_p][2]);
                    coordinate_in_poscar_out_before_transformation_sorting.push_back(position_before_transformation);
                    species_in_poscar_out_before_transformation_sorting.push_back(passive_components[passive_p][0]);
                }
            }
        }
    }
    
    //        cout<<"\n try printing the thing";
    //        for (int temp=0; temp<coordinate_in_poscar_out_before_transformation_sorting.size(); temp++) {
    //            cout<<endl;
    //            for (int i=0; i<3; i++) {
    //                cout<< coordinate_in_poscar_out_before_transformation_sorting[temp][i];
    //                cout<< " ";
    //            }
    //            cout<<species_in_poscar_out_before_transformation_sorting[temp];
    //        }
    
    vector<vector<double> > coordinate_in_poscar_out_after_transformation_before_sorting;
    vector<string> species_in_poscar_out_after_transformation_before_sorting=species_in_poscar_out_before_transformation_sorting;
    
    
    //        cout<<"\n try converting the thing";
    for (int temp=0; temp<coordinate_in_poscar_out_before_transformation_sorting.size(); temp++) {
        double x,y,z;
        x=coordinate_in_poscar_out_before_transformation_sorting[temp][0];
        y=coordinate_in_poscar_out_before_transformation_sorting[temp][1];
        z=coordinate_in_poscar_out_before_transformation_sorting[temp][2];
        vector<double> temp_vec;
        temp_vec.push_back(x/a0-y/a2*a1/a0+z*a1*a4/a0/a2/a5-z/a5*a3/a0);
        temp_vec.push_back(y/a2-z/a5*a4/a2);
        temp_vec.push_back(z/a5);
        coordinate_in_poscar_out_after_transformation_before_sorting.push_back(temp_vec);
    }
    //
    //        cout<<"\n try printing the thing after transformation";
    //        for (int temp=0; temp<coordinate_in_poscar_out_after_transformation_before_sorting.size(); temp++) {
    //            cout<<endl;
    //            for (int i=0; i<3; i++) {
    //                cout<< coordinate_in_poscar_out_after_transformation_before_sorting[temp][i];
    //                cout<< " ";
    //            }
    //            cout<<species_in_poscar_out_after_transformation_before_sorting[temp];
    //        }
    
    map<string, vector<vector<double> > > reverse_species_coordinate;
    for (int temp=0; temp<coordinate_in_poscar_out_after_transformation_before_sorting.size(); temp++) {
        
        reverse_species_coordinate[species_in_poscar_out_after_transformation_before_sorting[temp]].push_back(coordinate_in_poscar_out_after_transformation_before_sorting[temp]);
    }
    
    
    vector<vector<double> > coordinate_in_poscar_out_after_transformation_after_sorting;
    vector<string> species_in_poscar_out_after_transformation_after_sorting=species_in_poscar_out_before_transformation_sorting;
    
    
    
    if (reverse_species_coordinate.count("Va")>0) {
        reverse_species_coordinate.erase("Va");
    }
    
    if (reverse_species_coordinate.count("Vac")>0) {
        reverse_species_coordinate.erase("Vac");
    }
    
    for (map<string, vector<vector<double> > >::iterator it1=reverse_species_coordinate.begin(); it1!=reverse_species_coordinate.end(); it1++) {
        poscar_out<<it1->first<<" ";
    }
    poscar_out<<endl;
    for (map<string, vector<vector<double> > >::iterator it1=reverse_species_coordinate.begin(); it1!=reverse_species_coordinate.end(); it1++) {
        poscar_out<<(it1->second).size()<<" ";
    }
    poscar_out<<endl;
    poscar_out<<poscar_in_line[6];
    poscar_out<<endl;
    for (map<string, vector<vector<double> > >::iterator it1=reverse_species_coordinate.begin(); it1!=reverse_species_coordinate.end(); it1++) {
        vector<vector<double> > coordinates_here=it1->second;
        for (int temp=0; temp<coordinates_here.size(); temp++) {
            vector<double> vector_here=coordinates_here[temp];
            poscar_out<<vector_here[0]-floor(vector_here[0])<<" "<<vector_here[1]-floor(vector_here[1])<<" "<<vector_here[2]-floor(vector_here[2])<<" "<<it1->first<<endl;
        }
    }
}

std::string get_file_contents(const char *filename){
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (in)
    {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize(in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();
        return(contents);
    }
    throw(errno);
}


spin_no_periodic_struct spin_no_periodic_struct::remove_x_neg()
{
    spin_no_periodic_struct spin_struct_temp;
    map<tuple<int,int,int,int>,int> spin_temp;
    for (map<tuple<int,int,int,int>,int>::iterator it1=spin.begin();it1!=spin.end() ;it1++)
    {
        tuple<int,int,int,int> tuple_here=it1->first;
        int type_here=it1->second;
        int x=tuple_here.get<0>();
        int y=tuple_here.get<1>();
        int z=tuple_here.get<2>();
        int p=tuple_here.get<3>();
        if (x>1.5) {
            spin_temp[make_tuple(x-1,y,z,p)]=type_here;
        }
    }
    spin_struct_temp.spin=spin_temp;
    return spin_struct_temp;
    
}

spin_no_periodic_struct spin_no_periodic_struct::remove_x_pos()
{
    spin_no_periodic_struct spin_struct_temp;
    map<tuple<int,int,int,int>,int> spin_temp;
    
    int range=0;
    for (map<tuple<int,int,int,int>,int>::iterator it1=spin.begin();it1!=spin.end() ;it1++)
    {
        tuple<int,int,int,int> tuple_here=it1->first;
        int type_here=it1->second;
        int x=tuple_here.get<0>();
        int y=tuple_here.get<1>();
        int z=tuple_here.get<2>();
        int p=tuple_here.get<3>();
        if (x>range) {
            range=x;
        }
    }
    
    for (map<tuple<int,int,int,int>,int>::iterator it1=spin.begin();it1!=spin.end() ;it1++)
    {
        tuple<int,int,int,int> tuple_here=it1->first;
        int type_here=it1->second;
        int x=tuple_here.get<0>();
        int y=tuple_here.get<1>();
        int z=tuple_here.get<2>();
        int p=tuple_here.get<3>();
        if (x!=range) {
            spin_temp[make_tuple(x,y,z,p)]=type_here;
        }
    }
    spin_struct_temp.spin=spin_temp;
    return spin_struct_temp;
    
}

spin_no_periodic_struct spin_no_periodic_struct::remove_y_neg()
{
    spin_no_periodic_struct spin_struct_temp;
    map<tuple<int,int,int,int>,int> spin_temp;
    for (map<tuple<int,int,int,int>,int>::iterator it1=spin.begin();it1!=spin.end() ;it1++)
    {
        tuple<int,int,int,int> tuple_here=it1->first;
        int type_here=it1->second;
        int x=tuple_here.get<0>();
        int y=tuple_here.get<1>();
        int z=tuple_here.get<2>();
        int p=tuple_here.get<3>();
        if (y>1.5) {
            spin_temp[make_tuple(x,y-1,z,p)]=type_here;
        }
    }
    spin_struct_temp.spin=spin_temp;
    return spin_struct_temp;
    
}

spin_no_periodic_struct spin_no_periodic_struct::remove_y_pos()
{
    spin_no_periodic_struct spin_struct_temp;
    map<tuple<int,int,int,int>,int> spin_temp;
    
    int range=0;
    for (map<tuple<int,int,int,int>,int>::iterator it1=spin.begin();it1!=spin.end() ;it1++)
    {
        tuple<int,int,int,int> tuple_here=it1->first;
        int type_here=it1->second;
        int x=tuple_here.get<0>();
        int y=tuple_here.get<1>();
        int z=tuple_here.get<2>();
        int p=tuple_here.get<3>();
        if (y>range) {
            range=y;
        }
    }
    
    for (map<tuple<int,int,int,int>,int>::iterator it1=spin.begin();it1!=spin.end() ;it1++)
    {
        tuple<int,int,int,int> tuple_here=it1->first;
        int type_here=it1->second;
        int x=tuple_here.get<0>();
        int y=tuple_here.get<1>();
        int z=tuple_here.get<2>();
        int p=tuple_here.get<3>();
        if (y!=range) {
            spin_temp[make_tuple(x,y,z,p)]=type_here;
        }
    }
    spin_struct_temp.spin=spin_temp;
    return spin_struct_temp;
    
}


spin_no_periodic_struct spin_no_periodic_struct::remove_z_neg()
{
    spin_no_periodic_struct spin_struct_temp;
    map<tuple<int,int,int,int>,int> spin_temp;
    for (map<tuple<int,int,int,int>,int>::iterator it1=spin.begin();it1!=spin.end() ;it1++)
    {
        tuple<int,int,int,int> tuple_here=it1->first;
        int type_here=it1->second;
        int x=tuple_here.get<0>();
        int y=tuple_here.get<1>();
        int z=tuple_here.get<2>();
        int p=tuple_here.get<3>();
        if (z>1.5) {
            spin_temp[make_tuple(x,y,z-1,p)]=type_here;
        }
    }
    spin_struct_temp.spin=spin_temp;
    return spin_struct_temp;
    
}

spin_no_periodic_struct spin_no_periodic_struct::remove_z_pos()
{
    spin_no_periodic_struct spin_struct_temp;
    map<tuple<int,int,int,int>,int> spin_temp;
    
    int range=0;
    for (map<tuple<int,int,int,int>,int>::iterator it1=spin.begin();it1!=spin.end() ;it1++)
    {
        tuple<int,int,int,int> tuple_here=it1->first;
        int type_here=it1->second;
        int x=tuple_here.get<0>();
        int y=tuple_here.get<1>();
        int z=tuple_here.get<2>();
        int p=tuple_here.get<3>();
        if (z>range) {
            range=z;
        }
    }
    
    for (map<tuple<int,int,int,int>,int>::iterator it1=spin.begin();it1!=spin.end() ;it1++)
    {
        tuple<int,int,int,int> tuple_here=it1->first;
        int type_here=it1->second;
        int x=tuple_here.get<0>();
        int y=tuple_here.get<1>();
        int z=tuple_here.get<2>();
        int p=tuple_here.get<3>();
        if (z!=range) {
            spin_temp[make_tuple(x,y,z,p)]=type_here;
        }
    }
    spin_struct_temp.spin=spin_temp;
    return spin_struct_temp;
    
}


void calculate_cluster_type_and_energy_no_periodic(map<set<tuple<int,int,int,int,int> >, double> J_for_proof,double &energy,map<set<tuple<int,int,int,int,int> >, double> &cluster_type_here,map<tuple<int,int,int,int>,int>spin_now)
{
    cluster_type_here.clear();
    
    for (map<set<tuple<int,int,int,int,int> >, double>::iterator it1= J_for_proof.begin(); it1!=J_for_proof.end(); it1++) {
        set<tuple<int,int,int,int,int> > thisset = it1->first;
        cluster_type_here[thisset]=0;
        int i=1;
        int j=1;
        int k=1;
        {
            int all_matched_indicator = 0;
            for (set<tuple<int,int,int,int,int> >::iterator it1 = thisset.begin();
                 it1 != thisset.end(); it1++) {
                tuple<int,int,int,int,int> temp_tuple = *it1;
                int x = temp_tuple.get<0>();
                int y = temp_tuple.get<1>();
                int z = temp_tuple.get<2>();
                int position = temp_tuple.get<3>();
                int var = temp_tuple.get<4>();
                if (spin_now[make_tuple(i + x - 1, j + y - 1, k + z - 1, position)] != var) {
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
                cluster_type_here[thisset] += 1.0/(1 * 1 * 1);
            }
        }
    }
    
    energy = 0;
    
    for (map<set<tuple<int,int,int,int,int> >, double>::iterator it = cluster_type_here.begin();
         it != cluster_type_here.end(); it++) {
        energy += J_for_proof[it->first]*it->second;
    }
    
};


void calculate_cluster_type_and_energy_periodic(map<set<tuple<int,int,int,int,int> >, double> J_for_proof,tuple<int,int,int,int,int,int>periodicity ,double &energy,map<set<tuple<int,int,int,int,int> >, double> &cluster_type_here,map<tuple<int,int,int,int>,int>spin_now)
{
    cluster_type_here.clear();
    
    int a0=periodicity.get<0>();
    int a1=periodicity.get<1>();
    int a2=periodicity.get<2>();
    int a3=periodicity.get<3>();
    int a4=periodicity.get<4>();
    int a5=periodicity.get<5>();
    map<int,int> componentnumber;
    int x_range, y_range, z_range;
    calculate_range_from_J(J_for_proof, x_range, y_range, z_range, componentnumber);

    spin_now=extend_periodic_spin(periodicity, spin_now, J_for_proof);
    
    
    
    for (map<set<tuple<int,int,int,int,int> >, double>::iterator it = J_for_proof.begin();
         it != J_for_proof.end(); it++){
        set<tuple<int,int,int,int,int> > thisset = it->first;
        cluster_type_here[thisset] = 0;
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
                        assert(spin_now.count(make_tuple(i + x - 1, j + y - 1, k + z - 1, position))>0);
                        if (spin_now[make_tuple(i + x - 1, j + y - 1, k + z - 1, position)] != var) {
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
                        cluster_type_here[thisset] += 1.0/(a0*a2*a5);
                    }
                }
            }
        }
    }

    
    energy = 0;
    
    for (map<set<tuple<int,int,int,int,int> >, double>::iterator it = cluster_type_here.begin();
         it != cluster_type_here.end(); it++) {
        energy += J_for_proof[it->first]*it->second;
    }
    
}


map<tuple<int,int,int,int>,int> extend_periodic_spin(tuple<int,int,int,int,int,int>periodicity,map<tuple<int,int,int,int>,int>spin_now ,map<set<tuple<int,int,int,int,int> >, double> J)
{
    int a0=periodicity.get<0>();
    int a1=periodicity.get<1>();
    int a2=periodicity.get<2>();
    int a3=periodicity.get<3>();
    int a4=periodicity.get<4>();
    int a5=periodicity.get<5>();
    map<int,int> componentnumber;
    int x_range, y_range, z_range;
    calculate_range_from_J(J, x_range, y_range, z_range, componentnumber);
    map<tuple<int,int,int,int>,int> output_spin;
    for (int i = 1; i < a0 + x_range; i++) {
        for (int j = 1; j < a2 + y_range; j++) {
            for (int k = 1; k < a5 + z_range; k++) {
                for (map<int, int>::iterator cit = componentnumber.begin();
                     cit != componentnumber.end(); cit++) {
                    
                    int position = cit->first;
                    int i0 = positive_modulo(((i-1)-floor_int_division(j-(k-1)/a5*a4,a2)*a1-(k-1)/a5*a3), a0) + 1;
                    int j0 = positive_modulo(((j-1)-(k-1)/a5*a4), a2) + 1;
                    int k0 = ((k-1) % a5) + 1;
                    
                    output_spin[make_tuple(i,j,k,position)]=spin_now[make_tuple(i0,j0,k0,position)];

                }
            }
        }
    }
    return output_spin;
    
}

void construct_set_of_equivalent_from_prototype_set(set<set<tuple<int,int,int,int,int> > > &set_of_equivalent,set<tuple<int,int,int,int,int> > &prototype_set,int x_range,int y_range,int z_range)
{
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
}


void calculate_range_from_J(map<set<tuple<int,int,int,int,int> >, double> &J,
                            int &x_range,int &y_range,int &z_range,
                            map<int, int>&component)
{
    int xmax=0;
    int ymax=0;
    int zmax=0;
    for (map<set<tuple<int,int,int,int,int> >, double>::iterator it=J.begin(); it!=J.end(); it++) {
        set<tuple<int,int,int,int,int> > tuple_set_temp=it->first;
        for (set<tuple<int,int,int,int,int> >::iterator itv=tuple_set_temp.begin(); itv!=tuple_set_temp.end(); itv++) {
            if (xmax<(*itv).get<0>()) {
                xmax=(*itv).get<0>();
            }
            if (ymax<(*itv).get<1>()) {
                ymax=(*itv).get<1>();
            }
            if (zmax<(*itv).get<2>()) {
                zmax=(*itv).get<2>();
            }
            
            int location=(*itv).get<3>();
            int elementtype=(*itv).get<4>();
            
            if (component.count(location)==0) {
                component[location]=elementtype+1;
            }
            else if (component[location]<elementtype+1)
            {
                component[location]=elementtype+1;
            }
        }
    }
    x_range=xmax;
    y_range=ymax;
    z_range=zmax;
}


void  convert_J_input_tern_unaltered_to_J_tern_in(map<set< tuple<int,int,int,int,int,int> >,double > J_input_tern_unaltered,map<set< tuple<int,int,int,int,int> >,double >& J_tern_in)
{
    int all_terms_converted=0;
    
    while (all_terms_converted==0) {

        all_terms_converted=1;
        set< tuple<int,int,int,int,int,int> > temp_set_tuple_unaltered;
        double value;
        for (map<set< tuple<int,int,int,int,int,int> >,double >::iterator it1=J_input_tern_unaltered.begin() ; it1!=J_input_tern_unaltered.end(); it1++) {
            temp_set_tuple_unaltered=it1->first;
            value=it1->second;
            for (set< tuple<int,int,int,int,int,int> >::iterator it2=temp_set_tuple_unaltered.begin(); it2!=temp_set_tuple_unaltered.end(); it2++) {
                tuple<int,int,int,int,int,int> tuple_unaltered_now=*it2;
                int i0=tuple_unaltered_now.get<0>();
                int i1=tuple_unaltered_now.get<1>();
                int i2=tuple_unaltered_now.get<2>();
                int i3=tuple_unaltered_now.get<3>();
                int i4=tuple_unaltered_now.get<4>();
                int i5=tuple_unaltered_now.get<5>();
                
                if(i5!=-1)
                    
                {
                    all_terms_converted=0;
                    goto PORTAL_1324;
                }
                
            }
        }
    PORTAL_1324:
        
        if(all_terms_converted==0){
            for (set< tuple<int,int,int,int,int,int> >::iterator it2=temp_set_tuple_unaltered.begin(); it2!=temp_set_tuple_unaltered.end(); it2++) {
                tuple<int,int,int,int,int,int> tuple_unaltered_now=*it2;
                int i0=tuple_unaltered_now.get<0>();
                int i1=tuple_unaltered_now.get<1>();
                int i2=tuple_unaltered_now.get<2>();
                int i3=tuple_unaltered_now.get<3>();
                int i4=tuple_unaltered_now.get<4>();
                int i5=tuple_unaltered_now.get<5>();
                
                if (i5!=-1) {
                    set< tuple<int,int,int,int,int,int> > temp_set_tuple_unaltered_copy_1=temp_set_tuple_unaltered;
                    
                    set< tuple<int,int,int,int,int,int> > temp_set_tuple_unaltered_copy_2=temp_set_tuple_unaltered;
                    
                    if (i5==0) {
                        temp_set_tuple_unaltered_copy_1.erase(tuple_unaltered_now);
                        temp_set_tuple_unaltered_copy_1.insert(make_tuple(i0,i1,i2,i3,1,-1));
                        
                        temp_set_tuple_unaltered_copy_2.erase(tuple_unaltered_now);
                        temp_set_tuple_unaltered_copy_2.insert(make_tuple(i0,i1,i2,i3,2,-1));
                        
                        J_input_tern_unaltered.erase(temp_set_tuple_unaltered);
                        
                        J_input_tern_unaltered[temp_set_tuple_unaltered_copy_1]+=value;
                        J_input_tern_unaltered[temp_set_tuple_unaltered_copy_2]-=value;
                        
                        
                    }
                    if (i5==1) {
                        temp_set_tuple_unaltered_copy_1.erase(tuple_unaltered_now);
                        temp_set_tuple_unaltered_copy_1.insert(make_tuple(i0,i1,i2,i3,1,-1));
                        
                        temp_set_tuple_unaltered_copy_2.erase(tuple_unaltered_now);
                        temp_set_tuple_unaltered_copy_2.insert(make_tuple(i0,i1,i2,i3,2,-1));
                        
                        J_input_tern_unaltered.erase(temp_set_tuple_unaltered);
                        
                        J_input_tern_unaltered[temp_set_tuple_unaltered_copy_1]+=value;
                        J_input_tern_unaltered[temp_set_tuple_unaltered_copy_2]+=value;
                        
                        
                    }
                    
                    break;
                    
                }
                
                
            }
        }
        

        
        
        
    }
    
    for (map<set< tuple<int,int,int,int,int,int> >,double >::iterator it1=J_input_tern_unaltered.begin() ; it1!=J_input_tern_unaltered.end(); it1++) {
        set< tuple<int,int,int,int,int,int> > temp_set_tuple_unaltered_after_convert=it1->first;
        double value_after_convert=it1->second;
        
        
        set< tuple<int,int,int,int,int> > temp_set_tuple_after_convert;
        for (set< tuple<int,int,int,int,int,int> >::iterator it2=temp_set_tuple_unaltered_after_convert.begin(); it2!=temp_set_tuple_unaltered_after_convert.end(); it2++) {
            tuple<int,int,int,int,int,int> tuple_unaltered_now=*it2;
            int i0=tuple_unaltered_now.get<0>();
            int i1=tuple_unaltered_now.get<1>();
            int i2=tuple_unaltered_now.get<2>();
            int i3=tuple_unaltered_now.get<3>();
            int i4=tuple_unaltered_now.get<4>();
            int i5=tuple_unaltered_now.get<5>();
            
            assert(i5==-1);
            
            temp_set_tuple_after_convert.insert(make_tuple(i0,i1,i2,i3,i4));
            
        }
        
        J_tern_in[temp_set_tuple_after_convert]+=value_after_convert;
        
    }
    
    
    
    

}

void get_concentration_formation_energy_from_spin_and_J_in_ternary(tuple<int,int,int,int,int,int> periodicity_now,map< set<tuple<int,int,int,int,int> >, double > mu_first, map< set<tuple<int,int,int,int,int> >, double > mu_second ,map< set<tuple<int,int,int,int,int> >, double > J_fixed_part, map<tuple<int,int,int,int>,int>  spin_now, double constant, double mu_constant,double&formation_energy_out, tuple<double,double> &concentration_out,map< set<tuple<int,int,int,int,int> >, double > &clustertype_out)
{
    map<tuple<int,int,int,int>,int> new_spin_now=extend_periodic_spin(periodicity_now, spin_now, J_fixed_part);
    spin_now=new_spin_now;
    
    
    map< set<tuple<int,int,int,int,int> >, double > clustertype;
    
    int a0=periodicity_now.get<0>();
    int a1=periodicity_now.get<1>();
    int a2=periodicity_now.get<2>();
    int a3=periodicity_now.get<3>();
    int a4=periodicity_now.get<4>();
    int a5=periodicity_now.get<5>();
    
    for (map< set<tuple<int,int,int,int,int> >, double >::iterator it = J_fixed_part.begin();
         it != J_fixed_part.end(); it++){
        set<tuple<int,int,int,int,int> > thisset = it->first;
        clustertype[thisset] = 0;
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
                        if (spin_now[make_tuple(i + x - 1, j + y - 1, k + z - 1, position)] != var) {
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
                        clustertype[thisset] += 1.0/(a0*a2*a5);
                    }
                }
            }
        }
    }
    
    
    double formation_energy=0;
    for (map< set<tuple<int,int,int,int,int> >, double >::iterator it = J_fixed_part.begin();it != J_fixed_part.end(); it++) {
        formation_energy+=it->second*clustertype[it->first];
    }
    formation_energy+=constant;
    formation_energy-=mu_constant;
    
    //    cout<<"\ndebug here what is formation energy: "<<setw(20) << std::fixed << std::setprecision(8)<<formation_energy;
    
    double first_concentration_for_this_periodicity=0;
    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu_first.begin(); it1!=mu_first.end(); it1++) {
        if (clustertype.count(it1->first)==1) {
            first_concentration_for_this_periodicity+=clustertype[it1->first];
        }
        else{
            perror("error please check, error id 52919sdi");
            return ;
        }
    }
    first_concentration_for_this_periodicity=first_concentration_for_this_periodicity/mu_first.size();
    
    double second_concentration_for_this_periodicity=0;
    for (map< set<tuple<int,int,int,int,int> >, double>::iterator it1=mu_second.begin(); it1!=mu_second.end(); it1++) {
        if (clustertype.count(it1->first)==1) {
            second_concentration_for_this_periodicity+=clustertype[it1->first];
        }
        else{
            perror("error please check, error id 52919sdifqws");
            return ;
        }
    }
    second_concentration_for_this_periodicity=second_concentration_for_this_periodicity/mu_second.size();
    
    
    formation_energy_out=formation_energy;
    concentration_out=make_tuple(first_concentration_for_this_periodicity,second_concentration_for_this_periodicity);
    //    cout<<" concentration_for_this_periodicity is "<<concentration_for_this_periodicity;
    
    clustertype_out=clustertype;
}

int floor_int_division(int up,int down)
{
    double tmp=(double)up/down+1e-10;
    return floor(tmp);
}


