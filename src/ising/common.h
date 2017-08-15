//
//  common.h
//  ising 103 MPI conversion
//
//  Created by Wenxuan Huang on 10/27/15.
//  Copyright (c) 2015 wenxuan huang. All rights reserved.
//

#ifndef ising_103_MPI_conversion_common_h
#define ising_103_MPI_conversion_common_h


#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <string>
#include <stdio.h>
#include <vector>
#include <set>
#include <map>
#include "gurobi_c++.h"
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"
#include <boost/lexical_cast.hpp>
#include "solver.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <boost/thread.hpp>
//#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
//#include <boost/serialization/tuple.hpp>
#include <boost/serialization/utility.hpp>
//#include <boost/serialization/nvp.hpp>
//#include <boost/preprocessor/repetition.hpp>
#include <boost/mpi.hpp>
#include "boost/tuple/tuple.hpp"
#include <boost/serialization/string.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include "boost/tuple/tuple_io.hpp"
#include <math.h>       /* floor */

using namespace ::boost::tuples;
using namespace ::boost;

//extern boost::mpi::environment mpi_env;
//extern boost::mpi::communicator mpi_world;




namespace boost { namespace serialization {

    
    template<typename Archive, typename T1>
    
    void serialize(Archive & ar,
                   boost::tuple<T1> & t,
                   
                   const unsigned int)
    {
        ar & t.get<0>();
    }
    
    template<typename Archive, typename T1,typename T2>
    void serialize(Archive & ar,
                   
                   boost::tuple<T1,T2> & t,
                   
                   const unsigned int)
    {
        
        ar & t.get<0>();
        ar & t.get<1>();
    }
    
    
    template<typename Archive, typename T1,typename T2,typename T3>
    void serialize(Archive & ar,
                   
                   boost::tuple<T1, T2, T3> & t,
                   
                   const unsigned int)
    {
        
        ar & t.get<0>();
        ar & t.get<1>();
        ar & t.get<2>();
    }
    
    
    template<typename Archive, typename T1,typename T2,typename T3 ,typename T4>
    void serialize(Archive & ar,
                   
                   boost::tuple<T1, T2, T3 ,T4> & t,
                   
                   const unsigned int)
    {
        
        ar & t.get<0>();
        ar & t.get<1>();
        ar & t.get<2>();
        ar & t.get<3>();
    }
    
    template<typename Archive, typename T1,typename T2,typename T3 ,typename T4,typename T5>
    void serialize(Archive & ar,
                   
                   boost::tuple<T1, T2, T3 ,T4, T5> & t,
                   
                   const unsigned int)
    {
        
        ar & t.get<0>();
        ar & t.get<1>();
        ar & t.get<2>();
        ar & t.get<3>();
        ar & t.get<4>();
    }
    
    template<typename Archive, typename T1,typename T2,typename T3 ,typename T4,typename T5,typename T6 >
    void serialize(Archive & ar,
                   
                   boost::tuple<T1, T2, T3 ,T4, T5, T6> & t,
                   
                   const unsigned int)
    {
        
        ar & t.get<0>();
        ar & t.get<1>();
        ar & t.get<2>();
        ar & t.get<3>();
        ar & t.get<4>();
        ar & t.get<5>();
    }
    
    
    template<typename Archive, typename T1,typename T2,typename T3 ,typename T4,typename T5,typename T6 ,typename T7>
    void serialize(Archive & ar,
                   
                   boost::tuple<T1, T2, T3 ,T4, T5, T6 ,T7> & t,
                   
                   const unsigned int)
    {
        
        ar & t.get<0>();
        ar & t.get<1>();
        ar & t.get<2>();
        ar & t.get<3>();
        ar & t.get<4>();
        ar & t.get<5>();
        ar & t.get<6>();
    }

    
}}



#endif
