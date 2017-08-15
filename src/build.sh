#!/bin/bash

################################# SETUP #######################################
# Change these to the proper Gurobi INCLUDE and LIB directories
export ISING_GUROBI_INC=/home1/03631/key01027/ising_solver/gurobi604/linux64/include/
export ISING_GUROBI_LIB=/home1/03631/key01027/ising_solver/gurobi604/linux64/lib/
export BOOST_INC=/opt/apps/intel13/boost/1.55.0/x86_64/include/
export BOOST_LIB=/opt/apps/intel13/boost/1.55.0/x86_64/lib/
export BOOST_MPI_INC=/opt/apps/intel13/mvapich2_1_9/boost-mpi/1.55.0/include/
export BOOST_MPI_LIB=/opt/apps/intel13/mvapich2_1_9/boost-mpi/1.55.0/lib/
###############################################################################
echo "------------------------------------------------------"
echo "Cleaning old bins..."
test -d bin || mkdir -p bin
rm -f bin/ising
rm -f bin/akmaxsat*
rm -f bin/CCLS*
rm -f bin/polytope

echo "------------------------------------------------------"
echo "Building Ising solver..."
cd src/ising
make clean
make
cp bin/ising ../../bin/
cd ../../

echo "------------------------------------------------------"
echo "Building Polytope code..."
cd src/polytope
make clean
make
cp bin/polytope ../../bin/
cd ../../

echo "------------------------------------------------------"
echo "Building MAX-SAT solvers..."
cd src/ccls/akmaxsat_1.1_src
make clean
make
cp akmaxsat ../../../bin/
cd ../../../

cd src/ccls/ccls_2014_src
make cleanup
make
cp CCLS2014 ../../../bin/
cd ../../../

cd src/ccls/ccls_to_akmaxsat_src
make cleanup
make
cp CCLS_to_akmaxsat ../../../bin/
cd ../../../

cd src/ccls_lb/akmaxsat_1.1_src
make clean
make
cp akmaxsat_LB ../../../bin/
cd ../../../

cd src/ccls_lb/ccls_2014_src
make cleanup
make
cp CCLS2014_LB ../../../bin/
cd ../../../

cd src/ccls_lb/ccls_to_akmaxsat_src
make cleanup
make
cp CCLS_to_akmaxsat_LB ../../../bin/
cd ../../../
