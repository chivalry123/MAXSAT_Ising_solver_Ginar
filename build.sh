#!/bin/bash

################################# SETUP #######################################
# Change these to the proper Gurobi INCLUDE and LIB directories
export ISING_GUROBI_INC=/home/key01027/ising_solver/gurobi604/linux64/include
export ISING_GUROBI_LIB=/home/key01027/ising_solver/gurobi604/linux64/lib
export BOOST_INC=/home/key01027/ising_solver/boost_install/include
export BOOST_LIB=/home/key01027/ising_solver/boost_install/lib

###############################################################################
echo "------------------------------------------------------"
echo "Cleaning old bins..."
test -d bin || mkdir -p bin
rm -f bin/ising
rm -f bin/akmaxsat*
rm -f bin/CCLS*

echo "------------------------------------------------------"
echo "Building Ising solver..."
cd src/ising
make clean
make
cp bin/ising ../../bin/
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

echo "---------------------------------------------------"
echo "build CCEHC.."
cd src/CCEHC/
make cleanup
make
cp CCEHC ../../bin/
cd ../../../



