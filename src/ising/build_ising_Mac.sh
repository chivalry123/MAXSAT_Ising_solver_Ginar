rm ising


mpic++ main.cpp tools.cpp treeofdevil.cpp solver.cpp periodicfunction.cpp  -o ising -lgurobi_c++ -lgurobi65  -lpthread -lm -Wall -m64 -fPIE -fexceptions -frounding-math -O2 -I/Library/gurobi605/mac64/include/ -I/usr/local/include/ -I. -L/Library/gurobi605/mac64/lib/ -L/usr/local/lib/ -lboost_regex-mt -lboost_system-mt -lboost_mpi-mt -lboost_serialization-mt

cp ising 97_arena/
cp ising binary/
cp ising debug_mode/
cp ising binary_incomplete/
cp ising 97_arena_tern_incomplete/
cp ising binary_incomplete_major_periodicity/
cp ising binary_largecell_multi_mu_incomplete/
cp ising binary_20160801_constrain_x_concentration/
cp ising binary_incomplete_20160824_simulated_annealing/
cp ising binary_multiple_sites_20160909_Tina/

cp ising single_point_solver/

