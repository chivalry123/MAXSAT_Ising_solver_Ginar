WORKDIR = `pwd`
ISING_GUROBI_INC=~/ising_solver/gurobi650/linux64/include/
ISING_GUROBI_LIB=~/ising_solver/gurobi650/linux64/lib/
BOOST_INC=~/ising_solver/boost_install/include/
BOOST_LIB=~/ising_solver/boost_install/lib/

CC = gcc
CXX = mpiicpc
AR = ar
LD = mpiicpc

INC = -I $(ISING_GUROBI_INC) -I $(BOOST_INC) -I.
CFLAGS = -Wall -m64 -fPIE -fexceptions -frounding-math -O2 -fopenmp
LIB = -L $(ISING_GUROBI_LIB) -L $(BOOST_LIB) -parallel -lgurobi_c++ -lgurobi65 -lboost_regex -lboost_thread -lboost_system  -lboost_mpi -lboost_serialization

OBJDIR = obj
OUT_FINAL = bin/ising

OBJ = $(OBJDIR)/main.o $(OBJDIR)/solver.o $(OBJDIR)/periodicfunction.o $(OBJDIR)/tools.o $(OBJDIR)/treeofdevil.o

all: prepare build

prepare:
	test -d bin || mkdir -p bin
	test -d $(OBJDIR) || mkdir -p $(OBJDIR)

build: $(OBJ)
	$(LD) -o $(OUT_FINAL) $(OBJ) $(LDFLAGS) $(LIB)

$(OBJDIR)/main.o: main.cpp
	$(CXX) $(CFLAGS) $(INC) -c main.cpp -o $(OBJDIR)/main.o

$(OBJDIR)/solver.o: solver.cpp
	$(CXX) $(CFLAGS) $(INC) -c solver.cpp -o $(OBJDIR)/solver.o

$(OBJDIR)/periodicfunction.o: periodicfunction.cpp
	$(CXX) $(CFLAGS) $(INC) -c periodicfunction.cpp -o $(OBJDIR)/periodicfunction.o

$(OBJDIR)/tools.o: tools.cpp
	$(CXX) $(CFLAGS) $(INC) -c tools.cpp -o $(OBJDIR)/tools.o

$(OBJDIR)/treeofdevil.o: treeofdevil.cpp
	$(CXX) $(CFLAGS) $(INC) -c treeofdevil.cpp -o $(OBJDIR)/treeofdevil.o

clean:
	rm -f $(OBJDIR)/*
	rm -f $(OUT_FINAL)
