#MAXSAT Ising solver Ginar

## How to get everything start working on Ginar

This github repo represents the folder `/home/key01027/ising_solver/ising_119_some_debug/` in Ginar (which is already shared to all of you).

All the executables are from Wenxuan's directory in Ginar:

 `/home/key01027/ising_solver/ising_119_some_debug/bin/`

You will need to copy the binaries to the destination location with the corresponding input file to make it work.

## Examples 

In the following, I have preset the J_config.in file for different applications. Except the NSITES which simply means how large the supercell size do you want to go into searching for optimal solution. You don't need to set anything else.

### example 1, solution to fixed ECIs

**The single point solver**: Please look into: `examples/single_point_solver`. Given an ECIs in its 0,1 formulation , solve for the optimal solution up to super cell size of NSITES

To do this, you would need 3 files: `PRIM, J_config.in, J_in.in` 

`PRIM` indicates the lattice system in the CASM format. `J_config.in` includes options for the solver including NSITES, please also Check `README` for more inforamtion. If you are too busy for that, simply use the current settings. `J_in.in ` dictates the ECIs, in its 0/1 formulation. 



After requesting a compute node (with qrsh, or qsub), you could use to run the solver

	mpiexec.hydra -np 16 ./ising 

The results would be

	[key01027@compute-1-83 single_point_solver]$ mpiexec.hydra -np 8 ./ising 
	[MASTER] Sent all jobs.
	[MASTER] Handled all jobs, killed every process.
	
	 I am still alive at point: sjoqno2hcb02 
	
	printing POSCAROUT
	
	 tell me what is upper bound: -5.957034540000005e-01 what is lower bound: -1.000000000000000e+18

You could check the `POSCAROUT` for the predicted Ground state up to that supercell size. (For lower bound to run successfully, which we generally do not use in actual practice, we need to set `NLOOPS`>0 and have the gurobi server correctly run, which needs to be renewed yearly.)

### example 2, finding ground state phase diagram for binary system

Please look into: `examples/binary_system`. Perfom binary GS phase diagram search based on CASM output. 

CASM files that are needed: `eci.out  FCLUST  PRIM` . copy the `useful_scripts/casm2wenxuan.py` to the running directory (`examples/binary_system`). Run `python casm2wenxuan.py > J_in.in` to create `J_in.in`

Additional GS file needed include `J_config.in`, `mu_in.in`, `mu_config.in`. When you have a different system, you need to change mu_in.in manually. For example, when your PRIM is, 

	Mg2 Cr4 O8
	1.0
	4.163840 4.163840 0.000000
	4.163840 0.000000 -4.163840
	0.000000 4.163840 -4.163840
	2 4 8
	direct
	0.375000 0.875000 0.875000 Mg Va
	0.625000 0.125000 0.125000 Mg Va
	0.500000 0.500000 0.500000 Cr
	0.000000 0.500000 0.500000 Cr
	0.000000 0.500000 0.000000 Cr
	0.000000 0.000000 0.500000 Cr
	0.776100 0.741300 0.741300 O
	0.223900 0.258700 0.258700 O
	0.241300 0.741300 0.741300 O
	0.758700 0.258700 0.258700 O
	0.241300 0.741300 0.276100 O
	0.758700 0.258700 0.723900 O
	0.241300 0.276100 0.741300 O
	0.758700 0.723900 0.258700 O

you need to specified the point term about what x is in `mu_in.in`. Specifically, it would be 

	Constant 0
	
	Cluster 1
	1,1,1,1,1   
	J=0
	
	Cluster 2
	1,1,1,2,1   
	J=0

The two point term `1,1,1,1,1` and `1,1,1,2,1` belongs to the same concentration group. (note that the first `1,1,1` corresponds to the x,y,z, the forth term `1` and `2` corresponds to the sub-lattice sites of the two Mg, Va species. The fifth term, `1` corresponds to the first element, in this case it is the Mg).

`mu_config.in` only includes only one line `MODE_JPLUSMINUS=1`

You can then run GS solver with  

	mpiexec.hydra -np 16 ./ising

you will have the folder called `GS_solutions` containing the binary ground state phase diagram. The difference between `POSCAR_OUT` and `POSCAR` in `GS_solutions` is that `POSCAR_OUT` contains the explicit vacancy terms `Va` or `Vac`, whereas `POSCAR` does not. The direcotry tree of `GS_solutions` looks like:

	$tree GS_solutions/
	GS_solutions/
	|-- hull_debug.txt
	|-- hull.txt
	|-- x0
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.250000
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.500000
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	`-- x1
	    |-- POSCAR
	    `-- POSCAR_OUT


**more comments below if you want to edit J_config.in or** `mu_config.in`  
The difference compared to pointwise solver is ``MODE_JPLUSMINUS=1`` meaning that the J input from CASM is +1/-1 formulation. `work_with_mu=1` and `scan_chemical_potential=1` to ensure scanning (we are pin-pointing to be specific) of chemical potential for ground states.

The flexible flags include `binary_x_min=0`,`binary_x_max=1`,`NSITES=2`. `binary_x_min` and `binary_x_max` is used to indicate the minimum and maximum interested concentration. `NSITES` denotes the super cell size to search up to.

In terms of what is x, for example, in the PRIM file, `0.375000 0.875000 0.875000 Mg Va`. Since Mg corresponds to +1 and Va corresponds to -1, concentration fo Mg is x. 


### example 3 finding ground state phase diagram for ternary system

Please look into: `examples/ternary_system`.

CASM files that are needed: `eci.out  FBCLUST  PRIM` . Noted that we would need `FBCLUST` this time instead of `FCLUST` for practical reasons. copy the `useful_scripts/casm2wenxuan_ternary.py` to the running directory (`examples/binary_system`). Run `python casm2wenxuan_ternary.py > J_in_tern_casm.in` to create `J_in_tern_casm.in `

If you are interested in the internal detail of `J_in.in`, you would see that there is one more index than the previous five. This is due to some very practical issues CASM use 1/0/-1 CE system and the spin variables needs to have square term. So if the sixth index is 1, it correspodns to square term `s_i^2`, if the sixth index is 0 it corresponds to linear term `s_i`.

The PRIM looks like (note that we should place all the substitutions terms at front, otherwise it causes a lot more computations due to the current implementation...):

	$ cat PRIM 
	M2 Ni1 O1
	4.104210848570817
	0.0 0.5 0.5
	0.5 0.0 0.5
	0.5 0.5 0.0
	1  2  1 
	direct
	0.000000 0.000000 0.000000 Li Ni Vac 
	0.250000 0.250000 0.250000 Li Ni Vac
	0.750000 0.750000 0.750000 Li Ni Vac
	0.500000 0.500000 0.500000 O

Additional GS file needed include `J_config.in`, `mu_in_tern_casm.in`, `mu_config.in`. The `mu_in_tern_casm.in` needs to be constructed manually, in our case, this system is 

	$ cat mu_in_tern_casm.in 	
	Constant 0
	
	Cluster 1
	Group 1
	1,1,1,1,1   
	J=0
	
	Cluster 2
	Group 1
	1,1,1,2,1   
	J=0
	
	Cluster 3
	Group 1
	1,1,1,3,1   
	J=0
	
	Cluster 4
	Group 2
	1,1,1,1,2   
	J=0
	
	Cluster 5
	Group 2
	1,1,1,2,2   
	J=0
	
	Cluster 6
	Group 2
	1,1,1,3,2   
	J=0

Group 1 and Group 2 is used to denote the two types of species correspondingly. To enable ternary algorithm, you need to set in `J_config.in`

	work_with_mu=0
	scan_chemical_potential=0
	ternary_alg=1

The current solver is a bit stupid that if you set `work_with_mu=1` and `scan_chemical_potential=1`. It just activates the binary search algorithm even though you set `ternary_alg=1`. So, to run everything successfully you need to set `work_with_mu=0`, `scan_chemical_potential=0`.


Since in CASM for the specific PRIM Li is +1, Vac is -1, the GS code associate +1 with x and -1 with y autoamtically. If you know the physically interested range of concentration, you could set it here. This could save tons of computation time needed. Since the chemical potential pinpointing is very intense in ternary system:

	ternary_x_min=0
	ternary_x_max=0.17
	ternary_y_min=0
	ternary_y_max=1
	ternary_z_min=0.166666
	ternary_z_max=0.3334

This is also why I do not construct it to solve quaternary system since we could envision that it would be extremely unefficient in terms of chemical potential searching. 

After running with:

	mpiexec.hydra -np 16 ./ising

You could see the resulted directory tree:
	
	GS_solutions/
	|-- hull.txt
	|-- x0.000000y0.000000
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.000000y0.166667
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.000000y0.333333
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.000000y0.666667
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.000000y0.833333
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.000000y1.000000
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.166667y0.666667
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.166667y0.833333
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.333333y0.666667
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.666667y0.000000
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.666667y0.166667
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.666667y0.333333
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.833333y0.000000
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	|-- x0.833333y0.166667
	|   |-- POSCAR
	|   `-- POSCAR_OUT
	`-- x1.000000y0.000000
	    |-- POSCAR
	    `-- POSCAR_OUT


## More on installation detail 

If for some reason, you really need to compile everything somewhere else yourself (not recommended). Here are some important notes

### on Ginar 

simply use the build.sh to build

### on MAC 
remember to put (or remove everything up to it)

`/usr/local/bin`

on the top of PATH variable, or alternatively remove anaconda and generally locally installed package (from the PATH)

How to install GS solver on new mac:
Firstly check brew list, brew uninstall all gcc stuffs, mpich, boost-mpi, ensure the following is set in the environment

alias g++="g++-4.8"  
alias gcc="gcc-4.8"  

export PATH="/usr/local/sbin:$PATH"  
export HOMEBREW_CC=gcc-4.8  
export HOMEBREW_CXX=g++-4.8

the key is to do 

```
brew install gcc@4.8  
brew reinstall mpich --build-from-source  
brew reinstall boost-mpi --build-from-source  
```
 
then using `build_ising_Mac.sh` from `/home/key01027/ising_solver/ising_119_some_debug/src/ising/` should work

if there is some hostname error, look at `https://stackoverflow.com/questions/23112515/mpich2-gethostbyname-failed`

To change your /etc/hosts hostname type the following: `sudo nano /etc/hosts`

and then add the line   `127.0.0.1       my_new_hostname`

### on Linux (Ginar)

**important**! If you use anaconda, you need to temporarily disable Anaconda to use compile GS solver.

Also ensure your PATH looks something like like this in Ginar:
```
/share/apps/intel/composer_xe_2015.3.187/bin/intel64:/share/apps/intel/composer_xe_2015.3.187/debugger/gdb/intel64_mic/bin:/home/key01027/.rubies/ruby-2.3.1/bin/:/home/key01027/.local/bin/:/home/key01027/ising_solver/boost_build/bin/:/share/apps/bin:/share/apps/intel/bin:/share/apps/intel/advisor_xe_2015.1.10.380555/bin64:/share/apps/intel//impi/5.0.3.048/intel64/bin:/share/apps/intel//itac/9.0.3.051/intel64/bin:/share/apps/intel/composer_xe_2015.3.187/bin/intel64:/share/apps/intel/composer_xe_2015.3.187/debugger/gdb/intel64_mic/bin:/opt/openmpi/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/java/latest/bin:/opt/rocks/bin:/opt/rocks/sbin:/opt/gridengine/bin/linux-x64
```

Modify the Makefile accordingly to compile (it is not easy and very system dependent). If you use intel compiler this two links should be useful [link1](https://software.intel.com/en-us/articles/building-boost-with-intel-c-compiler-150), [link2](http://stackoverflow.com/questions/39502034/how-do-i-install-boost-with-intel-compiler-and-intel-mpi)
