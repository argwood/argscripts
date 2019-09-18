# argscripts
A collection of scripts that I have created for use during my PhD

**Wannier Wrap**: wraps Wannier centers from an xyz file into a periodic box 
                  - useful for Qbox output files, in which Wannier centers are 
				  unwrapped with respect to atom positions
				  
				  Requirements: g++, c++11
				  Compile: make
				  Usage: ./main input.xyz output.xyz bond_cutoff a_x a_y a_z  
