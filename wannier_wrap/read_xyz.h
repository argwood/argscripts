#ifndef READ_XYZ_H
#define READ_XYZ_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

struct Coordinates {					//each atom has and x y and z coordinate
	double x;
	double y;
	double z;
};

struct Snapshot  {
	vector<Coordinates> coords;			//each snapshot consists of a vector of all atom coordinates,
	vector<string> names;				//a vector of atom types,
	int snapcount;						//and a counter for which snap it is
};

struct Params {

	int n_atoms;				//total number of atoms, including wannier centers as atoms
	int n_atoms_NP;				//total number of atoms, excluding wannier centers
	int n_Pb; 					//number of Pb atoms
	int n_Cd; 					//number of Cd atoms
	int n_S;					//number of S atoms
	int n_Se;					//number of Se atoms
	int n_I;					//number of I atoms
	int n_H;					//number of H atoms
	int n_O;					//number of O atoms
	int n_Cl;					//number of Cl atoms
	int n_In;					//number of In atoms
	int n_Zn;					//number of Zn atoms
	int n_P;					//number of Te atoms

	int n_W;					//number of Wannier centers
	const int nw_S = 4;			//number of expected wannier centers around each atom type
	const int nw_Se = 4;	
	const int nw_Pb = 6;
	const int nw_Cd = 9;
	const int nw_I = 9;
	const int nw_H = 0;
	const int nw_O = 4;
	const int nw_Cl = 4;
	const int nw_In = 10;
	const int nw_Zn = 9;
	const int nw_P = 4;     

	/* Below are parameters that you might need to adjust for your system~ */
	
	double cutoff;
	double a_x;
	double a_y;
	double a_z;
};
void read_xyz(struct Params &params, vector<Snapshot> &traj, string input); 
void write_xyz(vector<Snapshot> &traj, struct Params params, string output);
#endif
