#include "read_xyz.h"

void read_xyz(struct Params &params, vector<Snapshot> &traj, string input) {

	string name;
	double xx;
	double yy;
	double zz;
	int snapcount = 1;

	int n_Pb_ = 0;
	int n_Cd_ = 0;
	int n_S_ = 0;
	int n_Se_ = 0;
	int n_I_ = 0;
	int n_H_ = 0;
	int n_O_ = 0;
	int n_Cl_ = 0;
	int n_In_ = 0;
	int n_Zn_ = 0;
	int n_P_ = 0;
	int n_W_ = 0;
	int n_atoms_NP_;


	ifstream infile;
	const char * c = input.c_str();
	infile.open(c); 		

	string line;
	getline(infile,line);
	params.n_atoms = atoi(line.c_str());
	infile.seekg (0,ios::beg);
	
	while (true) {
		stringstream ss;
		vector <string> names;
		vector <double> x;
		vector <double> y;
		vector <double> z;
		vector <Coordinates> coords;
		ss.clear ();
		ss.str ("");

		getline(infile,line);
		getline(infile,line);
		for (int i = 0; i < params.n_atoms; i++) {
			getline(infile,line);
			ss.clear ();
			ss.str ("");
			ss << line;
			ss >> name >> xx >> yy >> zz;
			names.push_back(name);
			x.push_back(xx);
			y.push_back(yy);
			z.push_back(zz);
			
		}
		
		if ( infile.eof() ) break;
	
		for (int i=0; i < params.n_atoms; i++) {
			coords.push_back(Coordinates());
			coords[i].x = x[i];
			coords[i].y = y[i];
			coords[i].z = z[i];
		}

		traj.push_back(Snapshot());
		traj[snapcount-1].names = names;
		traj[snapcount-1].coords = coords;
		traj[snapcount-1].snapcount = snapcount;
		
		snapcount++;
	}

	for (int i=0; i < params.n_atoms; i++) {
		if (traj[0].names[i] == "Pb") {
			n_Pb_ ++;
		}
		else if (traj[0].names[i] == "Cd") {
			n_Cd_ ++;
		}
		else if (traj[0].names[i] == "S") {
			n_S_ ++;
		}
		else if (traj[0].names[i] == "Se") {
			n_Se_ ++;
		}
		else if (traj[0].names[i] == "I") {
			n_I_ ++;
		}
		else if (traj[0].names[i] == "H") {
			n_H_ ++;
		}
		else if (traj[0].names[i] == "O") {
			n_O_ ++;
		}
		else if (traj[0].names[i] == "Cl") {
			n_Cl_ ++;
		}
		else if (traj[0].names[i] == "In") {
			n_In_ ++;
		}
		else if (traj[0].names[i] == "Zn") {
			n_Zn_ ++;
		}
		else if (traj[0].names[i] == "P") {
			n_P_ ++;
		}
		else if (traj[0].names[i] == "E") {
			n_W_ ++;
		}
	}
 	
	params.n_Pb = n_Pb_;
	params.n_Cd = n_Cd_;
	params.n_S = n_S_;
	params.n_Se = n_Se_;
	params.n_I = n_I_;
	params.n_H = n_H_;
	params.n_O = n_O_;
	params.n_Cl = n_Cl_;
	params.n_In = n_In_;
	params.n_Zn = n_Zn_;
	params.n_P = n_P_;
	params.n_W = n_W_;
	n_atoms_NP_ = n_Pb_ + n_S_ + n_I_ + n_Cd_ + n_Se_ + n_H_ + n_O_ + n_Cl_ + n_In_ + n_Zn_ + n_P_ ;
	params.n_atoms_NP = n_atoms_NP_;
	infile.close();
}


void write_xyz(vector<Snapshot> &traj, struct Params params, string output) {
	ofstream outfile;

	const char * c = output.c_str();
	outfile.open(c); 		

	for (int snap = 0; snap < traj.size(); snap++) {
		outfile << params.n_atoms << "\n";
		outfile << "MD snapshot # " << snap + 1 << ", wrapped WCs\n";
		for (int i = 0; i < params.n_atoms; i++) {
			outfile <<  traj[snap].names[i] << "\t" << traj[snap].coords[i].x << "\t" << traj[snap].coords[i].y << "\t" << traj[snap].coords[i].z << "\n";
		}	
	
	
	}
	outfile.close();
}
