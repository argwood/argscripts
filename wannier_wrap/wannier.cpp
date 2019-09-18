#include "read_xyz.h"
#include "wannier.h"


double apply_pbc_x(double dx, struct Params params) {

	return dx - round(dx/params.a_x)*params.a_x;	

}

double apply_pbc_y(double dy, struct Params params) {

	return dy - round(dy/params.a_y)*params.a_y;	

}

double apply_pbc_z(double dz, struct Params params) {

	return dz - round(dz/params.a_z)*params.a_z;	

}
int count_wannier(vector<Snapshot> &traj, struct Params params, int snap, int i, vector<int> &assigned) {
	double distance;
	int E_count = 0;
	for (int j = 0; j < params.n_atoms; j++) {
		if (traj[snap].names[j] == "E") {
			distance = sqrt(pow(traj[snap].coords[i].x - traj[snap].coords[j].x,2) + pow(traj[snap].coords[i].y - traj[snap].coords[j].y,2) + pow(traj[snap].coords[i].z - traj[snap].coords[j].z,2));
			if (distance < params.cutoff) {
				E_count++;
				assigned.push_back(j);
			}
		}
	}
	return E_count;
}

bool check_bonded(double dx, double dy, double dz, struct Params params) {
	double distance;
	distance = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
	//cout << distance << endl;
	if (distance < params.cutoff) {
		return true;
	}
	else {
		return false;
	}
}

void wrap_centers(vector<int> unassigned, vector<int> undercoordinated, vector<Snapshot> &traj, struct Params params, int snap) {	
	for (int i = 0; i < unassigned.size(); i++) {						// loop through the unassigned WCs 
		int W_index = unassigned[i];
			
		for (int j=0; j < undercoordinated.size(); j++) {							//then loop through each undercoordinated WC
			int atom_index = undercoordinated[j];
			double dx, dy, dz;
			double dx_new, dy_new, dz_new;

			dx = traj[snap].coords[W_index].x - traj[snap].coords[atom_index].x;
			dy = traj[snap].coords[W_index].y - traj[snap].coords[atom_index].y;
			dz = traj[snap].coords[W_index].z - traj[snap].coords[atom_index].z;
				
			dx_new = apply_pbc_x(dx, params);									//try to apply PBC and then check distance between WC and atom
			dy_new = apply_pbc_y(dy, params);
			dz_new = apply_pbc_z(dz, params);
			
			bool bonded = check_bonded(dx_new, dy_new, dz_new, params);

			if (bonded) {


				traj[snap].coords[W_index].x -=  round(dx/(params.a_x))*params.a_x; 
				traj[snap].coords[W_index].y -=  round(dy/(params.a_y))*params.a_y; 
				traj[snap].coords[W_index].z -=  round(dz/(params.a_z))*params.a_z; 

			}
		}
	}
}


void wannier(vector<Snapshot> &traj, struct Params params) {
 	int Pb_E_count;
 	int Cd_E_count;
	int S_E_count;
	int Se_E_count;
	int H_E_count;
	int O_E_count;
	int I_E_count;
	int Cl_E_count;
	int In_E_count;
	int Zn_E_count;
	int Te_E_count;
	
	for (int snap = 0; snap < traj.size(); snap++) {
		vector<int> assigned;						// a vector that collects indices for all "assigned" wannier functions
		vector<int> unassigned;						// a vector that collects indices for all "unassigned" wannier functions
		vector<int> undercoordinated;
		int n_unassigned = 0;

		for (int i = 0; i < params.n_atoms; i++) {
			if (traj[snap].names[i] == "Pb") {
				Pb_E_count = count_wannier(traj, params, snap, i, assigned);
				if (Pb_E_count < params.nw_Pb) {
					cout << "Undercoordinated Pb! Only " << Pb_E_count << " Wannier Centers" << endl;
					undercoordinated.push_back(i);
				}
				else if (Pb_E_count > params.nw_Pb) {
					cout << "Overcoordinated Pb!" << endl;			// This isn't necessarily an error -- you can get this if a WC is shared between two atoms, which is fine
				}
			}
			else if (traj[snap].names[i] == "Cd") {
				Cd_E_count = count_wannier(traj, params, snap, i, assigned);
				if (Cd_E_count < params.nw_Cd) {
					cout << "Undercoordinated Cd! Only " << Cd_E_count << " Wannier Centers" << endl;
					undercoordinated.push_back(i);
				}
				else if (Cd_E_count > params.nw_Cd) {
					cout << "Overcoordinated Cd!" << endl;
				}
			}
			else if (traj[snap].names[i] == "S") {
				S_E_count = count_wannier(traj, params, snap, i, assigned);
				if (S_E_count < params.nw_S) {
					cout << "Undercoordinated S! Only " << S_E_count << " Wannier Centers" << endl;
					undercoordinated.push_back(i);
				}
				else if (S_E_count > params.nw_S) {
					cout << "Overcoordinated S!" << endl;
				}
			}
			else if (traj[snap].names[i] == "Se") {
				Se_E_count = count_wannier(traj, params, snap, i, assigned);
				if (Se_E_count < params.nw_Se) {
					cout << "Undercoordinated Se! Only " << Se_E_count << " Wannier Centers" << endl;
					undercoordinated.push_back(i);
				}
				else if (Se_E_count > params.nw_Se) {
					cout << "Overcoordinated Se!" << endl;
				}
			}
			else if (traj[snap].names[i] == "I") {
				I_E_count = count_wannier(traj, params, snap, i, assigned);
				if (I_E_count < params.nw_I) {
					cout << "Undercoordinated I! Only " << I_E_count << " Wannier Centers"  << endl;
					undercoordinated.push_back(i);
				}
				else if (I_E_count > params.nw_I) {
					cout << "Overcoordinated I!" << endl;
				}
			}
			else if (traj[snap].names[i] == "O") {
				O_E_count = count_wannier(traj, params, snap, i, assigned);
				if (O_E_count < params.nw_O) {
					cout << "Undercoordinated O! Only " << O_E_count << " Wannier Centers"  << endl;
					undercoordinated.push_back(i);
				}
				else if (O_E_count > params.nw_O) {
					cout << "Overcoordinated O!" << endl;
				}
			}
			else if (traj[snap].names[i] == "H") {
				H_E_count = count_wannier(traj, params, snap, i, assigned);
				if (H_E_count < params.nw_H) {
					cout << "Undercoordinated H! Only " << H_E_count << " Wannier Centers"  << endl;
					undercoordinated.push_back(i);
				}
				else if (H_E_count > params.nw_H) {
					cout << "Overcoordinated H!" << endl;
				}
			}
			else if (traj[snap].names[i] == "Cl") {
				Cl_E_count = count_wannier(traj, params, snap, i, assigned);
				if (Cl_E_count < params.nw_Cl) {
					cout << "Undercoordinated Cl! Only " << Cl_E_count << " Wannier Centers"  << endl;
					undercoordinated.push_back(i);
				}
				else if (Cl_E_count > params.nw_Cl) {
					cout << "Overcoordinated Cl!" << endl;
				}
			}
			else if (traj[snap].names[i] == "In") {
				In_E_count = count_wannier(traj, params, snap, i, assigned);
				if (In_E_count < params.nw_In) {
					cout << "Undercoordinated In! Only " << In_E_count << " Wannier Centers"  << endl;
					undercoordinated.push_back(i);
				}
				else if (In_E_count > params.nw_In) {
					cout << "Overcoordinated In!" << endl;
				}
			}
			else if (traj[snap].names[i] == "Zn") {
				Zn_E_count = count_wannier(traj, params, snap, i, assigned);
				if (Zn_E_count < params.nw_Zn) {
					cout << "Undercoordinated Zn! Only " << Zn_E_count << " Wannier Centers"  << endl;
					undercoordinated.push_back(i);
				}
				else if (Zn_E_count > params.nw_Zn) {
					cout << "Overcoordinated Zn!" << endl;
				}
			}
			else if (traj[snap].names[i] == "Te") {
				Te_E_count = count_wannier(traj, params, snap, i, assigned);
				if (Te_E_count < params.nw_Te) {
					cout << "Undercoordinated Te! Only " << Te_E_count << " Wannier Centers"  << endl;
					undercoordinated.push_back(i);
				}
				else if (Te_E_count > params.nw_Te) {
					cout << "Overcoordinated Te!" << endl;
				}
			}

		}
		
		for (int i = params.n_atoms_NP; i < params.n_atoms; i++) {					//loop through the wannier centers ("atoms" after the real atoms) 
			if (find(assigned.begin(), assigned.end(), i) == assigned.end()) {		//if the wannier center is unassigned (not in "assigned" vector), record it's index #
				n_unassigned ++ ;
				unassigned.push_back(i);
			}
		}
		wrap_centers(unassigned, undercoordinated, traj, params, snap);	
	}
}
