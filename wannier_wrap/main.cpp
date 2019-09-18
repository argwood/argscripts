#include "read_xyz.h"
#include "wannier.h"


/**
 * The Wannier Wrap code will wrap Wannier centers into a periodic box for systems with the following elements:
 * Pb, Cd, S, Se, I, H, O, Cl, In, Zn, Te. 
 * Other elements may be added to suit the user's need by specifying the expected number of Wannier centers surrounded 
 * by that element in an ideal, neutral system (may be pseudopotential-dependent).
 * This script is intended for nanoparticles, nanoplatelets and other isolated or 2D systems with strong ionic character.
 *
 * Usage: ./main input_file.xyz output_file.xyz cutoff a_x a_y a_z
 *
 * input_file.xyz should include both the atoms and the Wannier centers, with Wannier centers labeled "E".
 * Cutoff should be tuned to be greater than the atom-WC bonds but smaller than atom-atom bonds. 
 * a_x, a_y, and a_z are the size of your box in Angstrom. 
 */


int main(int argc, char* argv[]) {

	if (argc < 7) {
		cerr << "Usage: " << argv[0] << " XYZ_IN XYZ_OUT cutoff a_x a_y a_z" << endl;
		return 1;
	}

	string input, output;//, cutoff_, a_x_, a_y_, a_z_;
	double cutoff_, a_x_, a_y_, a_z_;
	input = string(argv[1]);
	output = string(argv[2]);
	cutoff_ = stod(argv[3]);
	a_x_ = stod(argv[4]);
	a_y_ = stod(argv[5]);
	a_z_ = stod(argv[6]);
	struct Params params; 

	params.cutoff = cutoff_;
	params.a_x = a_x_;
	params.a_y = a_y_;
	params.a_z = a_z_;

	vector<Snapshot> traj; // this vector has all info for each snapshot, including names
	
	read_xyz(params, traj, input);
	
	wannier(traj, params);
	
	write_xyz(traj, params, output);

	return 0;
}
