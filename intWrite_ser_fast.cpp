#include <fstream>
#include <string>
#include <unordered_map>
#include <math.h>
#include <mutex>

#include "PoscarInfo.h"
#include "intWriteHead.h"

//File Constants
std::string POSCAR_LOC = "POSCAR";
std::string AFF_LOC = "AFF.csv";
std::string INT_OUTFILE_LOC = "ints.csv";
std::string GRAPH_OUTFILE_LOC = "plt.csv";

//Computation Constants
int HKL_MAX = 30;
float WAVELENGTH = 1.5f; ///angstroms
float BRAGG_TEMP_FACTOR = 0.2f; ///square angstroms (Generally between 0.2 & 0.9.  TODO: should be diff for each atom)

//Graphing Constants
float GRID_SPACING = 0.01f; //spacing of evaluations of 2*Theta for profile functions, in degrees
float GAUSS_WIDTH = 0.5f; //width of gaussian profile at half intensity, in degrees.  2 * theta
float GAUSS_RANGE_MULT = 1.25f; //specifies range of 2theta to apply gaussian profiling to.  Multiplies GAUSS_WIDTH

int main(int argc, char* argv[]) {
	//Command line arguments
	for (int i = 0; i < argc; i++) {
		//File io
		if (argv[i][0] == '-' && (argv[i][1] == 'p' || argv[i][1] == 'P')) //(in) poscar location
			POSCAR_LOC = argv[i + 1];
		if (argv[i][0] == '-' && (argv[i][1] == 'a' || argv[i][1] == 'A')) //(in) atomic force factor location
			AFF_LOC = argv[i + 1];
		if (argv[i][0] == '-' && (argv[i][1] == 'i' || argv[i][1] == 'I')) //(out) intensity table location
			INT_OUTFILE_LOC = argv[i + 1];
		if (argv[i][0] == '-' && (argv[i][1] == 'g' || argv[i][1] == 'G')) //(out) graph table location
			GRAPH_OUTFILE_LOC = argv[i + 1];
		//Computational 
		if (argv[i][0] == '-' && (argv[i][1] == 'h' || argv[i][1] == 'H' || 
								  argv[i][1] == 'k' || argv[i][1] == 'K' || 
							      argv[i][1] == 'l' || argv[i][1] == 'L')) //max magnitude of h, k, l index
			HKL_MAX = std::stoi(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'w' || argv[i][1] == 'W')) //wavelength in angstroms
			WAVELENGTH = std::stof(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 't' || argv[i][1] == 'T')) //bragg temp factor in sq. angstroms
			BRAGG_TEMP_FACTOR = std::stof(argv[i + 1]);
		//Graphing
		if (argv[i][0] == '-' && (argv[i][1] == 's' || argv[i][1] == 'S')) //spacing of 2 theta eval grid
			GRID_SPACING = std::stof(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'b' || argv[i][1] == 'B')) //breadth of gaussian profile about intensity spikes
			GAUSS_WIDTH = std::stof(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'm' || argv[i][1] == 'M')) //width multiplier about spikes to apply gauss profile to
			GAUSS_RANGE_MULT = std::stof(argv[i + 1]);
	}

	//Get poscar, convert to direct.  Get the Sij's
	Poscar p("readAll", POSCAR_LOC);
	p.convertToCartesian();
	p.convertToDirect();
	float Sij[6] = { 0.0f };
	GetSij(p, Sij);

	//Set up atomic form factors hash table
	std::unordered_map <std::string, AFF> affHash;
	SetHashTable(AFF_LOC, affHash);

	//Set up evaluation grid of 0 <= 2*theta <= highest 2*theta possible
	float highestAng = 360.0f / pi * GetMaxTheta(HKL_MAX, HKL_MAX, HKL_MAX, WAVELENGTH, p.volume, Sij); 
	float *twoTheta = (float*)malloc(int(highestAng / GRID_SPACING) * sizeof(float));
	float *intensity = (float*)malloc(int(highestAng / GRID_SPACING) * sizeof(float));
	float gp = 0.0f;
	for (int i = 0; i < int(highestAng / GRID_SPACING); i++) {
		twoTheta[i] = gp;
		intensity[i] = 0.0f;
		gp += GRID_SPACING;
	}

	//Get array of {h, k, l} that are inequivelent from (1,0,0) up to (HKL_MAX, HKL_MAX, HKL_MAX)
	unsigned int arrSize = HKL_MAX * HKL_MAX * HKL_MAX;
	hklSet *hklList = (hklSet*)malloc(arrSize * sizeof(hklSet));
	GetInequivHKLs(hklList, arrSize, HKL_MAX); ///also edits arrSize

	//Sum over inequiv array  Write outfile in-time
	std::ofstream outfile(INT_OUTFILE_LOC, 'w');
	outfile << "h,k,l,mult,2Theta,intensity\n";

	for(unsigned int i = 0; i < arrSize; i++){
		///Get evaluation angle of this (hkl)
		float bragAng = BraggAngle(hklList[i].h, hklList[i].k, hklList[i].l, WAVELENGTH, Sij, p.volume);
		if (bragAng == 80085.0f || bragAng < 0.0f) ///not intense enough wavelength - avoids problems with arcsine 
			continue;

		///Compute intensity of peak at this (hkl).  Write to raw intensity file
		float angDep = (1.0f + std::cosf(2.0f * bragAng) * std::cosf(2.0f * bragAng)) / (std::sinf(bragAng) * std::sinf(bragAng) * std::cosf(bragAng));
		float strucDep = StructFactMag(hklList[i].h, hklList[i].k, hklList[i].l, WAVELENGTH, p.volume, Sij, p.atomCoords, affHash);
		float tempDep = std::pow(eul, (-2.0f * BRAGG_TEMP_FACTOR * std::sinf(bragAng) * std::sinf(bragAng)) / (WAVELENGTH * WAVELENGTH));
		float thisIntensity = hklList[i].mult / (p.volume * p.volume) * angDep * strucDep * tempDep;
		outfile << int(hklList[i].h) << ',' << int(hklList[i].k) << ',' << int(hklList[i].l) << ',' << hklList[i].mult << ',' << 
				   bragAng * 360.0f / pi << ',' << thisIntensity << "\n";

		//Apply gaussian profile to important range of grid spacing
		for (int i = std::max(0, int((360.0f / pi * bragAng - GAUSS_RANGE_MULT * GAUSS_WIDTH) / GRID_SPACING));
			i < std::min(int(highestAng / GRID_SPACING), int((360.0f / pi * bragAng + GAUSS_RANGE_MULT * GAUSS_WIDTH) / GRID_SPACING)); i++)
			intensity[i] += thisIntensity * GaussianProfile(GAUSS_WIDTH, twoTheta[i], bragAng * 360.0f / pi);

	}
	outfile.close();

	//Normalize then write summed gaussian profiles to file.  Finish
	float norm = intensity[0];
	for (int i = 1; i < int(highestAng / GRID_SPACING); i++)
		if (intensity[i] > norm)
			norm = intensity[i];
	
	std::ofstream outfile2(GRAPH_OUTFILE_LOC);
	for (int i = 0; i < int(highestAng / GRID_SPACING); i++)
		outfile2 << twoTheta[i] << ',' << intensity[i] / norm << "\n";
	outfile2.close();

	return 0;
}

