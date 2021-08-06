#include <fstream>
#include <string>
#include <unordered_map>
#include <math.h>
#include <thread>
#include <mutex>

#include "PoscarInfo.h"
#include "intWriteHead.h"

//File Constants
std::string POSCAR_LOC = "POSCAR";
std::string AFF_LOC = "AFF.csv";
std::string INT_OUTFILE_LOC = "ints.csv";
std::string GRAPH_OUTFILE_LOC = "plt.csv";

//Computation Constants
uint8_t H_MAX = 6;
uint8_t K_MAX = 6;
uint8_t L_MAX = 6;
float WAVELENGTH = 1.5f; ///angstroms
float BRAGG_TEMP_FACTOR = 0.2f; ///square angstroms (Generally between 0.2 & 0.9.  TODO: should be diff for each atom)
unsigned int NUM_PROCS = std::thread::hardware_concurrency();

//Graphing Constants
float GRID_SPACING = 0.01f; //spacing of evaluations of 2*Theta for profile functions, in degrees
float GAUSS_WIDTH = 0.5f; //width of gaussian profile at half intensity, in degrees.  2 * theta
float GAUSS_RANGE_MULT = 1.25f; //specifies range of 2theta to apply gaussian profiling to.  Multiplies GAUSS_WIDTH

void ComputeIntSet(hklSet*, unsigned int, const float[6], Poscar&, std::unordered_map<std::string, AFF>&, float*, float*, unsigned int);


int main(uint8_t argc, char* argv[]) {
	//Command line arguments
	for (uint8_t i = 0; i < argc; i++) {
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
		if (argv[i][0] == '-' && (argv[i][1] == 'h' || argv[i][1] == 'H')) //max magnitude of h index
			H_MAX = std::stoi(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'k' || argv[i][1] == 'K')) //max magnitude of k index
			K_MAX = std::stoi(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'l' || argv[i][1] == 'L')) //max magnitude of l index
			L_MAX = std::stoi(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'w' || argv[i][1] == 'W')) //wavelength in angstroms
			WAVELENGTH = std::stof(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 't' || argv[i][1] == 'T')) //bragg temp factor in sq. angstroms
			BRAGG_TEMP_FACTOR = std::stof(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'n' || argv[i][1] == 'N')) //number of processors to use
			NUM_PROCS = std::stoi(argv[i + 1]);
		//Graphing
		if (argv[i][0] == '-' && (argv[i][1] == 's' || argv[i][1] == 'S')) //spacing of 2 theta eval grid
			GRID_SPACING = std::stof(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'b' || argv[i][1] == 'B')) //breadth of gaussian profile about intensity spikes
			GAUSS_WIDTH = std::stof(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'm' || argv[i][1] == 'M')) //width multiplier about spikes to apply gauss profile to
			GAUSS_RANGE_MULT = std::stof(argv[i + 1]);
	}

	//General constants
	unsigned int totMillerCount = (2 * H_MAX + 1) * (2 * K_MAX + 1) * (2 * L_MAX + 1) - 1;
	unsigned int jobsPerNode = unsigned int(std::floor(totMillerCount / NUM_PROCS));
	unsigned int jobsPerNodeEdge = unsigned int(std::floor(totMillerCount / NUM_PROCS) + (totMillerCount % NUM_PROCS));

	//Get poscar, convert to direct.  Get the Sij's
	Poscar p("readAll", POSCAR_LOC);
	p.convertToCartesian();
	p.convertToDirect();
	float Sij[6] = { 0.0f };
	GetSij(p, Sij);

	//Set up atomic form factors hash table
	std::unordered_map <std::string, AFF> affHash;
	SetHashTable(AFF_LOC, affHash);

	//Get array of hklSets, their sizes, and thread pointers
	hklSet** allSets = (hklSet**)malloc(sizeof(hklSet*) * NUM_PROCS);
	unsigned int* allSizes = (unsigned int*)malloc(sizeof(unsigned int) * NUM_PROCS);
	std::thread** threads = (std::thread**)malloc(sizeof(std::thread*) * NUM_PROCS);
	for (uint8_t i = 0; i < NUM_PROCS; i++) {
		allSets[i] = (hklSet*)malloc(sizeof(hklSet) * jobsPerNode);
		allSizes[i] = jobsPerNode;
	}
	allSets[NUM_PROCS - 1] = (hklSet*)malloc(sizeof(hklSet) * jobsPerNodeEdge);
	allSizes[NUM_PROCS - 1] = jobsPerNodeEdge;

	//Set up eval grids for gaussian profiles
	float highestAng = 360.0f / pi * GetMaxTheta(H_MAX, K_MAX, L_MAX, WAVELENGTH, p.volume, Sij);
	unsigned int maxIndex = unsigned int(std::floor(highestAng / GRID_SPACING));
	///theta
	float* twoTheta = (float*)malloc(sizeof(float) * maxIndex);
	float gp = 0.0f;
	for (unsigned int i = 0; i < maxIndex; i++) {
		twoTheta[i] = gp;
		gp += GRID_SPACING;
	}
	///intensity
	float** intensity = (float**)malloc(sizeof(float*) * NUM_PROCS);
	for (uint8_t i = 0; i < NUM_PROCS; i++) {
		intensity[i] = (float*)malloc(sizeof(float) * maxIndex);
		for (unsigned int j = 0; j < maxIndex; j++)
			intensity[i][j] = 0.0f;
	}

	//Initialize array of hklSets
	uint8_t currentProcCount = 0;
	unsigned int currentMillerCount = 0;
	for (int8_t h = -H_MAX; h < H_MAX + 1; h++) {
		for (int8_t k = -K_MAX; k < K_MAX + 1; k++) {
			for (int8_t l = -L_MAX; l < L_MAX + 1; l++) {
				if (h == 0 && k == 0 && l == 0)
					continue;

				allSets[currentProcCount][currentMillerCount] = *new hklSet(h, k, l);
				if (currentMillerCount == allSizes[currentProcCount] - 1) {
					threads[currentProcCount] = new std::thread(ComputeIntSet,
																allSets[currentProcCount], allSizes[currentProcCount], ///(hkl), bragAng, int 
																Sij, p, affHash, ///misc
																twoTheta, intensity[currentProcCount], maxIndex); ///gaussian profile
					currentMillerCount = -1;
					currentProcCount++;
				}

				currentMillerCount++;
			}
		}
	}

	//Wait for all to finish
	for (uint8_t i = 0; i < NUM_PROCS; i++) {
		threads[i]->join();
	}

	//Get normalization constant, prep final output
	float norm = 0.0f;
	float* intFin = (float*)malloc(sizeof(float) * maxIndex);
	for (unsigned int i = 0; i < maxIndex; i++) {
		intFin[i] = 0.0f;
		for (uint8_t j = 0; j < NUM_PROCS; j++)
			intFin[i] += intensity[j][i];

		if (intFin[i] > norm)
			norm = intFin[i];
	}

	//Write results
	///hkl table
	std::ofstream outfileHKL(INT_OUTFILE_LOC);
	for (uint8_t i = 0; i < NUM_PROCS; i++) 
		for (unsigned int j = 0; j < allSizes[i]; j++)
			if(allSets[i][j].valid)
				outfileHKL << int(allSets[i][j].h) << ',' << int(allSets[i][j].k) << ',' << int(allSets[i][j].l) << ',' 
							<< 360.0f / pi * allSets[i][j].bragAng << ',' << allSets[i][j].intensity << "\n";
	outfileHKL.close();
	///graph file
	std::ofstream outfileGraph(GRAPH_OUTFILE_LOC);
	for (unsigned int i = 0; i < maxIndex; i++)
		 outfileGraph << twoTheta[i] << ',' << intFin[i] / norm << "\n";
	outfileGraph.close();

	return 0;
}



void ComputeIntSet(hklSet* hkls, unsigned int hklsSize, const float* Sij, Poscar& p, std::unordered_map<std::string, AFF>& affHash,
				   float* thetaGrid, float* intGrid, unsigned int gridSize) {
	float mult = 1.0f;
	for (int i = 0; i < hklsSize; i++) {
		///Get evaluation angle of this (hkl)
		float bragAng = BraggAngle(hkls[i].h, hkls[i].k, hkls[i].l, WAVELENGTH, Sij, p.volume);
		if (bragAng == 80085.0f) ///not intense enough wavelength - avoids problems with arcsine 
			continue;

		///Get intensity at this (hkl)
		float sinAng = std::sinf(bragAng);
		float cosDoubleAng = std::cosf(2.0f * bragAng);
		float angDep = (1.0f + cosDoubleAng * cosDoubleAng) / (sinAng * sinAng * std::cosf(bragAng));
		float strucDep = StructFactMag(hkls[i].h, hkls[i].k, hkls[i].l, WAVELENGTH, p.volume, Sij, p.atomCoords, affHash);
		float tempDep = std::pow(eul, (-2.0f * BRAGG_TEMP_FACTOR * sinAng * sinAng) / (WAVELENGTH * WAVELENGTH));
		float intensity = mult / (p.volume * p.volume) * angDep * strucDep * tempDep;

		///Add info to arrays
		hkls[i].valid = true;
		hkls[i].bragAng = bragAng;
		hkls[i].intensity = intensity;

		int leftmostGauss = int((360.0f / pi * bragAng - GAUSS_RANGE_MULT * GAUSS_WIDTH) / GRID_SPACING);
		int righmostGauss = int((360.0f / pi * bragAng + GAUSS_RANGE_MULT * GAUSS_WIDTH) / GRID_SPACING);
		for (int j = std::max(0, leftmostGauss); j < std::min(int(gridSize), righmostGauss); j++) 
			intGrid[j] += intensity * GaussianProfile(GAUSS_RANGE_MULT, thetaGrid[j], 360.0f / pi * bragAng);
	}

	return;
}