#pragma once
#include <fstream>
#include <string>
#include <math.h>
#include <unordered_map>

#include "PoscarInfo.h"
#include "3DVect.h"

extern const float eul = 2.7182818f;
extern const float pi = 3.1415927f;

const float G_PROF_MULT = 0.93943727f; //2sqrt(ln(2))/sqrt(pi)
const float G_PROF_EXP_MULT = 2.7725887f; //4ln(2)

//Header file for powder diffraction pattern intensity writing
//Equations from SPREADSHEET SIMULATION OF X-RAY POWDER DIFFRACTION - H. W. G. SPRAGET 
//and
//https://en.wikipedia.org/wiki/Atomic_form_factor#cite_note-2
//with form factors taken directly from 
//http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php


//Atomic force factors for elements
struct AFF {
	std::string elem = "Empty";
	float a[4];
	float b[4];
	float c;

	void ElemPicker(int, float);
	inline AFF() { return; }
	AFF(std::string);
};
//appropriate csvLine is a string of "elem,a1,b1,a2,b2,a3,b3,a4,b4,c"
void AFF::ElemPicker(int num, float val) {
	if (num == 9)
		c = val;
	else if (num % 2 == 1)
		a[int(floor(num/2.0f))] = val;
	else if (num % 2 == 0)
		b[int(floor(num/2.0f)) - 1] = val;
	
	return;
}
AFF::AFF(std::string csvLine){
	csvLine += ',';

	std::string thisWord = "";
	for (int i = 0, initCounter = 0; i < csvLine.length(); i++) {
		if (csvLine[i] == ',') {
			if (initCounter == 0)
				elem = thisWord;
			else
				ElemPicker(initCounter, std::stof(thisWord));
			
			thisWord = "";
			initCounter++;
		}
		else
			thisWord += csvLine[i];
	}
}

void Sort(int8_t *arr, unsigned int n){ //optimized for small (i.e. length 3) lists
	///not usually part of a sorting algo.
	for(unsigned int i = 0; i < n; i++)
		arr[i] = std::abs(arr[i]);
	
	for(unsigned int i = 1; i < n; i++){
		uint8_t thisKey = arr[i];

		unsigned int j = i - 1;
		for(; j >= 0 && arr[j] > thisKey; j--)
			arr[j + 1] = arr[j];
		
		arr[j + 1] = thisKey;
	}
}

float Fact(unsigned int n){
	if(n == 0 || n == 1)
		return 1;
	return n*Fact(n - 1);
}
//variable for each (hkl) to make code somewhat readable
struct hklSet {
	int8_t h, k, l;
	float bragAng;
	float intensity;
	float mult;

	bool valid;

	hklSet() { return; }
	hklSet(int8_t, int8_t, int8_t);

	bool operator ==(const hklSet &) const;
};
hklSet::hklSet(int8_t h_, int8_t k_, int8_t l_) {
	h = h_;
	k = k_;
	l = l_;

	///Compute multiplicity
	unsigned int maxNonzeroReps = 0; unsigned int nNonzero = 0;
	int8_t tmpArr[3] = {h, k, l};
	for(unsigned int i = 0; i < 3; i++){
		unsigned int thisRep = 1;
		for(unsigned int j = i + 1; j < 3; j++)
			if(tmpArr[i] == tmpArr[j] && tmpArr[i]*tmpArr[j] != 0)
				thisRep ++;
		if(thisRep > maxNonzeroReps)
			maxNonzeroReps = thisRep;

		if(tmpArr[i] != 0)
			nNonzero ++;
	}
	
	mult = 1.0f;
	for(unsigned int i = 4 - nNonzero; i <= 3; i++)
		mult *= float(2.0f*i);
	mult = mult / Fact(maxNonzeroReps);

	valid = false;
	bragAng = 0.;
	intensity = 0.;
	return;
}
//Determines if two sets are equivalent
bool hklSet::operator ==(const hklSet &tst) const{
	int8_t cmpr1[3] = {h, k, l};
	int8_t cmpr2[3] = {tst.h, tst.k, tst.l};

	Sort(cmpr1, 3); ///also turns any negatives into positives
	Sort(cmpr2, 3);

	for(unsigned int i = 0; i < 3; i++){
		if(cmpr1[i] != cmpr2[i])
			return false;
	}
	return true;
}

void SetHashTable(std::string infileLoc, std::unordered_map <std::string, AFF> &map) {
	std::ifstream infile(infileLoc);
	for (int i = 0; ; i++) {
		if (infile.eof())
			break;

		std::string csvLine;
		infile >> csvLine;

		if (i == 0)
			continue;

		AFF thisAff(csvLine);
		map[thisAff.elem] = thisAff;
	}
}

bool Equiv(const int8_t arr1[], const int8_t arr2[]){
	int8_t cmpr1[3];
	int8_t cmpr2[3];
	for(unsigned int i = 0; i < 3; i++){
		cmpr1[i] = arr1[i];
		cmpr2[i] = arr2[i];
	}

	Sort(cmpr1, 3); ///also turns any negatives into positives
	Sort(cmpr2, 3);

	for (unsigned int i = 0; i < 3; i++) {
		if (cmpr1[i] != cmpr2[i])
			return false;
	}
	return true;
}

//arrSize is, at maximum, hklMax^3
void GetInequivHKLs(hklSet *hkls, unsigned int &arrSize, const unsigned int hklMax){
	
	//Setup arrays
	unsigned int currInd = 0;
	for (int8_t h = hklMax; h > 0; h--) {
		for (int8_t k = h; k >= 0; k--) {
			for (int8_t l = k; l >= 0; l--) {
				hkls[currInd] = hklSet(h, k, l);
				currInd ++;
			}
		}
	}

	arrSize = currInd;
	return;

	/*---Old extra safe alg: prob not necessary anymore...
	int8_t** hklArr = (int8_t**)malloc(sizeof(int8_t*) * arrSize);
	for (unsigned int i = 0; i < arrSize; i++) {
		hklArr[i] = (int8_t*)malloc(sizeof(int8_t) * 3);
		for (unsigned int j = 0; j < 3; j++)
			hklArr[i][j] = 0;
	}
	unsigned int currInd = 0;
	
	//Find all combos of inequiv (hkl)
	for(int8_t h = hklMax; h > 0; h--){
		for (int8_t k = hklMax; k >= 0; k--){
			for (int8_t l = hklMax; l >= 0; l--){
				int8_t test[3] = {h, k, l}; ///(hkl) to test against already existing ones

				bool add = true;
				for(unsigned int i = 0; i < currInd; i++){
					if(Equiv(test, hklArr[i])){
						add = false;
						break;
					}
				}
				if(add){
					for(unsigned int i = 0; i < 3; i++)
						hklArr[currInd][i] = test[i];
					currInd++;
				}
			}
		}
	}

	//Edit the actual hkl list by transforming arrays of (hkl) to actual (hkl) variables
	for(unsigned int i = 0; i < currInd; i++)
		hkls[i] = hklSet(hklArr[i][0], hklArr[i][1], hklArr[i][2]);
	arrSize = currInd;

	return;*/
}

float AtomStructFact(AFF &atomInfo, float theta, float wvlngth) {
	float sum = 0.0f;
	for (int i = 0; i < 4; i++) {
		float Qred = std::sinf(theta) / wvlngth; ///scattering vector reduced by 4pi 
		sum += atomInfo.a[i] * std::pow(eul, -atomInfo.b[i] * Qred * Qred) + atomInfo.c;
	}
	return sum;
}

//Sij = {S11, S22, S33, S12, S23, S31}
void GetSij(Poscar &pos, float *Sij) {
	vect aV(pos.superCellVectorA[0], pos.superCellVectorA[1], pos.superCellVectorA[2]);
	vect bV(pos.superCellVectorB[0], pos.superCellVectorB[1], pos.superCellVectorB[2]);
	vect cV(pos.superCellVectorC[0], pos.superCellVectorC[1], pos.superCellVectorC[2]);
	float al = angleBetween(bV, cV);
	float be = angleBetween(aV, cV);
	float ga = angleBetween(aV, bV);

	Sij[0] = (bV.len * cV.len * std::sinf(al)) * (bV.len * cV.len * std::sinf(al));
	Sij[1] = (aV.len * cV.len * std::sinf(be)) * (aV.len * cV.len * std::sinf(be));
	Sij[2] = (aV.len * bV.len * std::sinf(ga)) * (aV.len * bV.len * std::sinf(ga));
	Sij[3] = aV.len * bV.len * cV.len * cV.len * (std::cosf(al) * std::cosf(be) - std::cosf(ga));
	Sij[4] = aV.len * aV.len * bV.len * cV.len * (std::cosf(be) * std::cosf(ga) - std::cosf(al));
	Sij[5] = aV.len * bV.len * bV.len * cV.len * (std::cosf(ga) * std::cosf(al) - std::cosf(be));

	return;
}

float InterplSpacing(int h, int k, int l, const float Sij[6], float vol) {
	float dInvSq = (1.0f / (vol * vol)) * (Sij[0]*h*h + Sij[1]*k*k + Sij[2]*l*l + 2.0f*Sij[3]*h*k + 2.0f*Sij[4]*k*l + 2.0f*Sij[5]*h*l);
	return std::pow(dInvSq, -0.5);
}

int QuickGCD(int a, int b, bool init=false) {
	return (b == 0) ?  a : QuickGCD(b, a % b);
}
void NormalizeHKL(int &h, int &k, int &l) {
	int	norm = QuickGCD(QuickGCD(h, k), l);
	h /= std::abs(norm);
	k /= std::abs(norm);
	l /= std::abs(norm);
	return;
}

//Given a conditions, returns the max possible theta to evaluate intensities at
//The closer the |argument| to arcsine is to 1, the higher the theta that can be evaluated
float GetMaxTheta(int hLim, int kLim, int lLim, float wvlngth, float vol, const float Sij[6]) {
	float bestArg = 0.0f;
	for (int h = -hLim; h < hLim + 1; h++) {
		for (int k = -kLim; k < kLim + 1; k++) {
			for (int l = -lLim; l < lLim + 1; l++) {
				if (h == 0 && k == 0 && l == 0)
					continue;

				float frac = wvlngth / (2.0f * InterplSpacing(h, k, l, Sij, vol));
				if (std::abs(frac) > 1.0f)
					continue;
				if(std::abs(1.0f - frac) < std::abs(1.0f - bestArg))
					bestArg = frac;
			}
		}
	}
	return std::asinf(bestArg);
}

void NormalizeTheta(float &theta) {

	if (0.0f <= theta || theta <= 2.0f * pi)
		return;
	if (theta < 0.0f)
		theta += 2.0f * pi;
	else if (theta > 2.0f * pi)
		theta -= 2.0f * pi;

	NormalizeTheta(theta);
}
///returns -80085.0 if the angle is too great 
float BraggAngle(int h, int k, int l, float wvlngth, const float Sij[6], float vol) {
	float frac = wvlngth / (2.0f * InterplSpacing(h, k, l, Sij, vol));
	return (std::abs(frac) <= 1.0f) ? std::asinf(frac) : -80085.0f; 
}

//make sure atomCoords are in fractional coordinates before using this
float StructFactMag(int h, int k, int l, float wvlngth, float vol, const float Sij[6], 
					std::vector<Coords> &atomCoords, std::unordered_map <std::string, AFF> &affHash) {
	float sinPart = 0.0f;
	float cosPart = 0.0f;
	float theta = BraggAngle(h, k, l, wvlngth, Sij, vol);
	for (int i = 0; i < atomCoords.size(); i++) {
		float fi = AtomStructFact(affHash[atomCoords[i].atomType], theta, wvlngth);
		sinPart += fi * std::sinf(2.0f * pi * (h * atomCoords[i].a + k * atomCoords[i].b + l * atomCoords[i].c));
		cosPart += fi * std::cosf(2.0f * pi * (h * atomCoords[i].a + k * atomCoords[i].b + l * atomCoords[i].c));
	}
	return sinPart * sinPart + cosPart * cosPart;
}

float GaussianProfile(float width, float evalTheta, float bragTheta) {
	return (G_PROF_MULT / width) * std::pow(eul, -(G_PROF_EXP_MULT / (width * width)) *
											(2.0f * evalTheta - 2.0f * bragTheta) * (2.0f * evalTheta - 2.0f * bragTheta));
}