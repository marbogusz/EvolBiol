/*
 * PMatrix.hpp
 *
 *  Created on: Sep 17, 2014
 *      Author: root
 */

#ifndef PMATRIX_HPP_
#define PMATRIX_HPP_

#include "models\SubstitutionModelBase.hpp"
#include <vector>
#include <array>

using namespace std;

namespace EBC
{

class PMatrix
{
protected:

	SubstitutionModelBase* model;

	unsigned int matrixSize;

	double time;

	unsigned int rateCategories;

	vector<double*> ptMatrices;

	double* fastPairGammaPt;

	double matrixFullSize;

public:
	PMatrix(SubstitutionModelBase* m);
	virtual ~PMatrix();

void setTime(double t);

double getTransitionProb(unsigned int xi, unsigned int yi, unsigned int rateCat = 0);

inline double getEquilibriumFreq(unsigned int xi);

double getPairSitePattern(array<unsigned int, 2>& nodes);

double getPairSitePattern(unsigned int xi, unsigned int yi);

double getTripleSitePattern(unsigned int root,array<unsigned int, 3>& nodes, PMatrix* pm2, PMatrix* pm3);

void summarize();

};

} /* namespace EBC */

#endif /* PMATRIX_HPP_ */
