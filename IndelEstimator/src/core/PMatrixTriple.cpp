/*
 * PMatrixTriple.cpp
 *
 *  Created on: Sep 17, 2014
 *      Author: root
 */

#include <core/PMatrixTriple.hpp>

namespace EBC
{

PMatrixTriple::PMatrixTriple(SubstitutionModelBase* m) : PMatrix(m)
{

}

PMatrixTriple::~PMatrixTriple()
{
}

void PMatrixTriple::calculate()
{
	if (time != 0)
	{
		for(unsigned int i = 0; i< rateCategories; i++)
		{
			if (ptMatrices[i] != NULL)
				delete [] ptMatrices[i];
			ptMatrices[i] = this->model->calculatePt(time, i);

		}
	}
	else
		throw ProgramException("PMatrixTriple : attempting to calculate p(t) with t set to 0");
}


double PMatrixTriple::getTripleSitePattern(unsigned int root,
		const array<unsigned char, 3>& nodes, PMatrixTriple* pm2, PMatrixTriple* pm3)
{
	PMatrixTriple* pm1 = this;
	double prob = 0;
	for(int rt = 0; rt< rateCategories; rt++)
	{
		prob += getEquilibriumFreq(root) * pm1->getTransitionProb(root,nodes[0], rt) *
				pm2->getTransitionProb(root,nodes[1], rt) *
				pm3->getTransitionProb(root,nodes[2], rt) * model->gammaFrequencies[rt];
	}

	return prob;
}

double PMatrixTriple::getTransitionProb(unsigned int xi, unsigned int yi, unsigned int rateCat)
{
	if (yi >= matrixSize)
		return 1;
	return (ptMatrices[rateCat])[xi*matrixSize+yi];
}

void PMatrixTriple::summarize()
{
	cout << "P(t) matrix summary :" << endl;
	cout << "Divergence time : " << time << endl;

}

} /* namespace EBC */


