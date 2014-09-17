/*
 * PMatrix.cpp
 *
 *  Created on: Sep 17, 2014
 *      Author: root
 */

#include <core/PMatrix.hpp>

namespace EBC
{

PMatrix::PMatrix(SubstitutionModelBase* m) : model(m),  matrixSize(m-getMatrixSize()) ,time(0),
		rateCategories(model->getRateCategories()), ptMatrices(rateCategories, nullptr)
{
	// TODO Auto-generated constructor stub
	matrixFullSize = matrixSize*matrixSize;
	this->fastPairGammaPt = new double[matrixFullSize];
}

PMatrix::~PMatrix()
{
	// TODO Auto-generated destructor stub
	for (int i = 0; i < ptMatrices.size(); i++)
	{
		delete [] ptMatrices[i];
	}
	delete [] fastPairGammaPt;
}

void PMatrix::setTime(double t)
{
		this->time = t;

		std::fill(fastPairGammaPt, fastPairGammaPt+matrixFullSize, 0);

		for(unsigned int i = 0; i< rateCategories; i++)
		{
			if (ptMatrices[i] != NULL)
				delete [] ptMatrices[i];
			ptMatrices[i] = this->model->calculatePt(i);

			for (int j=0; j< matrixFullSize; j++)
			{
				fastPairGammaPt[j] += ptMatrices[i][j] * model->gammaFrequencies[i];
			}

		}
}


double PMatrix::getPairSitePattern(array<unsigned int, 2>& nodes)
{
	return getPairSitePattern(nodes[0],nodes[1]);
}

double PMatrix::getPairSitePattern(unsigned int xi, unsigned int yi)
{
	return fastPairGammaPt[xi*matrixSize+yi];
}

double PMatrix::getTripleSitePattern(unsigned int root,
		array<unsigned int, 3>& nodes, PMatrix* pm2, PMatrix* pm3)
{
	PMatrix* pm1 = this;
	double lnl = 0;
	for(int rt = 0; rt< rateCategories; rt++)
	{
		lnl += getEquilibriumFreq(root) * pm1->getTransitionProb(root,nodes[0], rt) *
				pm2->getTransitionProb(root,nodes[1], rt) *
				pm2->getTransitionProb(root,nodes[1], rt) * model->gammaFrequencies[rt];
	}

	return lnl;
}

double PMatrix::getTransitionProb(unsigned int xi, unsigned int yi, unsigned int rateCat)
{
	return (ptMatrices[rateCat])[xi*matrixSize+yi];
}

inline double PMatrix::getEquilibriumFreq(unsigned int xi)
{
	return model->getEquilibriumFrequencies(xi);
}

void PMatrix::summarize()
{

	return model->summarize();
}

} /* namespace EBC */

