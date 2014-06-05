/*
 * NucleotideSubstitutionMode.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "NucleotideSubstitutionModel.hpp"

namespace EBC
{

NucleotideSubstitutionModel::NucleotideSubstitutionModel(Dictionary* dict, Maths* alg, unsigned int alpha) :
	SubstitutionModelBase(dict,alg,alpha)
{
}

void NucleotideSubstitutionModel::doEigenDecomposition()
{
	std::fill(roots, roots+matrixSize, 0);
	std::fill(uMatrix, uMatrix+matrixFullSize, 0);
	std::fill(vMatrix, vMatrix+matrixFullSize, 0);
	std::fill(squareRoots, squareRoots+matrixFullSize, 0);
	this->maths->eigenQREV(qMatrix, piFreqs, matrixSize, roots, uMatrix, vMatrix, squareRoots);
}

void NucleotideSubstitutionModel::calculatePt()
{
	double *tmpRoots, *tmpUroots, *tmpPmatrix;
	unsigned int i;
	int matrixCount = rateCategories == 0 ? 1 : rateCategories;

	std::vector<double*> pMatVector(rateCategories);

	this->buildSmatrix();
	this->setDiagonalMeans();
	this->doEigenDecomposition();
	//qMatrix, roots, u and v matrices

	for(i = 0; i< matrixCount; i++)
	{

		tmpRoots = maths->expLambdaT(roots, time*gammaRates[i], matrixSize);
		tmpUroots = maths->matrixByDiagonalMultiply(uMatrix, tmpRoots, matrixSize);
		tmpPmatrix = this->maths->matrixMultiply(tmpUroots, vMatrix, matrixSize);
		pMatVector.push_back(tmpPmatrix);
		delete[] tmpRoots;
		delete[] tmpUroots;
	}
	//now we have a vector of p-mats, combine into 1 matrix since the probs are equal in this model
	for (i=1; i<matrixCount; i++)
	{
		maths->matrixAppend((double*)pMatVector[0],(double*)pMatVector[i],matrixSize);
		delete[] (double*)pMatVector[i];
	}
	this->pMatrix = pMatVector[0];
	if(rateCategories != 0)
	{
		this->maths->matrixScale(pMatrix,gammaFrequencies[0],matrixSize);
	}
}

NucleotideSubstitutionModel::~NucleotideSubstitutionModel()
{
}

} /* namespace EBC */


