/*
 * SubstitutionModelBase.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "models/SubstitutionModelBase.hpp"

namespace EBC
{

SubstitutionModelBase::SubstitutionModelBase(Dictionary* dict, Maths* alg, unsigned int rateCategories, unsigned int parameter_count)
	: dictionary(dict), maths(alg), rateCategories(rateCategories), paramsNumber(parameter_count),
	  parameterHiBounds(parameter_count), parameterLoBounds(parameter_count)
{
	this->matrixSize = dict->getAlphabetSize();
	this->matrixFullSize = matrixSize*matrixSize;
	this->allocateMatrices();
	this->meanRate=0.0;
	this-> roots = new double[matrixSize];
	//no alpha provided
	this->alpha = 0;
	//FIXME - maybe initial alpha set to a value ?
}


void SubstitutionModelBase::allocateMatrices()
{
	this->qMatrix = new double[matrixFullSize];
	this->uMatrix = new double[matrixFullSize];
	this->vMatrix = new double[matrixFullSize];
	this->squareRoots = new double[matrixFullSize];

	if (rateCategories > 0)
	{
		gammaFrequencies = new double[rateCategories];
		gammaRates = new double[rateCategories];
	}
	else
	{
		gammaFrequencies = new double[1];
		gammaRates = new double[1];
		gammaFrequencies[0] = gammaRates[0] = 1;
	}

	this->sitePatterns = new double*[matrixSize+1];
	for (int i =0; i<= matrixSize; i++ )
		sitePatterns[i] = new double[matrixSize+1];
}

void SubstitutionModelBase::calculateSitePatterns()
{
	//includes gaps - does not discard missing data!
	for (int i =0; i<= matrixSize; i++ )
		for (int j =0; j<= matrixSize; j++ )
		{
			if (i == matrixSize)
				sitePatterns[i][j]  = log(this->piFreqs[j]);
			else if (j == matrixSize)
				sitePatterns[i][j]  = log(this->piFreqs[i]);
			else
				sitePatterns[i][j] = log(this->getPXiYi(i,j));
		}
	sitePatterns[matrixSize][matrixSize] = 0;
}

void SubstitutionModelBase::destroyMatrices()
{
	//TODO delete all
}

void SubstitutionModelBase::doEigenDecomposition()
{
	std::fill(roots, roots+matrixSize, 0);
	std::fill(uMatrix, uMatrix+matrixFullSize, 0);
	std::fill(vMatrix, vMatrix+matrixFullSize, 0);
	std::fill(squareRoots, squareRoots+matrixFullSize, 0);
	this->maths->eigenQREV(qMatrix, piFreqs, matrixSize, roots, uMatrix, vMatrix, squareRoots);
}

void SubstitutionModelBase::calculateGammaPtMatrices()
{
	double *tmpRoots, *tmpUroots;
	unsigned int i;
	int matrixCount = rateCategories == 0 ? 1 : rateCategories;

	std::vector<double*> pMatVector(matrixCount);

	for(i = 0; i< matrixCount; i++)
	{

		tmpRoots = maths->expLambdaT(roots, time*gammaRates[i], matrixSize);
		tmpUroots = maths->matrixByDiagonalMultiply(uMatrix, tmpRoots, matrixSize);
		//tmpPmatrix = this->maths->matrixMultiply(tmpUroots, vMatrix, matrixSize);
		pMatVector[i] = this->maths->matrixMultiply(tmpUroots, vMatrix, matrixSize);
		delete[] tmpRoots;
		delete[] tmpUroots;
	}
	//now we have a vector of p-mats, combine into 1 matrix since the probs are equal in this model
	for (i=1; i<matrixCount; i++)
	{
		maths->matrixAppend(pMatVector[0],pMatVector[i],matrixSize);
		delete[] pMatVector[i];
	}
	this->pMatrix = pMatVector[0];
	if(rateCategories != 0)
	{
		this->maths->matrixScale(pMatrix,gammaFrequencies[0],matrixSize);
	}
}


void SubstitutionModelBase::setDiagonalMeans()
{
		int i,j;
		double sum;

		maths->matrixByDiagonalMultiplyMutable(qMatrix,piFreqs,matrixSize);
		meanRate = 0;

		for (i=0; i< this->matrixSize; i++)
		{
			qMatrix[i*matrixSize+i] = 0;
			sum = 0.0;
			for (j=0; j < this->matrixSize; j++)
			{
				sum -= qMatrix[i*matrixSize + j];
			}
			qMatrix[i*matrixSize + i] = sum;
			meanRate -= sum*piFreqs[i];
		}
		for (i=0; i< this->matrixSize; i++)
		{
				for (j=0; j < this->matrixSize; j++)
				{
					qMatrix[i*matrixSize + j] /= meanRate;
				}
		}

}

SubstitutionModelBase::~SubstitutionModelBase()
{
	destroyMatrices();
	if (this->gammaFrequencies !=NULL )
		delete[] gammaFrequencies;
	if (this->gammaRates !=NULL )
		delete[] gammaRates;

}

void SubstitutionModelBase::calculateGamma()
{
	int useMedian = 0;
	if(!(alpha <=0 || rateCategories == 0))
	{
		this->maths->DiscreteGamma(gammaFrequencies, gammaRates, alpha, alpha, rateCategories, useMedian);
		//DEBUGV(gammaRates,5);
		//DEBUGV(gammaFrequencies,5);

	}
}

void SubstitutionModelBase::setObservedFrequencies(double* observedFrequencies)
{
	//Size and order is clear
		this->piFreqs = observedFrequencies;

}

double SubstitutionModelBase::getPXiYi(unsigned int xi, unsigned int yi)
{
	//FIXME - RATES
		return piFreqs[xi]*pMatrix[(xi*matrixSize)+yi];
}


double SubstitutionModelBase::getQXi(unsigned int xi)
{
	return piFreqs[xi];
}

double SubstitutionModelBase::getPattern(unsigned int xi, unsigned int yi)
{
	return this->sitePatterns[xi][yi];
}

} /* namespace EBC */


