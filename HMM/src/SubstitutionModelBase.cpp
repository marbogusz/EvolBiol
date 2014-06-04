/*
 * SubstitutionModelBase.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "SubstitutionModelBase.hpp"

namespace EBC
{

SubstitutionModelBase::SubstitutionModelBase(Dictionary* dict, Maths* alg, unsigned int rateCategories)
	: dictionary(dict), algebra(alg), rateCategories(rateCategories)
{
	this->matrixSize = dict->getAlphabetSize();
	this->matrixFullSize = matrixSize*matrixSize;
	this->allocateMatrices();
	this->meanRate=0.0;
	this-> roots = new double[matrixSize];
}


void SubstitutionModelBase::allocateMatrices()
{
	this->qMatrix = new double[matrixFullSize];
	this->uMatrix = new double[matrixFullSize];
	this->vMatrix = new double[matrixFullSize];
	this->squareRoots = new double[matrixFullSize];
}

void SubstitutionModelBase::destroyMatrices()
{
	//TODO delete all
}


void SubstitutionModelBase::setDiagonalMeans()
{
		int i,j;
		double sum;

		algebra->matrixByDiagonalMultiply(qMatrix,piFreqs,matrixSize);
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

} /* namespace EBC */


