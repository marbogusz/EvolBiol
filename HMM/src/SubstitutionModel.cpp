/*
 * SubstitutionModel.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "SubstitutionModel.hpp"

namespace EBC
{

SubstitutionModel::SubstitutionModel(Dictionary* dict, Maths* alg)
{
	this->dictionary = dict;
	this->algebra = alg;
	this->matrixSize = dict->getAlphabetSize();
	this->matrixFullSize = matrixSize*matrixSize;
	this->allocateMatrices();
	this->meanRate=0.0;
	this-> roots = new double[matrixSize];
	//set the log mode to true ??
	this->logMode = false;
}


void SubstitutionModel::allocateMatrices()
{
	this->qMatrix = new double[matrixFullSize];
	this->uMatrix = new double[matrixFullSize];
	this->vMatrix = new double[matrixFullSize];
	this->squareRoots = new double[matrixFullSize];
}

void SubstitutionModel::destroyMatrices()
{
	//TODO delete all
}

void SubstitutionModel::doEigenDecomposition()
{
	//DEBUG("SubstitutionModel::doEigenDecomposition");
	std::fill(roots, roots+matrixSize, 0);
	std::fill(uMatrix, uMatrix+matrixFullSize, 0);
	std::fill(vMatrix, vMatrix+matrixFullSize, 0);
	std::fill(squareRoots, squareRoots+matrixFullSize, 0);

	this->algebra->eigenQREV(qMatrix, piFreqs, matrixSize, roots, uMatrix, vMatrix, squareRoots);
}

double* SubstitutionModel::calculatePt(double t)
{
	//DEBUG("SubstitutionModel::calculateSubstitutionModelPt");
	//exponentiate
	algebra->expLambdaT(roots, (logMode?exp(t):t), matrixSize);
	//algebra->expLambdaT(roots, logMode?log(t):t, matrixSize);
	//algebra->vectorMultiply(roots,t,matrixSize);
	//roots[0] = 0;
	//U.lambda exp
	algebra->matrixByDiagonalMultiply(uMatrix, roots, matrixSize);
	if (pMatrix != NULL)
	{
		delete[] pMatrix;
	}
	pMatrix = this->algebra->matrixMultiply(uMatrix, vMatrix, matrixSize);
	//FIXME -  make sure to delete pmatrix after using it!!!!
	return pMatrix;
}

void SubstitutionModel::calculatePt()
{
	this->calculatePt(this->time);
}


void SubstitutionModel::setDiagMeans()
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

SubstitutionModel::~SubstitutionModel()
{
	destroyMatrices();
}

void SubstitutionModel::setObservedFrequencies(double* observedFrequencies)
{
	//Size and order is clear
		this->piFreqs = this->observedFrequencies = observedFrequencies;

}

double SubstitutionModel::getPXiYi(unsigned int xi, unsigned int yi)
{
	if (logMode)
	{
		return  log(piFreqs[xi]) + log(pMatrix[(xi*matrixSize)+yi]);
	}
	else
	{
		return piFreqs[xi]*pMatrix[(xi*matrixSize)+yi];
	}
}


double SubstitutionModel::getQXi(unsigned int xi)
{
	return logMode? log(piFreqs[xi]) : piFreqs[xi];
}

} /* namespace EBC */


