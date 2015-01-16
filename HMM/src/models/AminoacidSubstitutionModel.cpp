/*
 * AminoacidSubstitutionMode.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "models/AminoacidSubstitutionModel.hpp"

namespace EBC
{

AminoacidSubstitutionModel::AminoacidSubstitutionModel(Dictionary* dict, Maths* alg, unsigned int alpha, Definitions::aaModelDefinition modelDef) :
	SubstitutionModelBase(dict,alg,alpha, Definitions::AAParamCount), eigenDecomposed(false)
{
	//FIXME - implement +F model with custom frequencies

	this->maxRate = modelDef.maxRate;
	this->piFreqs = new double[matrixSize];
	this->piLogFreqs = new double[matrixSize];

	std::copy(modelDef.aaRates,modelDef.aaRates + this->matrixFullSize, this->qMatrix);
	std::copy(modelDef.aaFreqs,modelDef.aaFreqs + this->matrixSize, this->piFreqs);

	for(unsigned int i = 0; i< this->matrixSize; i++)
	{
			piLogFreqs[i] = log(piFreqs[i]);
	}

	this->setDiagonalMeans();
	this->doEigenDecomposition();
	eigenDecomposed=true;
}


void AminoacidSubstitutionModel::calculateModel()
{
	if(this->eigenDecomposed == false)
	{
		this->setDiagonalMeans();
		this->doEigenDecomposition();
		eigenDecomposed=true;
	}
	//calculateGammaPtMatrices();
}

AminoacidSubstitutionModel::~AminoacidSubstitutionModel()
{
	//delete[] piFreqs;
	//delete[] piLogFreqs;
}

void AminoacidSubstitutionModel::summarize()
{
	summarizeRates();
}

} /* namespace EBC */


