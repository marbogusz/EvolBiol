/*
 * NucleotideSubstitutionMode.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "models/CodonModel.hpp"

namespace EBC
{

CodonModel::CodonModel(Dictionary* dict, Maths* alg, unsigned int alpha) :
	SubstitutionModelBase(dict,alg,alpha,Definitions::CodonM0ParamCount)
{
	this->parameters = new double[Definitions::CodonM0ParamCount];
}

void CodonModel::buildInitialQmatrix()
{

}

void CodonModel::calculateModel()
{
	//this->buildSmatrix();
	//this->setDiagonalMeans();
	//this->doEigenDecomposition();
	//qMatrix, roots, u and v matrices

	//calculateGammaPtMatrices();
}

CodonModel::~CodonModel()
{
	if (this->parameters)
		delete[] this->parameters;
}

} /* namespace EBC */

