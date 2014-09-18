/*
 * NucleotideSubstitutionMode.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "models/NucleotideSubstitutionModel.hpp"

namespace EBC
{

NucleotideSubstitutionModel::NucleotideSubstitutionModel(Dictionary* dict, Maths* alg, unsigned int alpha, unsigned int pcnt) :
	SubstitutionModelBase(dict,alg,alpha,pcnt)
{
}


void NucleotideSubstitutionModel::calculateModel()
{
	this->buildSmatrix();
	this->setDiagonalMeans();
	this->doEigenDecomposition();
	//qMatrix, roots, u and v matrices

	//calculateGammaPtMatrices();
}

NucleotideSubstitutionModel::~NucleotideSubstitutionModel()
{
}

} /* namespace EBC */

