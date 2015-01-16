/*
 * HKY85Model.cpp
 *
 *  Created on: Jan 17, 2014
 *      Author: root
 */

#include "models/HKY85Model.hpp"

namespace EBC
{

HKY85Model::HKY85Model(Dictionary* dict, Maths* alg, unsigned int rates)
	: NucleotideSubstitutionModel(dict, alg, rates, Definitions::HKY85ParamCount)
{

	this->parameters = new double[this->paramsNumber];

	this->parameterHiBounds[0] = 5;
	this->parameterLoBounds[0] = 0.000001;

	//TODO - set within the model
}


void HKY85Model::setParameters(const vector<double>& par)
{
	for (unsigned int i = 0; i< paramsNumber; i++)
	{
		this->parameters[i] = par[i];
	}
}

void HKY85Model::buildSmatrix() {
	this->k = &parameters[0];

	int s = this->matrixSize;

	for (int i=0; i<s; i++)
		for (int j=0; j<s; j++)
		{
			if (i!=j)
			{
				qMatrix[(i*s)+j] = 1;
			}
		}
	qMatrix[1] = qMatrix[4]= qMatrix[11] = qMatrix[14] = parameters[0];

}

void HKY85Model::summarize()
{
	INFO("HKY85 model summary:");
	INFO("kappa " << *k );
	INFO("Frequencies (T C A G)");
	INFO(this->piFreqs[0] << '\t' << this->piFreqs[1] << '\t' << this->piFreqs[2] << '\t' << this->piFreqs[3] << '\t');
	summarizeRates();
}

} /* namespace EBC */
