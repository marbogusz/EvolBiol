/*
 * HKY85Model.cpp
 *
 *  Created on: Jan 17, 2014
 *      Author: root
 */

#include "HKY85Model.hpp"

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


void HKY85Model::setParameters(vector<double>& par)
{
	for (int i = 0; i< paramsNumber; i++)
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
	cout << endl << "HKY85 model summary:" << endl;
	cout << "kappa\ttime" << endl;
	cout << *k << "\t" << time << endl;
}

} /* namespace EBC */
