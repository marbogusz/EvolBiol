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
	: NucleotideSubstitutionModel(dict, alg, rates)
{

	//6 parameters - 4 frequencies + kappa + time
	//2 parameters to estimate - k and t
	this->paramsNumber = 2;
	this->parameters = new double[this->paramsNumber];
	//TODO - set within the model

}

void HKY85Model::buildSmatrix() {
	this->k = &parameters[0];
	this->time = parameters[1];

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
