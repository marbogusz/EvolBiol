/*
 * HKY85Model.cpp
 *
 *  Created on: Jan 17, 2014
 *      Author: root
 */

#include "HKY85Model.hpp"

namespace EBC
{

HKY85Model::HKY85Model(Dictionary* dict, Maths* alg) : SubstitutionModel(dict, alg)
{
	//6 elements to estimate - a, b, c, d, e, f=1
	//6th element is divergence time;
	this->paramsNumber = 2;
	this->parameters = new double[this->paramsNumber];
	//TODO - set within the model
	this->modelSetsOwnInitials = false;

}


void HKY85Model::setParametersInMatrix() {
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

		

	//this->qMatrix[(s-1)*s+2] = this->qMatrix[(s-2)*s+3] =  1.0;//scale;  /* r_AG = r_GA = 1. */ //FIXME - exponent
	//for(i=0,k=0; i<s-1; i++) for (j=i+1; j<s; j++)
	//	if(i*s+j != 2*s+3)
	//		this->qMatrix[i*s+j] = this->qMatrix[j*s+i] = parameters[k++];


}

} /* namespace EBC */
