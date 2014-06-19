/*
 * GTRModel.cpp
 *
 *  Created on: Jan 17, 2014
 *      Author: root
 */

#include "GTRModel.hpp"

namespace EBC
{

GTRModel::GTRModel(Dictionary* dict, Maths* alg, unsigned int rates)
	: NucleotideSubstitutionModel(dict, alg, rates, Definitions::GTRParamCount)
{

	this->parameters = new double[this->paramsNumber];

	for (int i=0;i<5;i++)
	{
		this->parameterLoBounds[i] = 0.000001;
		this->parameterHiBounds[i] = 5;
	}
}

void GTRModel::setParameters(vector<double>& par)
{
	for (int i = 0; i< paramsNumber; i++)
	{
		this->parameters[i] = par[i];
	}
}


void GTRModel::buildSmatrix() {
	this->a = &parameters[0];
	this->b = &parameters[1];
	this->c = &parameters[2];
	this->d = &parameters[3];
	this->e = &parameters[4];
	this->f = &scale;

	//FIXME - deal with alpha somehow
	//this->time = parameters[5];

	int s = this->matrixSize;
	int i,j,k;


	this->qMatrix[(s-1)*s+2] = this->qMatrix[(s-2)*s+3] =  1.0;//scale;  /* r_AG = r_GA = 1. */ //FIXME - exponent
	for(i=0,k=0; i<s-1; i++) for (j=i+1; j<s; j++)
		if(i*s+j != 2*s+3)
			this->qMatrix[i*s+j] = this->qMatrix[j*s+i] = parameters[k++];

}

void GTRModel::summarize()
{
	cout << endl << "REV model summary:" << endl;
	cout << "a\tb\tc\td\te\tttime" << endl;
	cout << *a << "\t" << *b << "\t"<< *c << "\t"<< *d << "\t"<< *e << "\t"<< time << endl;
}

} /* namespace EBC */
