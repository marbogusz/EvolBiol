/*
 * GTRModel.cpp
 *
 *  Created on: Jan 17, 2014
 *      Author: root
 */

#include "models/GTRModel.hpp"

namespace EBC
{

GTRModel::GTRModel(Dictionary* dict, Maths* alg, unsigned int rates)
	: NucleotideSubstitutionModel(dict, alg, rates, Definitions::GTRParamCount)
{
	for (int i=0;i<5;i++)
	{
		//FIXME - move magic numbers to definitions
		this->parameterLoBounds[i] = 0.000001;
		this->parameterHiBounds[i] = 5;
	}
}

void GTRModel::setParameters(const vector<double>& par)
{
	for (unsigned int i = 0; i< paramsNumber; i++)
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

	int s = this->matrixSize;
	int i,j,k;


	this->qMatrix[(s-1)*s+2] = this->qMatrix[(s-2)*s+3] =  1.0;
	for(i=0,k=0; i<s-1; i++) for (j=i+1; j<s; j++)
		if(i*s+j != 2*s+3)
			this->qMatrix[i*s+j] = this->qMatrix[j*s+i] = parameters[k++];

}

void GTRModel::summarize()
{
	cout << endl << "REV model summary:" << endl;
	cout << "a\tb\tc\td\te" << endl;
	cout << *a << "\t" << *b << "\t"<< *c << "\t"<< *d << "\t"<< *e << endl;
	//cout << "Frequencies" << endl;
	//cout << this->piFreqs[0] << '\t' << this->piFreqs[1] << '\t' << this->piFreqs[2] << '\t' << this->piFreqs[3] << '\t' << endl;
	//cout << "Eigenvalues" << endl;
	//cout << roots[0] << ' ' << roots[1] << ' ' << roots[2] << ' ' << roots[3] << endl << endl;

	summarizeRates();
}

} /* namespace EBC */
