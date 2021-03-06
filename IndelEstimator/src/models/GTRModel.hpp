/*
 * GTRModel.h
 *
 *  Created on: Jan 17, 2014
 *      Author: root
 */

#ifndef GTRMODEL_H_
#define GTRMODEL_H_

#include "models/NucleotideSubstitutionModel.hpp"

namespace EBC
{

class GTRModel: public EBC::NucleotideSubstitutionModel
{
protected:

double *a, *b, *c, *d, *e, *f;
double scale;


public:
	GTRModel(Dictionary* dict, Maths* alg,unsigned int rates);

	void buildSmatrix();

	void summarize();

	void setParameters(const vector<double>&);
};

} /* namespace EBC */
#endif /* GTRMODEL_H_ */
