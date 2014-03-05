/*
 * GTRModel.h
 *
 *  Created on: Jan 17, 2014
 *      Author: root
 */

#ifndef GTRMODEL_H_
#define GTRMODEL_H_

#include "SubstitutionModel.hpp"

namespace EBC
{

class GTRModel: public EBC::SubstitutionModel
{
protected:

double *a, *b, *c, *d, *e, *f;
double scale;


public:
	GTRModel(Dictionary* dict, Maths* alg);

	void setParametersInMatrix();

	void summarize();
};

} /* namespace EBC */
#endif /* GTRMODEL_H_ */
