/*
 * HKY85Model.h
 *
 *  Created on: Jan 17, 2014
 *      Author: root
 */

#ifndef HKY85MODEL_H_
#define HKY85MODEL_H_

#include "SubstitutionModel.hpp"

namespace EBC
{

class HKY85Model: public EBC::SubstitutionModel
{
protected:

double *k;


public:
	HKY85Model(Dictionary* dict, Maths* alg);


	void setParametersInMatrix();


	void summarize();
};

} /* namespace EBC */
#endif /* HKY85MODEL_H_ */
