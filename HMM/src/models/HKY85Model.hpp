/*
 * HKY85Model.h
 *
 *  Created on: Jan 17, 2014
 *      Author: root
 */

#ifndef HKY85MODEL_H_
#define HKY85MODEL_H_

#include "models/NucleotideSubstitutionModel.hpp"

namespace EBC
{

class HKY85Model: public EBC::NucleotideSubstitutionModel
{
protected:

double *k;


public:
	HKY85Model(Dictionary* dict, Maths* alg, unsigned int);


	void buildSmatrix();


	void summarize();

	void setParameters(const vector<double>&);
};

} /* namespace EBC */
#endif /* HKY85MODEL_H_ */
