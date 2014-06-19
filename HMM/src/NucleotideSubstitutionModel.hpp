/*
 * NucleotideSubstitutionModel.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef NUCL_MODEL_H_
#define NUCL_MODEL_H_

#include "Dictionary.hpp"
#include "Definitions.hpp"
#include "SubstitutionModelBase.hpp"
#include "Maths.hpp"
#include <cmath>

namespace EBC
{

class NucleotideSubstitutionModel : public EBC::SubstitutionModelBase
{
protected:

	virtual void buildSmatrix()=0;

public:

	NucleotideSubstitutionModel(Dictionary*, Maths*i, unsigned int, unsigned int);

	virtual ~NucleotideSubstitutionModel();

	void calculatePt();

};

} /* namespace EBC */
#endif /* MODEL_H_ */
