/*
 * NucleotideSubstitutionModel.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef NUCL_MODEL_H_
#define NUCL_MODEL_H_

#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
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

	void calculateModel();

};

} /* namespace EBC */
#endif /* MODEL_H_ */
