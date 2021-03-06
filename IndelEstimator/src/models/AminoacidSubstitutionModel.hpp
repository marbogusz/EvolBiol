/*
 * AminoacidSubstitutionModel.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef AA_MODEL_H_
#define AA_MODEL_H_

#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
#include <cmath>

namespace EBC
{

class AminoacidSubstitutionModel : public EBC::SubstitutionModelBase
{
protected:

	bool eigenDecomposed;

	double maxRate;

public:

	AminoacidSubstitutionModel(Dictionary*, Maths*, unsigned int, Definitions::aaModelDefinition);

	virtual ~AminoacidSubstitutionModel();

	void calculateModel();

	void summarize();

	void setParameters(const vector<double>&){}

	void setObservedFrequencies(double* observedFrequencies) {};
};

} /* namespace EBC */
#endif /* MODEL_H_ */
