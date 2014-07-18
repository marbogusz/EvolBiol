/*
 * AminoacidSubstitutionModel.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef AA_MODEL_H_
#define AA_MODEL_H_

#include "Dictionary.hpp"
#include "Definitions.hpp"
#include "SubstitutionModelBase.hpp"
#include "Maths.hpp"
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

	void calculatePt();

	void summarize();

	void setParameters(const vector<double>&){}
};

} /* namespace EBC */
#endif /* MODEL_H_ */
