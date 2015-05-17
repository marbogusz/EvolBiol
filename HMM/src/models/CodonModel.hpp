/*
 * NucleotideSubstitutionModel.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef CODON_MODEL_H_
#define CODON_MODEL_H_

#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
#include <cmath>

namespace EBC
{

class CodonModel : public EBC::SubstitutionModelBase
{

protected:

	/*
	vector<unsigned int> piCodons;
	vector<unsigned int> kCodons;
	vector<unsigned int> wCodons;
	vector<unsigned int> kwCodons;
	*/
	//params[0] - kappa
	//params[1] - omega

	void buildInitialQmatrix();


public:

	CodonModel(Dictionary*, Maths*i, unsigned int);

	virtual ~CodonModel();

	void calculateModel();

};

} /* namespace EBC */
#endif /* MODEL_H_ */
