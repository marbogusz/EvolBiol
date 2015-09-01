/*
 * PMatrix.hpp
 *
 *  Created on: Sep 17, 2014
 *      Author: root
 */

#ifndef PMATRIXDB_HPP_
#define PMATRIXDB_HPP_

#include "models/SubstitutionModelBase.hpp"
#include "core/ProgramException.hpp"
#include "core/PMatrix.hpp"
#include <vector>
#include <array>

using namespace std;

namespace EBC
{

//TODO - split into 2 regular pt and pairwise
class PMatrixDouble : public PMatrix
{
protected:

	double* fastPairGammaPt;
	double* fastLogPairGammaPt;

	double ** sitePatterns;

	void calculatePairSitePatterns();

public:
	PMatrixDouble(SubstitutionModelBase* m);
	virtual ~PMatrixDouble();

	void calculate();

	inline double getPairSitePattern(unsigned int xi, unsigned int yi)
	{
		return sitePatterns[xi][yi];
	}

	double getPairTransition(array<unsigned int, 2>& nodes);

	double getPairTransition(unsigned int xi, unsigned int yi);

	double getLogPairTransition(unsigned int xi, unsigned int yi);

	void summarize();

};

} /* namespace EBC */

#endif /* PMATRIX_HPP_ */
