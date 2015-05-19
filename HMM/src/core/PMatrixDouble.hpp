/*
 * PMatrix.hpp
 *
 *  Created on: Sep 17, 2014
 *      Author: root
 */

#ifndef PMATRIXDB_HPP_
#define PMATRIXDB_HPP_

#include "models/SubstitutionModelBase.hpp"
#include "core/HmmException.hpp"
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

	inline double getPairSitePattern(int xi, int yi)
	{
		if( xi < 0 || yi <0 )
			return 0;
		return sitePatterns[xi][yi];
	}

	double getPairTransition(array<int, 2>& nodes);

	double getPairTransition(int xi, int yi);

	double getLogPairTransition(int xi, int yi);

	void summarize();

};

} /* namespace EBC */

#endif /* PMATRIX_HPP_ */
