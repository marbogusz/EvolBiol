/*
 * OptimizedModelParameters.cpp
 *
 *  Created on: Jun 17, 2014
 *      Author: root
 */

#include "OptimizedModelParameters.h"

namespace EBC
{

OptimizedModelParameters::OptimizedModelParameters(unsigned int iParCount, unsigned int sParCount, unsigned int dCount,
		bool ie, bool se, bool ae) : indelParameters(iParCount), substParameters(sParCount), divergenceTimes(dCount),
				estimateIndelParams(ie), estimateSubstParams(se), estimateAlpha(ae)
{
	indelCount = iParCount;
	substCount = sParCount;
	distCount = dCount;
	optCount = estimateSubstParams ? substCount : 0 +
			estimateIndelParams ? indelCount : 0 + distCount +
			estimateAlpha ? 1 : 0;
}

} /* namespace EBC */

unsigned int EBC::OptimizedModelParameters::optParamCount()
{
	return optCount;
}

void EBC::OptimizedModelParameters::toDlibVector(column_vector& cv)
{
	for (unsigned int i = 0; i< optCount; i++)
	{
		cv(i) =
	}
}

void EBC::OptimizedModelParameters::fromDlibVector(const column_vector& cv)
{
}
