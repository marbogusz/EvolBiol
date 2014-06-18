/*
 * OptimizedModelParameters.cpp
 *
 *  Created on: Jun 17, 2014
 *      Author: root
 */

#include "OptimizedModelParameters.hpp"

namespace EBC
{

OptimizedModelParameters::OptimizedModelParameters(SubstitutionModelBase* sm, IndelModel* im, unsigned int dCount,
		bool ie, bool se, bool ae, Maths* m) : sm(sm), im(im), maths(m), indelParameters(im->getParamsNumber()),
		substParameters(sm->getParamsNumber()), divergenceTimes(dCount),
		estimateIndelParams(ie), estimateSubstParams(se), estimateAlpha(ae)
{
	indelCount = indelParameters.size();
	substCount = substParameters.size();
	distCount = dCount;
	optCount = (estimateSubstParams ? substCount : 0) +
			(estimateIndelParams ? indelCount : 0) + (distCount +
			estimateAlpha ? 1 : 0);
	this->divergenceBound  = 4;

	if(estimateIndelParams)
		generateInitialIndelParameters();
	if(estimateSubstParams)
		generateInitialSubstitutionParameters();
	generateInitialDistanceParameters();
}

} /* namespace EBC */

unsigned int EBC::OptimizedModelParameters::optParamCount()
{
	return optCount;
}

void EBC::OptimizedModelParameters::toDlibVector(column_vector& vals, column_vector& lbounds, column_vector& hbounds)
{
	unsigned int i;
	unsigned int ptr = 0;

	if(estimateSubstParams)
	{
		for (i=0; i < substCount; i++)
		{
			vals(i) = substParameters[i];
			//default probs bounds
			lbounds(i) = sm->getLoBound(i);
			//FIXME - provide external bounds????????????
			hbounds(i) = sm->getHiBound(i);
		}
		ptr += substCount;
	}
	if(estimateIndelParams)
	{
		for (i=0; i < indelCount; i++)
		{
			vals(i+ptr) = indelParameters[i];
			//default probs bounds
			lbounds(i+ptr) = im->getLoBound(i);
			//FIXME - provide external bounds????????????
			hbounds(i+ptr) = im->getHiBound(i);
		}
		ptr += indelCount;
	}
	if(estimateAlpha)
	{
		vals(ptr) = alpha;
		lbounds(ptr) = 0.000001;
		hbounds(ptr) = 99.999999;
		ptr++;
	}
	for (i=0; i < distCount; i++)
	{
		vals(i+ptr) = divergenceTimes[i];
		lbounds(i+ptr) = 0.000001;
		hbounds(i+ptr) = divergenceBound;
	}
}

void EBC::OptimizedModelParameters::generateInitialIndelParameters()
{
	for(unsigned int i=0; i< substCount; i++)
	{
		substParameters[i] = 0.2 + 0.1*maths->rndu();
	}
}

void EBC::OptimizedModelParameters::generateInitialSubstitutionParameters()
{
	for(unsigned int i=0; i< indelCount; i++)
	{
		indelParameters[i] = 0.2 + 0.1*maths->rndu();
	}
}

void EBC::OptimizedModelParameters::generateInitialDistanceParameters()
{
	for(unsigned int i=0; i< distCount; i++)
	{
		divergenceTimes[i] = 0.2 + 0.1*maths->rndu();
	}
}

void EBC::OptimizedModelParameters::fromDlibVector(const column_vector& vals)
{
	unsigned int i;
	unsigned int ptr = 0;
	if(estimateSubstParams)
	{
		for (i=0; i < substCount; i++)
		{
			substParameters[i] = vals(i);
		}
		ptr += substCount;
	}
	if(estimateIndelParams)
	{
		for (i=0; i < indelCount; i++)
		{
			indelParameters[i] = vals(i+ptr);
		}
		ptr += indelCount;
	}
	if(estimateAlpha)
	{
		alpha = vals(ptr);
		ptr++;
	}
	for (i=0; i < distCount; i++)
	{
		divergenceTimes[i] = vals(i+ptr);
	}
}

void EBC::OptimizedModelParameters::setAlpha(double a)
{
	alpha = a;
}

void EBC::OptimizedModelParameters::setUserIndelParams(
		vector<double> allocator)
{
	indelParameters = allocator;
}

void EBC::OptimizedModelParameters::setUserSubstParams(
		vector<double> allocator)
{
	substParameters = allocator;
}
