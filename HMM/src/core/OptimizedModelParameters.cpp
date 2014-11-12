/*
 * OptimizedModelParameters.cpp
 *
 *  Created on: Jun 17, 2014
 *      Author: root
 */

#include "core/OptimizedModelParameters.hpp"

namespace EBC
{

OptimizedModelParameters::OptimizedModelParameters(SubstitutionModelBase* sm, IndelModel* im, unsigned int sCount, unsigned int dCount,
		bool se, bool ie, bool ae, bool de, Maths* m) : maths(m), sm(sm), im(im), indelParameters(im != NULL ? im->getParamsNumber() : 0),
		substParameters(sm != NULL ? sm->getParamsNumber() : 0), divergenceTimes(dCount),
		estimateIndelParams(ie), estimateSubstParams(se), estimateAlpha(ae), estimateDivergence(de)
{
	indelCount = indelParameters.size(); //-1;  //FIXME - hack!!!
	substCount = substParameters.size();
	seqCount = sCount;
	distCount = dCount;
	optCount = (estimateSubstParams ? substCount : 0) +
			(estimateIndelParams ? indelCount : 0) + (estimateDivergence ? distCount : 0) + (estimateAlpha ? 1 : 0);

	this->divergenceBound  = Definitions::divergenceBound;

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
	if(estimateDivergence)
	{
		for (i=0; i < distCount; i++)
		{
			vals(i+ptr) = divergenceTimes[i];
			lbounds(i+ptr) = 0.000001;
			hbounds(i+ptr) = divergenceBound;
		}
	}
}

void EBC::OptimizedModelParameters::generateInitialSubstitutionParameters()
{
	for(unsigned int i=0; i< substCount; i++)
	{
		substParameters[i] = 0.2 + 0.1*maths->rndu();
	}
	FileLogger::DebugLogger() << "Initial substitution parameters:\n";
	FileLogger::DebugLogger() << substParameters;
}

void EBC::OptimizedModelParameters::generateInitialIndelParameters()
{
	for(unsigned int i=0; i< indelCount; i++)
	{
		indelParameters[i] = 0.05 + 0.1*maths->rndu();
		//indelParameters[i+1] = 0.5; //FIXME
	}
	FileLogger::DebugLogger() << "Initial indel parameters:\n";
	FileLogger::DebugLogger() << indelParameters;
}

void EBC::OptimizedModelParameters::generateInitialDistanceParameters()
{
	for(unsigned int i=0; i< distCount; i++)
	{
		divergenceTimes[i] = 0.2 + 0.1*maths->rndu();
	}

	FileLogger::DebugLogger() << "Initial divergence times:\n";
	FileLogger::DebugLogger() << divergenceTimes;
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
			//indelParameters[i+1] = 0.5; //FIXME
		}
		ptr += indelCount;
	}
	if(estimateAlpha)
	{
		alpha = vals(ptr);
		ptr++;
	}
	if(estimateDivergence)
	{
		for (i=0; i < distCount; i++)
		{
			divergenceTimes[i] = vals(i+ptr);
		}
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

void EBC::OptimizedModelParameters::setUserDivergenceParams(
		vector<double> allocator)
{
	divergenceTimes = allocator;
}

void EBC::OptimizedModelParameters::setUserSubstParams(
		vector<double> allocator)
{
	substParameters = allocator;
}

void EBC::OptimizedModelParameters::logParameters()
{
	//for (auto p : substParameters)
	//		std::cerr << p  << '\t';
	//for (auto p : indelParameters)
	//		std::cerr << p  << '\t';
	//for (auto p : divergenceTimes)
	//	std::cerr << p  << '\t';
	//if(this->estimateAlpha)
	//	std::cerr << alpha;
	//std::cerr << std::endl;
	FileLogger::InfoLogger() << substParameters;
	FileLogger::InfoLogger() << indelParameters;
	FileLogger::InfoLogger() << divergenceTimes;
	if(this->estimateAlpha)
		FileLogger::InfoLogger() << "Alpha: " << alpha << "\n";

}

void EBC::OptimizedModelParameters::outputToConsole()
{
	for (auto p : substParameters)
			std::cout << p  << '\t';
	for (auto p : indelParameters)
			std::cout << p  << '\t';
	for (auto p : divergenceTimes)
		std::cout << p  << '\t';
	if(this->estimateAlpha)
		std::cout << alpha;
	std::cout << std::endl;
}


double EBC::OptimizedModelParameters::getDistanceBetween(unsigned int i,
		unsigned int j)
{
	unsigned int total  = seqCount;
	unsigned int a,b;
	a=i;
	b=j;
	if (i>j)
	{
		b=i;
		a=j;
	}
	unsigned int index;
	if (i==j) return 0;
	else
	{
		index = (b - a - 1) + (a*total) - (((1+a)/2.0)*(a*1.0));
		return divergenceTimes[index];
	}
}
