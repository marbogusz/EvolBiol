/*
 * OptimizedModelParameters.h
 *
 *  Created on: Jun 17, 2014
 *      Author: root
 */

#ifndef OPTIMIZEDMODELPARAMETERS_H_
#define OPTIMIZEDMODELPARAMETERS_H_

#include <vector>
#include <dlib/optimization.h>
#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
#include "models/IndelModel.hpp"

typedef dlib::matrix<double,0,1> column_vector;

using namespace std;

namespace EBC
{

class OptimizedModelParameters
{
protected:

	Maths* maths;

	SubstitutionModelBase* sm;
	IndelModel* im;

	vector<double> indelParameters;
	vector<double> substParameters;
	vector<double> divergenceTimes;
	double alpha;

	vector<double> indelHiBounds;

	bool estimateIndelParams;
	bool estimateSubstParams;
	bool estimateAlpha;
	bool estimateDivergence;

	unsigned int indelCount;
	unsigned int substCount;
	unsigned int distCount;
	unsigned int seqCount;
	unsigned int optCount;

	//FIXME - static bound vectors!



public:
	double divergenceBound;

	OptimizedModelParameters(SubstitutionModelBase*, IndelModel*, unsigned int, unsigned int, bool, bool, bool, bool, Maths*);

	unsigned int optParamCount();

	void useSubstitutionModelInitialParameters();

	void useIndelModelInitialParameters();

	void boundDivergenceBasedOnLambda(double lambda);

	void boundLambdaBasedOnDivergence(double time);

	void toDlibVector(column_vector&,column_vector&,column_vector&);

	void fromDlibVector(const column_vector&);

	void setAlpha(double);

	void setUserIndelParams(vector<double>);

	void setUserSubstParams(vector<double>);

	void setUserDivergenceParams(vector<double>);

	void setSingleDivergenceParam(unsigned int pos, double val){
		divergenceTimes[pos] = val;
	}

	void logParameters();

	void outputToConsole();

	double getDistanceBetween(unsigned int i, unsigned int j);

	void generateInitialIndelParameters();
	void generateInitialSubstitutionParameters();
	void generateInitialDistanceParameters();

	double getAlpha() const
	{
		return alpha;
	}

	double getDivergenceTime(unsigned int index)
	{
		return divergenceTimes[index];
	}

	vector<double> getIndelParameters()
	{
		return indelParameters;
	}

	vector<double> getSubstParameters()
	{
		return substParameters;
	}

	vector<double> getDivergenceTimes()
	{
		return divergenceTimes;
	}
};

} /* namespace EBC */

#endif /* OPTIMIZEDMODELPARAMETERS_H_ */
