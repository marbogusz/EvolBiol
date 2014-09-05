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
	double divergenceBound;

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


	void generateInitialIndelParameters();
	void generateInitialSubstitutionParameters();
	void generateInitialDistanceParameters();

public:
	OptimizedModelParameters(SubstitutionModelBase*, IndelModel*, unsigned int, unsigned int, bool, bool, bool, bool, Maths*);

	unsigned int optParamCount();

	void toDlibVector(column_vector&,column_vector&,column_vector&);

	void fromDlibVector(const column_vector&);

	void setAlpha(double);

	void setUserIndelParams(vector<double>);

	void setUserSubstParams(vector<double>);

	void setUserDivergenceParams(vector<double>);

	void outputParameters();

	double getDistanceBetween(unsigned int i, unsigned int j);

	double getAlpha() const
	{
		return alpha;
	}

	double getDivergenceTime(unsigned int index)
	{
		return divergenceTimes[index];
	}

	const vector<double>& getIndelParameters() const
	{
		return indelParameters;
	}

	const vector<double>& getSubstParameters() const
	{
		return substParameters;
	}

	vector<double>& getDivergenceTimes()
	{
		return divergenceTimes;
	}
};

} /* namespace EBC */

#endif /* OPTIMIZEDMODELPARAMETERS_H_ */
