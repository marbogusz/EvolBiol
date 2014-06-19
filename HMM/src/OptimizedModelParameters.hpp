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
#include "SubstitutionModelBase.hpp"
#include "Maths.hpp"
#include "IndelModel.hpp"

typedef dlib::matrix<double,0,1> column_vector;

using namespace std;

namespace EBC
{

class OptimizedModelParameters
{
protected:

	vector<double> indelParameters;
	vector<double> substParameters;
	vector<double> divergenceTimes;
	double alpha;
	double divergenceBound;

	bool estimateIndelParams;
	bool estimateSubstParams;
	bool estimateAlpha;

	unsigned int indelCount;
	unsigned int substCount;
	unsigned int distCount;
	unsigned int optCount;

	//FIXME - static bound vectors!
	SubstitutionModelBase* sm;
	IndelModel* im;

	Maths* maths;

	void generateInitialIndelParameters();
	void generateInitialSubstitutionParameters();
	void generateInitialDistanceParameters();

public:
	OptimizedModelParameters(SubstitutionModelBase*, IndelModel*, unsigned int, bool, bool, bool, Maths*);

	unsigned int optParamCount();

	void toDlibVector(column_vector&,column_vector&,column_vector&);

	void fromDlibVector(const column_vector&);

	void setAlpha(double);

	void setUserIndelParams(vector<double>);

	void setUserSubstParams(vector<double>);

	void outputParameters();

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
};

} /* namespace EBC */

#endif /* OPTIMIZEDMODELPARAMETERS_H_ */