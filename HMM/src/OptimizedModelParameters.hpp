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

	bool estimateIndelParams;
	bool estimateSubstParams;
	bool estimateAlpha;

	unsigned int indelCount;
	unsigned int substCount;
	unsigned int distCount;
	unsigned int optCount;

public:
	OptimizedModelParameters(unsigned int, unsigned int, unsigned int, bool, bool, bool);

	unsigned int optParamCount();

	void toDlibVector(column_vector&);
	void fromDlibVector(const column_vector&);

};

} /* namespace EBC */

#endif /* OPTIMIZEDMODELPARAMETERS_H_ */
