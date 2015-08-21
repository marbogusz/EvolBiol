/*
 * BrentOptimizer.hpp
 *
 *  Created on: Aug 5, 2015
 *  Based on John D. Cook's implementation
 *  http://www.codeproject.com/Articles/30201/Optimizing-a-Function-of-One-Variable
 */

#ifndef CORE_BRENTOPTIMIZER_HPP_
#define CORE_BRENTOPTIMIZER_HPP_

#include "core/OptimizedModelParameters.hpp"
#include "core/IOptimizable.hpp"


namespace EBC {

class BrentOptimizer
{
protected:

	OptimizedModelParameters* omp;
	IOptimizable* target;
	double accuracy;

public:
	BrentOptimizer(OptimizedModelParameters* mp, IOptimizable* opt, double accuracy=Definitions::accuracyBFGS);
	double optimize();
	void setTarget(IOptimizable* opt);
	double objectiveFunction(double x);

	double getAccuracy() const {
		return accuracy;
	}

	void setAccuracy(double accuracy) {
		this->accuracy = accuracy;
	}
};

} /* namespace EBC */

#endif /* CORE_BRENTOPTIMIZER_HPP_ */
