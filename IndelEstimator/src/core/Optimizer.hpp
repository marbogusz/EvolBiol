/*
 * Optimizer.hpp
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#ifndef OPTIMIZER_HPP_
#define OPTIMIZER_HPP_

#include "core/OptimizedModelParameters.hpp"
#include "core/IOptimizable.hpp"

namespace EBC
{

class Optimizer
	{
	protected:
		column_vector initParams;
		column_vector lowerBounds;
		column_vector upperBounds;

		unsigned int paramsCount;

		OptimizedModelParameters* omp;
		IOptimizable* target;
		double accuracy;

		Definitions::OptimizationType optimizationType;

	public:
		Optimizer(OptimizedModelParameters* mp, IOptimizable* opt, Definitions::OptimizationType ot, double accuracy=Definitions::accuracyBFGS);
		virtual ~Optimizer();

		double optimize();

		void setTarget(IOptimizable* opt);

		double objectiveFunction(const column_vector& m);

		const column_vector objectiveFunctionDerivative(const column_vector& m);
	};

} /* namespace EBC */

#endif /* OPTIMIZER_HPP_ */
