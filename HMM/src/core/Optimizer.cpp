/*
 * Optimizer.cpp
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#include <core/Optimizer.hpp>

namespace EBC
{

Optimizer::Optimizer(OptimizedModelParameters* mp, IOptimizable* opt, Definitions::OptimizationType ot) :
		optimizationType(ot), omp(mp), target(opt)
{


	paramsCount = omp->optParamCount();
	this->initParams.set_size(paramsCount);
	this->lowerBounds.set_size(paramsCount);
	this->upperBounds.set_size(paramsCount);

	omp->toDlibVector(initParams,lowerBounds,upperBounds);

	DEBUG("Numeric optimizer with " << paramsCount << " parameter(s) created");
	//cerr << "DLIB optimizer init with " << paramsCount << " parameters" << endl;
}

void Optimizer::setTarget(IOptimizable* opt)
{
	target = opt;
}

Optimizer::~Optimizer()
{
}

double Optimizer::objectiveFunction(const column_vector& bfgsParameters)
{
	//FIXME - address an issue of vector copying!
	omp->fromDlibVector(bfgsParameters);
	return target->runIteration();
}


const column_vector Optimizer::objectiveFunctionDerivative(const column_vector& bfgsParameters)
{
	column_vector results(this->paramsCount);
	return results;
}


double Optimizer::optimize()
{
	using std::placeholders::_1;
	std::function<double(const column_vector&)> f_objective= std::bind( &Optimizer::objectiveFunction, this, _1 );
	double likelihood;

	switch(optimizationType)
	{
		case Definitions::OptimizationType::BFGS:
		{
			likelihood = dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
					dlib::objective_delta_stop_strategy(1e-8),
					f_objective,
					derivative(f_objective),
					initParams,
					lowerBounds,
					upperBounds);
			break;
		}
		case Definitions::OptimizationType::BOBYQA:
		{
			likelihood = dlib::find_min_bobyqa(f_objective, initParams, paramsCount+4,
					lowerBounds,upperBounds, 0.05, 1e-7, 20000 );
			break;
		}
	}
	omp->fromDlibVector(initParams);
	return likelihood;
	//omp->outputParameters();
	//cout  << likelihood << "\n";

}


} /* namespace EBC */
