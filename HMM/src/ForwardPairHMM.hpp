/*
 * ForwardPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef FORWARDPAIRHMM_HPP_
#define FORWARDPAIRHMM_HPP_

#include "EvolutionaryPairHMM.hpp"
#include <dlib/optimization.h>

typedef dlib::matrix<double,0,1> column_vector;



namespace EBC
{

class ForwardPairHMM: public EBC::EvolutionaryPairHMM
{

private:

	//BFGS optimization wrapper for dlib
	class BFGS
	{
	protected:
		column_vector initParams;
		column_vector lowerBounds;
		column_vector upperBounds;

		unsigned int paramsCount;

		ForwardPairHMM* parent;

		Definitions::OptimizationType optimizationType;

	public:
		BFGS(ForwardPairHMM* enclosing, Definitions::OptimizationType ot);
		virtual ~BFGS();
		void optimize();

		double objectiveFunction(const column_vector& m);

		const column_vector objectiveFunctionDerivative(const column_vector& m);
	};



protected:

	void initializeModels();

	void initializeStates();

	BFGS* bfgs;

	//Bound scale
	unsigned int bandFactor;
	unsigned int bandSpan;

	bool bandingEnabled;

	bool estimateSubstitutionParams;
	bool estimateIndelParams;
	bool estimateDivergence;

	vector<double> userIndelParameters;
	vector<double> userSubstParameters;

	double* optParameters;
	unsigned int optParametersCount;

	void setParameters();

	void getBandWidth()
	{
		this->bandSpan = ySize/(bandFactor);
		DEBUG("Band span " << bandSpan);
	}

	inline bool withinBand(unsigned int line, int position, unsigned int width)
	{
		int low = line - width;
		int high = line + width;
		bool result = ((position >= low) && (position <= high));

		//if(result == false)
		//{
		//	DEBUG("FALSE RESULT FOR l :" << line << " position " << position << " low " << low << " high " << high);
		//}

		return result;
	}


public:
	ForwardPairHMM(Sequences* inputSeqs, Definitions::ModelType model,std::vector<double> indel_params,
			std::vector<double> subst_params, Definitions::OptimizationType ot, bool banding,
			unsigned int bandPercentage, double evolDistance);

	virtual ~ForwardPairHMM();

	double runForwardAlgorithm();

	double runForwardIteration(const column_vector& m);

	inline double* getMlParameters()
	{
		return this->mlParameters;
	}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
