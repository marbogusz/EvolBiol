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

	public:
		BFGS(ForwardPairHMM* enclosing);
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
	ForwardPairHMM(Sequences* inputSeqs, bool optimize);

	virtual ~ForwardPairHMM();

	double runForwardAlgorithm();

	double runForwardIteration(const column_vector& m);
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
