/*
 * ForwardPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef FORWARDPAIRHMM_HPP_
#define FORWARDPAIRHMM_HPP_

#include "EvolutionaryPairHMM.hpp"

namespace EBC
{

class ForwardPairHMM: public EBC::EvolutionaryPairHMM
{

protected:

	vector<double> userIndelParameters;
	vector<double> userSubstParameters;

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
		return result;
	}


public:
	ForwardPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict,  Definitions::ModelType model, bool banding,
			unsigned int bandPercentage, unsigned int rateCategories, double alpha, Maths*);

	virtual ~ForwardPairHMM();

	double runForwardAlgorithm();

	inline double* getMlParameters()
	{
		return this->mlParameters;
	}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
