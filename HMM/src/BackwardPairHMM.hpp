/*
 * BackwardPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef BACKWARDPAIRHMM_HPP_
#define BACKWARDPAIRHMM_HPP_

#include "EvolutionaryPairHMM.hpp"

namespace EBC
{

class BackwardPairHMM: public EBC::EvolutionaryPairHMM
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
	BackwardPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict,  Definitions::ModelType model, bool banding,
			unsigned int bandPercentage, unsigned int rateCategories, Maths*);

	virtual ~BackwardPairHMM();

	double runBackwardAlgorithm();

	inline double* getMlParameters()
	{
		return this->mlParameters;
	}
};

} /* namespace EBC */
#endif /* BACKWARDPAIRHMM_HPP_ */
