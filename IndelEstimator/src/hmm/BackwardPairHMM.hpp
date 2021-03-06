/*
 * BackwardPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef BACKWARDPAIRHMM_HPP_
#define BACKWARDPAIRHMM_HPP_

#include "hmm/EvolutionaryPairHMM.hpp"
#include "hmm/ForwardPairHMM.hpp"


namespace EBC
{

class BackwardPairHMM: public EBC::EvolutionaryPairHMM
{

protected:

	//maximum posterior state;
	PairwiseHmmStateBase* MPstate;

	inline bool withinBand(unsigned int line, int position, unsigned int width)
	{
		int low = line - width;
		int high = line + width;
		bool result = ((position >= low) && (position <= high));
		return result;
	}


public:
	BackwardPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* smdl, IndelModel* imdl,
			Definitions::DpMatrixType mt, Band* bandObj = nullptr);

	virtual ~BackwardPairHMM();

	double runAlgorithm();

	void calculatePosteriors(ForwardPairHMM* fwd);

	double getAlignmentLikelihood(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, bool post, vector<vector<double> >& posteriors);

	void calculateMaximumPosteriorMatrix();

	//maximum posteriori alignment
	pair<string, string> getMPAlignment();

	pair<vector<double>*, pair<vector<unsigned char>*, vector<unsigned char>*> >
	getMPDWithPosteriors();

};

} /* namespace EBC */
#endif /* BACKWARDPAIRHMM_HPP_ */
