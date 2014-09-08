/*
 * TripletAligner.cpp
 *
 *  Created on: Sep 8, 2014
 *      Author: root
 */

#include <heuristics/TripletAligner.hpp>

namespace EBC
{

void TripletAligner::align()
{
	GotohAlgorithm* algo = new GotohAlgorithm(inputSeqs->getRawSequenceAt(s2),inputSeqs->getRawSequenceAt(s1));
	//scoring matrices etc....
	algo->run();
}

TripletAligner::TripletAligner(Sequences* iSeq, array<unsigned int, 3>& triplet) : inputSeqs(iSeq)
{
	s1 = triplet[0];
	s2 = triplet[1];
	s3 = triplet[2];
}

} /* namespace EBC */
