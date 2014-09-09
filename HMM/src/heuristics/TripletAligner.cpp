/*
 * TripletAligner.cpp
 *
 *  Created on: Sep 8, 2014
 *      Author: root
 */

#include <heuristics/TripletAligner.hpp>
#include <vector>
#include <array>

namespace EBC
{

array<vector<SequenceElement>, 3> TripletAligner::align()
{
	string& seq1 = inputSeqs->getRawSequenceAt(s1);
	string& seq2 = inputSeqs->getRawSequenceAt(s2);
	string& seq3 = inputSeqs->getRawSequenceAt(s3);

	GotohAlgorithm* algo = new GotohAlgorithm(seq2,seq1);

	algo->setDistance(distMat.getDistance(s1,s2));
	//scoring matrices etc....
	algo->run();
	auto firstPair = algo->getAlignment();

	algo->setDistance(distMat.getDistance(s2,s3));
	algo->setSequences(seq2,seq3);
	algo->run();
	auto secondPair = algo->getAlignment();

	//align
	string& anch1 = firstPair.first;
	string& anch2 = secondPair.first;

	string& p1al = firstPair.second;
	string& p2al = secondPair.second;

	DEBUG("Triplet aligner sequence 1 " << seq1);
	DEBUG("Triplet aligner sequence 2 " << seq2);
	DEBUG("Triplet aligner sequence 3 " << seq3);

	unsigned int alSize = std::max(anch1.size(), anch2.size());
	alSize *= 0.2;

	string tr1, tr2, tr3;
	tr1.reserve(alSize);
	tr3.reserve(alSize);
	tr2.reserve(alSize);


	unsigned int ctr1;
	unsigned int ctr2;

	ctr1 = ctr2 = 0;

	while(ctr1 < anch1.size() || ctr2 < anch2.size())
	{
		if(anch1[ctr1] == anch2[ctr2])
		{
			tr1 += p1al[ctr1];
			tr2 += anch1[ctr1];
			tr3 += p2al[ctr2];
			ctr1 ++;
			ctr2 ++;
		}
		else if (anch1[ctr1] == '-')
		{
			tr1 += p1al[ctr1];
			//gap
			tr2 += anch1[ctr1];
			tr3 += anch1[ctr1];
			ctr1 ++;

		}
		else
		{
			tr1 += anch2[ctr2];
			//gap
			tr2 += anch2[ctr2];
			tr3 += p2al[ctr2];
			ctr2 ++;
		}
	}

	array<vector<SequenceElement>, 3> triAlignment;
	triAlignment[0] = inputSeqs->getDictionary()->translate(tr1, false);
	triAlignment[1] = inputSeqs->getDictionary()->translate(tr2, false);
	triAlignment[2] = inputSeqs->getDictionary()->translate(tr3, false);

	DEBUG("Triplet alignment : ");
	DEBUG(tr1);
	DEBUG(tr2);
	DEBUG(tr3 << endl);

	return triAlignment;

}

TripletAligner::TripletAligner(Sequences* iSeq, array<unsigned int, 3> triplet, DistanceMatrix& dm) : inputSeqs(iSeq), distMat(dm)
{
	s1 = triplet[0];
	s2 = triplet[1];
	s3 = triplet[2];
}

} /* namespace EBC */
