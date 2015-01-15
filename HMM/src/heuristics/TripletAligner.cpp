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

void TripletAligner::assembleFromPairs(pair<string, string>& p1,
		pair<string, string>& p2)
{
	string& anch1 = p1.second;
	string& anch2 = p2.first;

	string& p1al = p1.first;
	string& p2al = p2.second;

	unsigned int alSize = std::max(anch1.size(), anch2.size());
	alSize *= 1.2;

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

	//TODO - remove

	triAlignment[0] = inputSeqs->getDictionary()->translate(tr1, false);
	triAlignment[1] = inputSeqs->getDictionary()->translate(tr2, false);
	triAlignment[2] = inputSeqs->getDictionary()->translate(tr3, false);

	DUMP("Triplet alignment : ");
	DUMP(tr1);
	DUMP(tr2);
	DUMP(tr3);

}

array<vector<SequenceElement>, 3> TripletAligner::align(
		pair<string, string>& p1, pair<string, string>& p2)
{
	//FIXME - implement
	assembleFromPairs(p1,p2);
	return this->triAlignment;
}

array<vector<SequenceElement>, 3> TripletAligner::align(pair<vector<SequenceElement>, vector<SequenceElement> >& p1, pair<vector<SequenceElement>, vector<SequenceElement> >& p2)
{
	vector<SequenceElement> tr1;
	vector<SequenceElement> tr2;
	vector<SequenceElement> tr3;

	vector<SequenceElement>& anch1 = p1.second;
	vector<SequenceElement>& anch2 = p2.first;

	vector<SequenceElement>& p1al = p1.first;
	vector<SequenceElement>& p2al = p2.second;

	unsigned int alSize = std::max(anch1.size(), anch2.size());
	alSize *= 1.2;  // make it 20% longer to avoid reallocations

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
			tr1.push_back(p1al[ctr1]);
			tr2.push_back(anch1[ctr1]);
			tr3.push_back(p2al[ctr2]);
			ctr1 ++;
			ctr2 ++;
		}
		else if (anch1[ctr1].isIsGap())
		{
			tr1.push_back(p1al[ctr1]);
			//gap
			tr2.push_back(anch1[ctr1]);
			tr3.push_back(anch1[ctr1]);
			ctr1 ++;
		}
		else
		{
			tr1.push_back(anch2[ctr2]);
			//gap
			tr2.push_back(anch2[ctr2]);
			tr3.push_back(p2al[ctr2]);
			ctr2 ++;
		}
	}

	triAlignment[0] = tr1;
	triAlignment[1] = tr2;
	triAlignment[2] = tr3;

	return triAlignment;
}

array<vector<SequenceElement>, 3> TripletAligner::align(array<unsigned int, 3> triplet)
{
	s1 = triplet[0];
	s2 = triplet[1];
	s3 = triplet[2];

	string& seq1 = inputSeqs->getRawSequenceAt(s1);
	string& seq2 = inputSeqs->getRawSequenceAt(s2);
	string& seq3 = inputSeqs->getRawSequenceAt(s3);

	GotohAlgorithm* algo = new GotohAlgorithm(inputSeqs->getDictionary()->getAlphabetSize(),inputSeqs->getDictionary());
	algo->setDistance(distMat->getDistance(s1,s2));
	algo->setSequences(seq1,seq2);
	//scoring matrices etc....
	algo->run();
	auto firstPair = algo->getAlignment();

	algo->setDistance(distMat->getDistance(s2,s3));
	algo->setSequences(seq2,seq3);
	algo->run();
	auto secondPair = algo->getAlignment();

	assembleFromPairs(firstPair,secondPair);

	return this->triAlignment;

}

TripletAligner::TripletAligner(Sequences* iSeq, DistanceMatrix* dm) : inputSeqs(iSeq), distMat(dm)
{
	DEBUG("Starting TripletAligner");


}

} /* namespace EBC */


