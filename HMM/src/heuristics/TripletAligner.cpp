/*
 * TripletAligner.cpp
 *
 *  Created on: Sep 8, 2014
 *      Author: root
 */

#include <heuristics/TripletAligner.hpp>
#include <vector>
#include <array>

#include <sstream>

namespace EBC
{
/*
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
*/
array<vector<SequenceElement*>*, 3>* TripletAligner::align(pair<vector<SequenceElement*>*, vector<SequenceElement*>* >* p1, pair<vector<SequenceElement*>*, vector<SequenceElement*>* >* p2)
{
	stringstream sp1, sp2, sp3, sp4;


/*
	for (auto i1 : p1.first)
	{
		sp1 << i1.getSymbol();
	}
	for (auto i2 : p1.second)
	{
		sp2 << i2.getSymbol();
	}
	for (auto i3 : p2.first)
	{
		sp3 << i3.getSymbol();
	}
	for (auto i4 : p2.second)
	{
		sp4 << i4.getSymbol();
	}
	//cout << "TRIPLET ALIGNER " << endl;
	//cout << sp1.str() << endl;
	//cout << sp2.str() << endl;
	//cout << sp3.str() << endl;
	//cout << sp4.str() << endl<<endl;

	DUMP("TRIPLET ALIGNER ");
	DUMP(sp1.str());
	DUMP(sp2.str());
	DUMP(sp3.str());
	DUMP(sp4.str());
*/
	vector<SequenceElement*>* tr1 = new vector<SequenceElement*>();
	vector<SequenceElement*>* tr2 = new vector<SequenceElement*>();
	vector<SequenceElement*>* tr3 = new vector<SequenceElement*>();

	vector<SequenceElement*>* anch1 = p1->second;
	vector<SequenceElement*>* anch2 = p2->first;

	vector<SequenceElement*>* p1al = p1->first;
	vector<SequenceElement*>* p2al = p2->second;

	int c;
	/*
	for (c=0; c < p1al->size(); c++)
	{
		cerr << (*p1al)[c]->getSymbol();
	}
	cerr << endl;
	for (c=0; c < anch1->size(); c++)
	{
		cerr << (*anch1)[c]->getSymbol();
	}
	cerr << endl;
	for (c=0; c < anch2->size(); c++)
	{
		cerr << (*anch2)[c]->getSymbol();
	}
	cerr << endl;
	for (c=0; c < p2al->size(); c++)
	{
		cerr << (*p2al)[c]->getSymbol();
	}
	cerr << endl;
	cerr << anch1->size() << endl << anch2->size() << endl;
*/
	unsigned int alSize = std::max(anch1->size(), anch2->size());
	alSize *= 1.2;  // make it 20% longer to avoid reallocations

	tr1->reserve(alSize);
	tr3->reserve(alSize);
	tr2->reserve(alSize);

	unsigned int ctr1;
	unsigned int ctr2;

	ctr1 = ctr2 = 0;

	while(ctr1 < anch1->size() && ctr2 < anch2->size())
	{
		if(((*anch1)[ctr1]->getMatrixIndex()) == ((*anch2)[ctr2]->getMatrixIndex()))
		{
			tr1->push_back((*p1al)[ctr1]);
			tr2->push_back((*anch1)[ctr1]);
			tr3->push_back((*p2al)[ctr2]);
			//cerr << (*p1al)[ctr1]->getSymbol() << (*anch1)[ctr1]->getSymbol() << (*p2al)[ctr2]->getSymbol() << endl;
			ctr1 ++;
			ctr2 ++;
		}
		else if ((*anch1)[ctr1]->isIsGap())
		{
			tr1->push_back((*p1al)[ctr1]);
			//gap
			tr2->push_back((*anch1)[ctr1]);
			tr3->push_back((*anch1)[ctr1]);
			//cerr << (*p1al)[ctr1]->getSymbol() << (*anch1)[ctr1]->getSymbol() << (*anch1)[ctr1]->getSymbol() << endl;
			ctr1 ++;
		}
		else
		{
			tr1->push_back((*anch2)[ctr2]);
			//gap
			tr2->push_back((*anch2)[ctr2]);
			tr3->push_back((*p2al)[ctr2]);
			//cerr << (*anch2)[ctr2]->getSymbol() << (*anch2)[ctr2]->getSymbol() << (*p2al)[ctr2]->getSymbol() << endl;
			ctr2 ++;
		}
	}

	//triAlignment[0] = tr1;
	//triAlignment[1] = tr2;
	//triAlignment[2] = tr3;

	//int c;
	for (c=0; c < tr1->size(); c++)
	{
		sp1 << (*tr1)[c]->getSymbol();
	}
	for (c=0; c < tr2->size(); c++)
	{
		sp2 << (*tr2)[c]->getSymbol();
	}
	for (c=0; c < tr3->size(); c++)
	{
		sp3 << (*tr3)[c]->getSymbol();
	}
/*
	cerr << sp1.str() <<endl;
	cerr << sp2.str() <<endl;
	cerr << sp3.str() <<endl;
*/
	DUMP(sp1.str());
	DUMP(sp2.str());
	DUMP(sp3.str());
	return new array<vector<SequenceElement*>*, 3>({tr1,tr2,tr3});
}


TripletAligner::TripletAligner(Sequences* iSeq, DistanceMatrix* dm) : inputSeqs(iSeq), distMat(dm)
{
	DEBUG("Starting TripletAligner");


}

} /* namespace EBC */


