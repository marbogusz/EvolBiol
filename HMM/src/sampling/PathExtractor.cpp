/*
 * PathExtractor.cpp
 *
 *  Created on: Feb 11, 2015
 *      Author: marcin
 */

#include <sampling/PathExtractor.hpp>

namespace EBC {

PathExtractor::PathExtractor(string sq1, string sq2, vector<SequenceElement*>* vp1, vector<SequenceElement*>* vp2) :
		s1(sq1), s2(sq2), v1(*vp1), v2(*vp2), l1(s1.size()-1), l2(s2.size()-1) {
	// TODO Auto-generated constructor stub

}
double PathExtractor::getSumLnl(EvolutionaryPairHMM* phmm, Maths* mt){
	double res = Definitions::minMatrixLikelihood;
	cerr << "Path extractor " << vPaths.size() << " paths" << endl;
	for (auto pr : vPaths)
	{
		res = mt->logSum(res, phmm->getAlignmentLikelihood(&(pr.first), &(pr.second)));
	}
	return res;
}

void PathExtractor::genAllSeqs(vector<SequenceElement*> c1, vector<SequenceElement*> c2, SequenceElement* gap, unsigned int i1, unsigned int i2, unsigned int s)
{
	//s == 0 M
	//s == 1 I
	//s == 2 D
	//s == 3 - fresh start
	unsigned int id1, id2;

	if (s == 0 ){
		//M
		c1.push_back(v1[i1++]);
		c2.push_back(v2[i2++]);
	}
	else if (s == 1){
		//I
		c1.push_back(v1[i1++]);
		c2.push_back(gap);

	}
	else if (s == 2){
		//D
		c1.push_back(gap);;
		c2.push_back(v2[i2++]);
	}

	if (i1 > l1 && i2 > l2){
		vPaths.push_back(make_pair(c1,c2));
		return;
	}
	else{
		if(i1 > l1){
			this->genAllSeqs(c1,c2, gap, i1, i2,2);
		}
		else if(i2 > l2){
			this->genAllSeqs(c1,c2, gap, i1, i2,1);
		}
		else{
			this->genAllSeqs(c1,c2, gap, i1, i2,0);
			this->genAllSeqs(c1,c2, gap, i1, i2,1);
			this->genAllSeqs(c1,c2, gap, i1, i2,2);
		}
	}
}

void PathExtractor::genAllSeqs(string c1, string c2, unsigned int i1, unsigned int i2, unsigned int s)
{
	//s == 0 M
	//s == 1 I
	//s == 2 D
	//s == 3 - fresh start
	unsigned int id1, id2;
/*
	if(i1 > l1){
		//fill the rest and return
		while(i2 <= l2){
			c1 += '-';
			c2 += s2[i2];
			i2++;
		}
		paths.push_back(make_pair(c1,c2));
		cerr << c1 << endl << c2 << endl << endl;
		return;
	}
	if(i2 > l2){
		while(i1 <= l1){
			c2 += '-';
			c1 += s1[i1];
			i1++;
		}
		paths.push_back(make_pair(c1,c2));
		cerr << c1 << endl << c2 << endl << endl;
		return;
	}
*/
	if (s == 0 ){
		//M
		c1 += s1[i1++];
		c2 += s2[i2++];
	}
	else if (s == 1){
		//I
		c1 += s1[i1++];
		c2 += "-";

	}
	else if (s == 2){
		//D
		c1 += '-';
		c2 += s2[i2++];
	}

	if (i1 > l1 && i2 > l2){
		sPaths.push_back(make_pair(c1,c2));
		cerr << c1 << endl << c2 << endl << endl;
		return;
	}
	else{
		if(i1 > l1){
			this->genAllSeqs(c1,c2, i1, i2,2);
		}
		else if(i2 > l2){
			this->genAllSeqs(c1,c2, i1, i2,1);
		}
		else{
			this->genAllSeqs(c1,c2, i1, i2,0);
			this->genAllSeqs(c1,c2, i1, i2,1);
			this->genAllSeqs(c1,c2, i1, i2,2);
		}
	}
}

PathExtractor::~PathExtractor() {
	// TODO Auto-generated destructor stub
}

} /* namespace EBC */
