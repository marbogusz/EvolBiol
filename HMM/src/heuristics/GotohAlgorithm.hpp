/*
 * GotohAlgorithm.h
 *
 *  Created on: Sep 23, 2013
 *      Author: mbogusz
 */

#ifndef GOTOHALGORITHM_H_
#define GOTOHALGORITHM_H_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "heuristics/GotohScoringMatrix.hpp"

#include "heuristics/GotohSMatrix.hpp"
#include "heuristics/GotohHMatrix.hpp"
#include "heuristics/GotohVMatrix.hpp"

using namespace std;

namespace EBC
{

class GotohAlgorithm
{

protected:

	string seq_a;
	string seq_b;

	GotohScoringMatrix* scores;

	GotohHMatrix* H;
	GotohVMatrix* V;
	GotohSMatrix* S;


	unsigned int seqALen, seqBLen;
	unsigned int xSize, ySize;

	//void traceback();
	void processMatrices();


public:

	GotohAlgorithm(string& a, string& b);

	virtual ~GotohAlgorithm();

	void setSequences(string& a, string& b);

	void run();
};

} /* namespace EBC */



#endif /* GOTOHALGORITHM_H_ */
