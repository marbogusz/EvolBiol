/*
 * InputParser.cpp
 *
 *  Created on: Sep 20, 2013
 *      Author: mbogusz
 */

#include "core/TextInputParser.hpp"

using namespace std;

namespace EBC {

TextInputParser::TextInputParser(char* s1, char* s2, bool debug) throw (ProgramException&) : sequenceA (s1), sequenceB (s2)
{
	this->dbg = debug;
	this->moreData = true;

	if(this->sequenceA.length() == 0 || this->sequenceB.length() == 0)
	{
		throw ProgramException("Zero length sequence. Quitting");
	}

	if (this->dbg)
	{
		cerr << this->sequenceA << endl;
		cerr << this->sequenceB << endl;
	}
}

bool TextInputParser::getSequencePair(string& seqA, string& seqB)
{
	seqA = this->sequenceA;
	seqB = this->sequenceB;
	this->moreData = false;
	return true;
}

bool TextInputParser::hasMoreData()
{
	return this->moreData;
}

TextInputParser::~TextInputParser()
{
}

} /* namespace EBC */


