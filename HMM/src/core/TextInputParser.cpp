//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

/*
 * InputParser.cpp
 *
 *  Created on: Sep 20, 2013
 *      Author: mbogusz
 */

#include "core/TextInputParser.hpp"

using namespace std;

namespace EBC {

TextInputParser::TextInputParser(char* s1, char* s2, bool debug) throw (HmmException&) : sequenceA (s1), sequenceB (s2)
{
	this->dbg = debug;
	this->moreData = true;

	if(this->sequenceA.length() == 0 || this->sequenceB.length() == 0)
	{
		throw HmmException("Zero length sequence. Quitting");
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


