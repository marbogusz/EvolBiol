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
 * InputParser.hpp
 *
 *  Created on: Sep 20, 2013
 *      Author: mbogusz
 */


#ifndef INPUTPARSER_H_
#define INPUTPARSER_H_

#include <iostream>
#include <string>
#include "core/HmmException.hpp"
#include "core/IParser.hpp"

namespace EBC {

class TextInputParser : public IParser {

protected:



	std::string sequenceA;
	std::string sequenceB;
	bool dbg;			//debug mode
	bool moreData;

public:

	TextInputParser(char* s1, char* s2, bool debug) throw (HmmException&);
	bool hasMoreData();
	bool getSequencePair(string& seqA, string& seqB);
	~TextInputParser();

};

} /* namespace EBC */
#endif /* INPUTPARSER_H_ */
