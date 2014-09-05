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
