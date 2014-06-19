/*
 * FileParser.hpp
 *
 *  Created on: Oct 1, 2013
 *      Author: mbogusz
 */

#ifndef FILEPARSER_H_
#define FILEPARSER_H_

#include <fstream>
#include "HmmException.hpp"
#include "IParser.hpp"
#include <vector>

using namespace std;

namespace EBC
{

class FileParser : public IParser
{
private:

	string filename;
	ifstream infile;
	vector<string> sequences;
	vector<string>::iterator it;

public:

	FileParser(const char* filename);

	string getNextSequence();

	unsigned int getSequenceCount();

	string getSequenceAt(unsigned int);

	bool isDefinitionLine(string&);

	void trimWsChars(string&);

	virtual ~FileParser();
};

} /* namespace EBC */
#endif /* FILEPARSER_H_ */
