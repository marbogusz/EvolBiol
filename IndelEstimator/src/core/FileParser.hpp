/*
 * FileParser.hpp
 *
 *  Created on: Oct 1, 2013
 *      Author: mbogusz
 */

#ifndef FILEPARSER_H_
#define FILEPARSER_H_

#include <fstream>
#include "core/HmmException.hpp"
#include "core/IParser.hpp"
#include <vector>

using namespace std;

namespace EBC
{

class FileParser : public IParser
{
private:

	string filename;
	ifstream infile;
	vector<string>* sequences;
	vector<string>* names;
	vector<string>::iterator it;
	vector<string>::iterator itN;

public:

	FileParser(const char* filename);

	string getNextSequence();

	string getNextName();

	unsigned int getSequenceCount();

	string getSequenceAt(unsigned int);

	string getSequenceNameAt(unsigned int);

	bool isDefinitionLine(string&);

	string getSequenceName(string&);

	void trimWsChars(string&);

	virtual ~FileParser();

	inline vector<string>* getNames() {
		return names;
	}

	inline vector<string>* getSequences() {
		return sequences;
	}
};

} /* namespace EBC */
#endif /* FILEPARSER_H_ */
