/*
 * Sequences.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef SEQUENCES_H_
#define SEQUENCES_H_

#include <vector>
#include "Definitions.hpp"
#include "IParser.hpp"
#include "Dictionary.hpp"
#include "ParseException.hpp"
#include "SequenceElement.hpp"

namespace EBC
{


class Sequences
{
private:

	vector<string> rawSequences;
	vector<vector<SequenceElement> > translatedSequences;
	Dictionary* dict;
	unsigned int sequenceCount;
	double* observedFrequencies;

public:

	//Input from file or console
	Sequences(IParser*) throw (ProgramException&);

	~Sequences();

	//Return the dictionary for the input set
	Dictionary* getDictionary();

	//Get then in a dictionary order eg. T C A G, check the definitions in constants
	//FIXME - change to A T C G
	double* getElementFrequencies();

	//void getSequencePair(vector<SequenceElement> s1, vector<SequenceElement> s2 );
	vector<SequenceElement> getSequencesAt(int pos);

	string getRawSequenceAt(unsigned int pos);
private:

	void buildDictionary();

};

} /* namespace EBC */



#endif /* SEQUENCES_H_ */
