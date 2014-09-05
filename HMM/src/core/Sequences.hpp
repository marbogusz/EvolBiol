/*
 * Sequences.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef SEQUENCES_H_
#define SEQUENCES_H_

#include <vector>
#include "core/Definitions.hpp"
#include "core/IParser.hpp"
#include "core/Dictionary.hpp"
#include "core/HmmException.hpp"
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"

namespace EBC
{


class Sequences
{
private:

	vector<string> rawSequences;
	vector<vector<SequenceElement> > translatedSequences;
	vector<std::pair<unsigned int, unsigned int> > pairs;

	vector<std::pair<unsigned int, unsigned int> >::iterator pairIterator;

	Dictionary* dict;
	unsigned int sequenceCount;
	double* observedFrequencies;

	bool fixedAlignment;

public:

	//Input from file or console
	Sequences(IParser*, Definitions::SequenceType, bool fixedAlignment=false) throw (HmmException&);

	~Sequences();

	//Return the dictionary for the input set
	Dictionary* getDictionary();

	//Get then in a dictionary order eg. T C A G, check the definitions in constants
	//FIXME - change to A T C G
	double* getElementFrequencies();

	//void getSequencePair(vector<SequenceElement> s1, vector<SequenceElement> s2 );
	vector<SequenceElement> getSequencesAt(int pos);

	unsigned int getPairCount()
	{
		unsigned int ct = translatedSequences.size();
		return (ct*(ct-1))/2;
	}

	inline unsigned int getSequenceCount()
	{
		return translatedSequences.size();
	}

	string& getRawSequenceAt(unsigned int pos);

	std::pair<unsigned int, unsigned int> getPairOfSequenceIndices(unsigned int idx)
	{
		return pairs[idx];
	}
private:

	void buildDictionary(Definitions::SequenceType);

};

} /* namespace EBC */



#endif /* SEQUENCES_H_ */
