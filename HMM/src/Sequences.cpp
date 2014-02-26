/*
 * Sequences.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "Sequences.hpp"

namespace EBC
{

Sequences::Sequences(IParser* iParser) throw (ProgramException&)
{
	//use the file parser to get sequences and build the dictionary
	unsigned int size = iParser->getSequenceCount();
	if (size <= 0)
	{
		throw ProgramException("No sequences found. Quitting");
	}

	this->buildDictionary();

	this->sequenceCount = size;

	while(size > 0)
	{
		this->rawSequences.push_back(iParser->getNextSequence());
		this->translatedSequences.push_back(dict->translate(*(rawSequences.end()-1),false));
		size--;
	}


}

/*void Sequences::getSequencePair(vector<SequenceElement> s1,
		vector<SequenceElement> s2) {
	//FIXME - a mockup, only get the first pair
	//vector<vector<SequenceElement> >::iterator it = translatedSequences.begin();
	//TODO - consider friend with Algorithm
	s1 = translatedSequences[0];
	s2 = translatedSequences[1];
}*/

vector<SequenceElement> Sequences::getSequencesAt(int pos) {
		return translatedSequences[pos];
}

Sequences::~Sequences()
{
	delete dict;
	delete observedFrequencies;
}

Dictionary* Sequences::getDictionary()
{
	return dict;
}

unsigned int Sequences::getSequenceCount()
{

}

//Assuming the longest size
unsigned int Sequences::getSequenceSize()
{

}

double* Sequences::getElementFrequencies()
{
	//TODO implement
	//FIXME !!!!
	this->observedFrequencies = new double[dict->getAlphabetSize()];

	observedFrequencies[0] = 0.26089;
	observedFrequencies[1] = 0.32737;
	observedFrequencies[2] = 0.30782;
	observedFrequencies[3] = 0.10391;


	//observedFrequencies[0] = 0.25;
	//observedFrequencies[1] = 0.25;
	//observedFrequencies[2] = 0.25;
	//observedFrequencies[3] = 0.25;

	return observedFrequencies;
}

void Sequences::buildDictionary()
{
	//TODO differentiate based on the contents of the file!!!
	dict = new NucleotideDictionary();
}



} /* namespace EBC */

