/*
 * Sequences.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "core/Sequences.hpp"

namespace EBC
{

Sequences::Sequences(IParser* iParser,Definitions::SequenceType st, bool fa) throw (HmmException&)
{
	//use the file parser to get sequences and build the dictionary
	fixedAlignment = fa;
	unsigned int size = iParser->getSequenceCount();
	if (size <= 0)
	{
		throw HmmException("No sequences found. Quitting");
	}

	this->buildDictionary(st);

	this->sequenceCount = size;

	pairs.reserve(this->getPairCount());

	for(unsigned int i=0; i< size;i++)
		for(unsigned int j=i+1; j<size;j++)
			pairs.push_back(std::make_pair(i,j));

		pairIterator = pairs.begin();

	while(size > 0)
	{
		this->rawSequences.push_back(iParser->getNextSequence());
		this->translatedSequences.push_back(dict->translate(*(rawSequences.end()-1),fixedAlignment==false));
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

//FIXME - RETURN A REFERENCE!
string& Sequences::getRawSequenceAt(unsigned int pos)
{
	return this->rawSequences[pos];
}

double* Sequences::getElementFrequencies()
{
	this->observedFrequencies = new double[dict->getAlphabetSize()];
	int i;

	for (i=0; i<dict->getAlphabetSize(); i++)
		this->observedFrequencies[i] = 0;

	unsigned int count =0;

	for (vector<vector<SequenceElement> >::iterator it1 = translatedSequences.begin() ; it1 != translatedSequences.end(); ++it1)
	{
		for(vector<SequenceElement>::iterator it2 = it1->begin(); it2 != it1->end(); ++it2)
		{
			if (!(it2->isIsGap()))
			{
				count++;
				observedFrequencies[it2->getMatrixIndex()]++;
			}
		}
	}

	for (i=0; i<dict->getAlphabetSize(); i++)
		this->observedFrequencies[i]/=count;


	DEBUGV(observedFrequencies,4);
	return observedFrequencies;
}

double* Sequences::getElementFrequencies(array<unsigned int, 3> triplet)
{
		this->observedFrequencies = new double[dict->getAlphabetSize()];
		int i;

		for (i=0; i<dict->getAlphabetSize(); i++)
			this->observedFrequencies[i] = 0;

		unsigned int count =0;

		for (auto it1 : triplet)
		{
			for(vector<SequenceElement>::iterator it2 = translatedSequences[it1].begin(); it2 != translatedSequences[it1].end(); ++it2)
			{
				if (!(it2->isIsGap()))
				{
					count++;
					observedFrequencies[it2->getMatrixIndex()]++;
				}
			}
		}

		for (i=0; i<dict->getAlphabetSize(); i++)
			this->observedFrequencies[i]/=count;

		DEBUGV(observedFrequencies,4);
		return observedFrequencies;
}

void Sequences::buildDictionary(Definitions::SequenceType st)
{
	//TODO differentiate based on the contents of the file!!!
	switch(st)
	{
	case (Definitions::SequenceType::Aminoacid):
		dict = new AminoacidDictionary();
		break;

	case (Definitions::SequenceType::Nucleotide):
		dict = new NucleotideDictionary();
		break;
	}
}



} /* namespace EBC */

