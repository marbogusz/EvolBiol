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
	observedFrequencies = NULL;

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

	this->rawSequences = iParser->getSequences();
	this->sequenceNames = iParser->getNames();

	for (auto it = rawSequences->begin(); it != rawSequences->end(); it++)
		this->translatedSequences.push_back(dict->translate(*it,fixedAlignment==false));

}

vector<SequenceElement*>* Sequences::getSequencesAt(unsigned int pos){
		return translatedSequences[pos];
}

Sequences::~Sequences()
{
	delete dict;
	delete[] observedFrequencies;
}

Dictionary* Sequences::getDictionary()
{
	return dict;
}

//FIXME - RETURN A REFERENCE!
string& Sequences::getRawSequenceAt(unsigned int pos)
{
	return (*rawSequences)[pos];
}

string& Sequences::getSequenceName(unsigned int pos)
{
	return (*sequenceNames)[pos];
}

void Sequences::calculateObservedFrequencies()
{
	this->observedFrequencies = new double[dict->getAlphabetSize()];

	unsigned int i;

	for (i=0; i<dict->getAlphabetSize(); i++)
		this->observedFrequencies[i] = 0;

	unsigned int count =0;

	for (auto it1 = translatedSequences.begin() ; it1 != translatedSequences.end(); ++it1)
	{
		for(auto it2 = (*it1)->begin(); it2 != (*it1)->end(); ++it2)
		{
			if (!((*it2)->isIsGap()))
			{
				count++;
				observedFrequencies[(*it2)->getMatrixIndex()]++;
			}
		}
	}

	for (i=0; i < dict->getAlphabetSize(); i++){
		//FIXME -change  asap!!!!
		this->observedFrequencies[i] = 1.0/61;// /= count;

	}
}

double* Sequences::getElementFrequencies()
{
	if(observedFrequencies == NULL)
		calculateObservedFrequencies();
	//DEBUGV(observedFrequencies,4);
	return observedFrequencies;
}

double* Sequences::getElementFrequencies(array<unsigned int, 3>& triplet)
{

	if(observedFrequencies == NULL)
		this->observedFrequencies = new double[dict->getAlphabetSize()];
	int i;

	for (i=0; i<dict->getAlphabetSize(); i++)
		this->observedFrequencies[i] = 0;

	unsigned int count =0;

	for (auto it1 : triplet)
	{
		for(auto it2 = (translatedSequences[it1])->begin(); it2 != (translatedSequences[it1])->end(); ++it2)
		{
			if (!((*it2)->isIsGap()))
			{
				count++;
				observedFrequencies[(*it2)->getMatrixIndex()]++;
			}
		}
	}

	for (i=0; i < dict->getAlphabetSize(); i++)
		this->observedFrequencies[i] /= count;

	//DEBUGV(observedFrequencies,4);
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

	case (Definitions::SequenceType::Codon):
		dict = new CodonDictionary();
		break;
	}
}



} /* namespace EBC */


