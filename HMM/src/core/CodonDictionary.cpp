/*
 * Dictionary.cpp
 *
 *  Created on: Jan 14, 2014
 *      Author: root
 */

#include "core/CodonDictionary.hpp"
#include "core/Definitions.hpp"
#include "core/HmmException.hpp"
#include <algorithm>

namespace EBC
{


//FIXME  - U and T equivalence!!!

CodonDictionary::CodonDictionary(Definitions::FrequencyScheme fs) : Dictionary(fs)
{
	DUMP("Codon dictionary: construction ");
	//this->setGeneticCode(Dictionary::geneticCode);
	unsigned short codonNo = 64;
	//this->setAlphabet(Dictionary::codons,64);

	unsigned short i,j;
	unsigned short actualNo = 0;

	for (i = 0; i < codonNo ;i++){
		if (Dictionary::geneticCode[i] >= 0){
			alphabet.push_back(Dictionary::codons[i]);
			reducedGenCode[actualNo] = Dictionary::geneticCode[i];
			actualNo++;
		}
	}
	//push the gap
	alphabet.push_back(Dictionary::codons[codonNo]);
	//push the gap as well
	gapId = actualNo;
	alphabetSize = actualNo;

	i=0;
	for(auto cod : alphabet)
	{
			//this->alphabet.push_back(string(1,dict[i]));
		this->translator.insert(std::make_pair(cod,new SequenceElement(i==gapId, i, NULL, alphabet[i])));
		i++;
	}

	elementFrequencies = new double[alphabetSize];

	for (unsigned int i=0; i<alphabetSize; i++)
			this->elementFrequencies[i] = 0;

	if (fs == Definitions::FrequencyScheme::Equal){
		setEqualFrequencies();
	}
}


int CodonDictionary::getAminoacidId(unsigned int codonId)
{
	return reducedGenCode[codonId];
}

bool CodonDictionary::isPurine(char base)
{
	return (base == 'A' || base == 'G');
}

bool CodonDictionary::isPyramidine(char base)
{
	return (base == 'T' || base == 'C');
}


void CodonDictionary::calculateFreqeuencies(){
	if (fScheme == Definitions::FrequencyScheme::Equal)
		return;
	else if (fScheme == Definitions::FrequencyScheme::Empirical){
		for (unsigned int i=0; i < alphabetSize; i++){
			elementFrequencies[i] /= simpleCount;
		}
	}
	else if (fScheme == Definitions::FrequencyScheme::F1x4){

	}

}

unsigned int CodonDictionary::getNumberOfDifferentPositions(unsigned int codon1,
		unsigned int codon2, bool& synonymous, bool& transition)
{
	synonymous = false;
	transition = true;

	string& c1 = alphabet[codon1];
	string& c2 = alphabet[codon2];

	//iterate over positions
	unsigned int counter = 0;
	int changedPos = -1;
	for(unsigned int i = 0; i < 3; i++)
	{
		if (c1[i] != c2[i]){
			counter++;
			changedPos = i;
		}
	}
	if(counter == 0 || counter > 1)
		return counter;

	if (getAminoacidId(codon1) ==  getAminoacidId(codon2))
		synonymous = true;

	bool c1Purine = isPurine(c1[changedPos]);
	bool c2Purine = isPurine(c2[changedPos]);

	if (c1Purine != c2Purine)
		transition = false;

	return counter;
}

vector<SequenceElement*>* CodonDictionary::translate(string& sequence, bool disregardIndels)
{
	if(sequence.size() % 3 != 0)
		throw HmmException("Sequence length should be divisible by 3 for codon models! Quitting");

	vector<SequenceElement*> *translatedVector = new vector<SequenceElement*>(sequence.size()/3);
	unsigned short currentEl;

	unsigned int pos = 0;
	unsigned int strpos = 0;
	SequenceElement* se;


	while(strpos < sequence.size()){
		string cstrg = sequence.substr(strpos,3);
		se = getSequenceElement(cstrg);
		(*translatedVector)[pos] = se;
		strpos +=3;
		pos++;
		if (fScheme == Definitions::FrequencyScheme::Empirical){
			simpleCount++;
			elementFrequencies[se->getMatrixIndex()] += 1;
		}
		//else if(fScheme == Definitions::FrequencyScheme::F1x4)

	}
//	DEBUGN(std::endl);
	return translatedVector;

}

}//Namespace definition
