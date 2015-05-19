/*
 * Dictionary.cpp
 *
 *  Created on: Jan 14, 2014
 *      Author: root
 */

#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include "core/HmmException.hpp"
#include <algorithm>

namespace EBC
{


//FIXME  - U and T equivalence!!!

void Dictionary::setAlphabet(const string dict[], unsigned short size)
{
	unsigned short i;

	for (i=0; i<=size;i++){
		alphabet.push_back(dict[i]);
	}

	for(unsigned short i=0; i<=size; i++)
	{
		//this->alphabet.push_back(string(1,dict[i]));
		this->translator.insert(std::make_pair(alphabet[i],new SequenceElement(i==gapId, i, NULL, alphabet[i])));
	}

	//alphabet size does not include gap e.g. size is 4 for nucleotides
	this->alphabetSize = size;
}

void Dictionary::outputAlphabet()
{
	cout << "Model dictionary: " << endl;
	for(auto sym : alphabet){
		cout << sym << ",";
	}
}

string& Dictionary::getSymbolAt(unsigned char i)
{
	return (alphabet[i]);
}

SequenceElement* Dictionary::getSequenceElement(string& symbol)
{
	return translator[symbol];
}

vector<SequenceElement*>* Dictionary::translate(string& sequence, bool disregardIndels)
{

	vector<SequenceElement*> *translatedVector = new vector<SequenceElement*>(sequence.size());
	unsigned short currentEl;

//	DEBUG("Transled: ");
	unsigned int pos = 0;
	for(string::iterator it = sequence.begin(); it < sequence.end(); it++)
	{
		string cstrg = string(1, *it);
		(*translatedVector)[pos] = getSequenceElement(cstrg);
		pos++;
		//currentEl = getSymbolIndex(*it);
		//if (currentEl == alphabetSize && disregardIndels)
		//	continue;
//		DEBUGN(currentEl);
		//translatedVector.push_back(SequenceElement((currentEl == alphabetSize), currentEl,NULL, getSymbolAt(currentEl)));
		//translatedVector.push_back(SequenceElement((currentEl == alphabetSize), currentEl,NULL, getSymbolAt(currentEl)));
	}
//	DEBUGN(std::endl);
	return translatedVector;

}

const string Dictionary::nucleotides[] = {"T", "C", "A", "G", "-"};

//FIXME - T and U equivalence needs to be addressed
//const string Dictionary::UracilSymbol = "U";
//const unsigned int Dictionary::UracilId = 0;

//const char Dictionary::nucleotides[] = {'A', 'C', 'G', 'T'};
const string Dictionary::aminoacids[] = {"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-"};
									//	  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20

const string Dictionary::codons[] = {"TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG","TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
									 "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG","CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
									 "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG","AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
									 "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG","GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG", "---"};


//aminoacid Id to codon conversion
const int Dictionary::geneticCode[] = {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,-1,17,
									   10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
                                        9, 9, 9,-1,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
                                       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7 };


const string Dictionary::gapChar = "-";

unsigned short Dictionary::getAlphabetSize()
{
	return this->alphabetSize;
}


NucleotideDictionary::NucleotideDictionary()
{
	gapId = 4;
	this->setAlphabet(Dictionary::nucleotides,4);

	this->translator.insert(std::make_pair("U",new SequenceElement(false, 0, NULL, "U")));

}


AminoacidDictionary::AminoacidDictionary()
{
	gapId = 20;
	this->setAlphabet(Dictionary::aminoacids,20);

}

CodonDictionary::CodonDictionary()
{
	this->setGeneticCode(Dictionary::geneticCode);
	unsigned short codonNo = 64;
	//this->setAlphabet(Dictionary::codons,64);

	unsigned short i,j;
	unsigned short actualNo = 0;

	for (i = 0; i <= codonNo ;i++){
		if (genCode[i] > 0){
			actualNo++;
		}
		alphabet.push_back(Dictionary::codons[i]);
	}
	//push the gap as well
	gapId = codonNo;
	alphabetSize = actualNo;

	i=0;
	j=0;
	short negIdx = -1;
	for(auto cod : alphabet)
	{
			//this->alphabet.push_back(string(1,dict[i]));

		this->translator.insert(std::make_pair(cod,new SequenceElement(i==gapId, genCode[i] < 0 ? negIdx-- : j++, NULL, alphabet[i])));
		i++;
	}
		//alphabet size does not include gap e.g. size is 4 for nucleotides
}

void CodonDictionary::setGeneticCode(const int gc[])
{
	std::copy(gc,gc + this->alphabetSize, this->genCode);
}

int CodonDictionary::getAminoacidId(unsigned int codonId)
{
	return genCode[codonId];
}

bool CodonDictionary::isPurine(char base)
{
	return (base == 'A' || base == 'G');
}

bool CodonDictionary::isPyramidine(char base)
{
	return (base == 'T' || base == 'C');
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

//	DEBUG("Transled: ");


	unsigned int pos = 0;
	unsigned int strpos = 0;


	while(strpos < sequence.size()){
		string cstrg = sequence.substr(strpos,3);
		(*translatedVector)[pos] = getSequenceElement(cstrg);
		strpos +=3;
		pos++;

	}
//	DEBUGN(std::endl);
	return translatedVector;

}

}//Namespace definition
