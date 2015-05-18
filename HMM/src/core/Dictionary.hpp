/*
 * Dictionary.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef DICTIONARY_H_
#define DICTIONARY_H_

#include <vector>
#include <string>
#include <map>
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"

using namespace std;

namespace EBC
{
	class Dictionary
	{
	protected:
		unsigned short alphabetSize;
		unsigned char gapId;
		vector<string> alphabet;
		map<string,SequenceElement*> translator;

	public:

		static const string nucleotides[5];
		static const string aminoacids[21];
		static const string codons[65];
		static const int geneticCode[64];
		static const string gapChar;

		virtual vector<SequenceElement*>* translate(string &sequence, bool disregardIndels = false);

		virtual unsigned short getAlphabetSize();

		SequenceElement* getSequenceElement(string& symbol);

		virtual string& getSymbolAt(unsigned char i);

		virtual void outputAlphabet();

		inline unsigned char getGapID()
		{
			return gapId;
		}

	protected:
		virtual void setAlphabet(const string alphabet[], unsigned short size);

	};

	class NucleotideDictionary : public Dictionary
	{
	public:
		NucleotideDictionary();
	};

	class AminoacidDictionary : public Dictionary
	{
	public:
		AminoacidDictionary();
	};

	class CodonDictionary : public Dictionary
	{
	protected:
		int genCode[64];

		bool isPurine(char base);
		bool isPyramidine(char base);

	public:
		CodonDictionary();

		vector<SequenceElement*>* translate(string &sequence, bool disregardIndels = false);

		void setGeneticCode(const int[]);

		int getAminoacidId(unsigned int codonId);

		unsigned int getNumberOfDifferentPositions(unsigned int codon1, unsigned int codon2,
				bool& synonymous, bool& transition);
	};
}


#endif /* DICTIONARY_H_ */
