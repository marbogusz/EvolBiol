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
		string alphabet;
		map<char,SequenceElement*> translator;

	public:

		static const char nucleotides[5];
		static const char aminoacids[21];
		static const char gapChar;

		//virtual short getSymbolIndex(string &symbol);
		virtual unsigned char getSymbolIndex(char symbol);
		virtual vector<SequenceElement*>* translate(string &sequence, bool disregardIndels = false);
		virtual unsigned short getAlphabetSize();

		SequenceElement* getSequenceElement(char symbol);

		virtual char getSymbolAt(unsigned char i);

		virtual void outputAlphabet();

		inline unsigned char getGapID()
		{
			return gapId;
		}

	protected:
		virtual void setAlphabet(char alphabet[], unsigned short size);

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
}


#endif /* DICTIONARY_H_ */
