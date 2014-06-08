/*
 * Dictionary.hpppp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef DICTIONARY_H_
#define DICTIONARY_H_

#include <vector>
#include <string>
#include <map>
#include "SequenceElement.hpp"
#include "Definitions.hpp"

using namespace std;

namespace EBC
{
	class Dictionary
	{
	protected:
		unsigned short alphabetSize;
		vector<string> alphabet;
		map<string,short> translator;

	public:

		static const char nucleotides[4];
		static const char aminoacids[20];

		virtual short getSymbolIndex(string &symbol);
		virtual short getSymbolIndex(char symbol);
		virtual vector<SequenceElement> translate(string &sequence, bool disregardIndels);
		virtual unsigned short getAlphabetSize();

		virtual string getSymbolAt(short i);

		virtual void outputAlphabet();

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
