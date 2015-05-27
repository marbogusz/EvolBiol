/*
 * Dictionary.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef AADICTIONARY_H_
#define AADICTIONARY_H_

#include <vector>
#include <string>
#include <map>
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"
#include "core/Dictionary.hpp"

using namespace std;

namespace EBC
{
	class AminoacidDictionary : public Dictionary
	{
	public:
		AminoacidDictionary(Definitions::FrequencyScheme fs);
	};
}


#endif /* DICTIONARY_H_ */
