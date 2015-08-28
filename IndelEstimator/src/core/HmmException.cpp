/*
 * ParseException.cpp
 *
 *  Created on: Sep 23, 2013
 *      Author: mbogusz
 */

#include "core/HmmException.hpp"
#include "core/Definitions.hpp"

namespace EBC {

HmmException::HmmException()
{
	// General message
	msg = "Parse Exception";
}

HmmException::~HmmException() throw()
{
}

HmmException::HmmException(string message) : msg(message)
{
	ERROR(message);
}

const char* HmmException::what() const throw ()
{
	return this->msg.c_str();
}

} /* namespace EBC */
