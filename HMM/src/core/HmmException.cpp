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
	// TODO Auto-generated constructor stub
	// General message
	msg = "Parse Exception";
}

HmmException::~HmmException() throw()
{
	// TODO Auto-generated destructor stub
}

HmmException::HmmException(string message) : msg(message)
{
	FileLogger::ErrorLogger() <<  message << '\n';
}

const char* HmmException::what() const throw ()
{
	return this->msg.c_str();
}

} /* namespace EBC */
