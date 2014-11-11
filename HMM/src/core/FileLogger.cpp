/*
 * FileLogger.cpp

 *
 *  Created on: Nov 11, 2014
 *      Author: marcin
 */

#include "core/FileLogger.hpp"

namespace EBC
{

FileLogger FileLogger::instance;


FileLogger& FileLogger::getLogger()
{
	return instance;
}



}



