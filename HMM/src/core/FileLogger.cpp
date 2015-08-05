/*
 * FileLogger.cpp

 *
 *  Created on: Nov 11, 2014
 *      Author: marcin
 */

#include "core/FileLogger.hpp"

namespace EBC
{

std::ofstream FileLogger::logFile;

FileLogger FileLogger::errL;
FileLogger FileLogger::wrnL;
FileLogger FileLogger::dbgL;
FileLogger FileLogger::dmpL;
FileLogger FileLogger::infL;
FileLogger FileLogger::clnL;


//FIXME Separate default logger?
FileLogger& FileLogger::Logger()
{
	return clnL;
}

FileLogger& FileLogger::DumpLogger()
{
	return dmpL;
}
FileLogger& FileLogger::DebugLogger()
{
	return dbgL;
}
FileLogger& FileLogger::ErrorLogger()
{

	return errL;
}
FileLogger& FileLogger::InfoLogger()
{
	return infL;
}
FileLogger& FileLogger::WarningLogger()
{
	return wrnL;
}

}



