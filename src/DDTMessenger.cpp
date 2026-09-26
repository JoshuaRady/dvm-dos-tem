/***********************************************************************************************//**
 * @file DDTMessenger.cpp
 * \author Joshua M. Rady
 * Woodwell Climate Research Center
 * \date 2026
 *
 * @brief This derived class is used to connect the Fireweed messaging system into DVM-DOS-TEM's
 * Boost based system.  It inherits from the FWMessenger class and overrides the core functions that
 * determine the way in which the messages are recorded and fatal errors are handled.
 *
 **************************************************************************************************/

#include "./include/DDTMessenger.h"
#include "../include/TEMLogger.h"

#include <stdexcept>

extern src::severity_logger< severity_level > glg;

/** Post neutral log messages at the info level.
 *
 * @param message A message to log.
 */
void DDTMessenger::Log(const std::string& message) const
{
	BOOST_LOG_SEV(glg, info) << message;
}

/** Post a non-fatal warning.
 *
 * @param message A warning message.
 */
void DDTMessenger::Warning(const std::string& message) const
{
	BOOST_LOG_SEV(glg, warn) << message;
}

/** Post the passed message and shutdown.
 *
 * This is used in for fatal errors that can't be recovered from. The message is recorded and an 
 * exception it thrown to terminate the current grid cell without stopping simulation of other
 * cells.
 *
 * @param message An error message.
 */
void DDTMessenger::Stop(const std::string& message) const
{
	*errorStream << "Error: " << message << std::endl;//Temporary: Leave in the leaky messaging for now.
	//The code catching the exception should record the message but we do it manually to be safe.
	BOOST_LOG_SEV(glg, fatal) << message;
	throw std::runtime_error(message);
}
