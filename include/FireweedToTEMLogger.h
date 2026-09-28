/***********************************************************************************************//**
 * @file FireweedToTEMLogger.h
 * \author Joshua M. Rady
 * Woodwell Climate Research Center
 * \date 2026
 *
 * @brief This derived class is used to connect the Fireweed messaging system into DVM-DOS-TEM's
 * Boost based system.  It inherits from the FWMessenger class and overrides the core functions that
 * determine the way in which the messages are recorded and fatal errors are handled.
 *
 **************************************************************************************************/

#ifndef FIREWEEDTOTEMLOGGER_H
#define FIREWEEDTOTEMLOGGER_H

#include "FireweedMessaging.h"

class FWtoTEMMessenger : public FWMessenger {
  public:
    void Log(const std::string& message) const;
    void Warning(const std::string& message) const;
    void Stop(const std::string& message) const;
};

#endif //FIREWEEDTOTEMLOGGER_H
