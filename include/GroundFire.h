/**
 * @file GroundFire.h
 * \author Joshua M. Rady
 * Woodwell Climate Research Center
 * \date 2025
 *
 * @brief This header file declares functions to simulate smoldering ground fire using a moderate
 * complexity model where combustion is predicted layer by layer.
 *
 */

#ifndef FW_GROUNDFIRE_H
#define FW_GROUNDFIRE_H

#include <vector>

#include "GFProfile.h"

double SoilHeatOfCombustion(const double inorganicPct, const double organicHOC = 21.39148);
double SoilHeatSink(const double moistureContentPct, const double T_int, const double T_ig, const double c_s);
std::vector <double> HeatDistributionLinear(const int numLayers);
double DominoGroundFire(GFProfile& soilCol, const double fireHeatInput = 0.0,//With no heat input nothing will burn.
                        const double heatLossFactor = 0.83, const double d_max = 5.0);

#endif //FW_GROUNDFIRE_H
