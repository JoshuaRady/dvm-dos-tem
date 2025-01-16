/**
 * @file FW_Interface.h
 * \author Joshua M. Rady
 * Woodwell Climate Reseach Center
 * \date 2024
 *
 * @brief This header file defines functions to connect the DVM-DOS-TEM WildFire class to the
 * revised wildfire model components.
 *
 * This is currently a very rough draft, possibly a placeholder.  Everything including the name
 * should be expected to change.
 */

#ifndef FW_INTERFACE_H
#define FW_INTERFACE_H

#include "FireweedRAFireSpread.h"

int GetMatchingFuelModel(int cmt);
void GetDeadFuelSizeDistribution(const FuelModel& fm, std::vector <double>& distribSAVs,
                                 std::vector <double>& distribWts);
bool IsShrub(const int cmtNumber, const int pftIdx);
void CalculateFuelBedDepth(FuelModel& fm, bool dynamic = true);
void SimulateSurfaceCombustion(const FuelModel& fm, SpreadCalcs raData, double tempAir, double windSpeed);
double SimulateGroundFire();

#endif //FW_INTERFACE_H
