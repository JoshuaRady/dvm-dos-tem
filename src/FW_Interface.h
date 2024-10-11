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

void RevisedFire(const Cohort& thisCohort, int monthIndex);

int GetMatchingFuelModel(int cmt);
void CohortStatesToFuelLoading(const Cohort& thisCohort, FuelModel& fm, bool treatMossAsDead);
void GetDeadFuelSizeDistribution(const FuelModel& fm, std::vector <double>& distribSAVs,
                                 std::vector <double>& distribWts)
bool IsShrub(const Cohort& thisCohort, int pftIdx);

double GetMidflameWindSpeed();

std::vector <double> CalculateFuelMoisture(const Cohort& thisCohort, int monthIndex);

#endif //FW_INTERFACE_H
