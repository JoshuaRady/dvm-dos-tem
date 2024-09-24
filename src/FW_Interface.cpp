/**
 * @file FW_Interface.cpp
 * \author Joshua M. Rady
 * Woodwell Climate Reseach Center
 * \date 2024
 *
 * @brief This file contains functions to connect the DVM-DOS-TEM WildFire class to the revised
 * wildfire model components.
 *
 * This is currently a very rough draft, possibly a placeholder.  Everything including the name
 * should be expected to change.
 */

#include "FireweedRAFireSpread.h"
#include "FireweedDeadFuelMoistureFosberg.h"
#include "FireweedLiveFuelMoistureGSI.h"
#include "FireweedMetUtils.h"
#include "TEMUtilityFunctions.h"//For length_of_day().

/** Calculate wildfire behavior and effects using the modeled vegetation, fuels, and meteorology,
 *
The function takes the ecosystem state (cohort pointer? WildFire pointer) and meteorology...

Output:
The output should be stored in the cohort object...
The soil burn depth should be returned/recorded to take advantage of the existing code.

 */
void RevisedFire(WildFire* wf)//The name will definitely change.
{

  //Determine the fuel model matching the location's CMT:------------------
  
  //Get the CMT for the grid cell:
  //The CMT number is in the CohortData member of the cohort.  The WildFire object maintains a
  //pointer to it's parent cohort data but it is  private.  This code needs to be in the class,
  //a friend or have access to the cohort.
  int theCMTumber = wf->cd.cmttype;//Private data access!!!!!

  //Crosswalk from CMT to fuel model:
  fuelModelNumber = GetMatchingFuelModel(theCMTumber);

  //The fuel model table file needs to be added to the config file and be loaded:
  std::string fuelModelTablePath = "/Some/Path/Dropbox/StandardFuelModelTableFileName.csv";//Or tab delimited.
  FuelModel fm = GetFuelModelFromCSV(fuelModelTablePath, fuelModelNumber);

  //

  //Determine the surface fuels from the model vegetation and soil states and update the fuel
  //loadings from their default values:
  CohortStatesToLoading(wf->cd, fm);//Private data access!!!!!

  //Do we calculate the fuel bed depth or is it fixed?
  //We could calculate it from af fixed fuel fuel bed density as FATES does.
  //fm.delta = X;
  
  
  //Gather weather conditions:
  double tempAir = wf->edall.d_atms.ta;//Daily air temp (at surface).  Private data access!!!!!
  //Curten (daily if not hourly) (relative) humidity is needed.  A recent time history of
  //double humidity = ?????
  
  double windSpeed = GetMidflameWindSpeed();
  
  //The wildfire object contains the slope:
  //It is also in the CohortData object.
  double slope = wf->slope;//What at the units?  May have to convert to fractional slope.
  //Shortwave radiation may be needed for the moisture calculations.


  //Fuel moisture must be calculated:-------------
  //Note: It is better to calculate fuel moisture after calculating fuel loadings since that process
  //might change the fuel sizes.
  std::vector <double> M_f_ij = CalculateFuelMoisture(fm, tempAir, slope);

  //Add the moisture to the fuel model possibly computing dynamic fuel moisture:
  if (UseDynamicFuelMoisture)//Add switch for dynamic moisture!!!!!
  {
    fm.CalculateDynamicFuelCuring(M_f_ij);
  }
  else
  {
    fm.SetFuelMoisture(M_f_ij);
  }


  //Feed fuels and weather conditions into surface fire models:

  //First to Rothermel & Albini.  We use the calculations to get some component values.
  //This interface is under development.  It takes a fuel model (and attendant data) and returns the
  //calculation details.
  //M_f_ij does not need to be included if it is adde to the fuel model object above.
  SpreadCalcs raData = SpreadCalcsRothermelAlbini_Het(fm,
                                                      windSpeed,//U
                                                      slope)//,//slopeSteepness
                                                      //M_f_ij,//fuel moisture
                                                      //false,//useWindLimit
                                                      //FALSE);//debug


  //Use the fire energy flux into the soil along with the fuel model as input into Burnup:
  //There are additional fuel and environmental parameters needed here!!!!!
  //raData.XXXXX



  //Use the energy flux from the aboveground fire into the soil surface (RA + Burnup + crown) and
  //the fire air temp (Burnup fire environmental temperature?) as input to the ground fire model:
  //...
  
  //Also pass in the soil profile conditions...
  
  //SimulateGroundFire();
  //wf->fd->fire_soid.burnthick = burnDepth;
  //Include messaging of original code?????

}

/** Determine the fuel model (number) matching a given community type:
 *
 * Each CMT has [will have] a predetermined fuel model assigned to it via a parameter file.
 * It needs to be determined if this will be the CMT parmameter file or an additional lookup table.
 */
//FuelModel GetMatchingFuelModel(int cmt)
int GetMatchingFuelModel(int cmt)//Or could return fuel model code.
{
  //Get the number of the fuel model from the crosswalk in the parameter files.
  //This crosswalk needs to be made!!!!!
  int fuelModelNumber = 1;//Temporarily hardwired
  
  
  //If no match either throw an error or warn and return a default fuel model.
  
  //Get the fuel model data:
  //fuelModel = GetFuelModel(fuelModelNumber);
  //return fuelModel;
  return fuelModelNumber;
}

/** Get the wind speed at midflame height in m/min.
 *
 * The Rothermel Albini spread model takes midflame wind speed, something we are unlikely to ever
 * have.  If fact the midflame height is not even known prior to running the model, which results
 * in an unresolvable circularity.  In practice no one seems to worry about this much.  An estimate
 * of near surface wind speed is what we need.  As a first pass we are going to approximate it as
 * ~2m wind speed.  The input may not match this so we need to compute our best estimate.
 *
 We need a sub-daily value ~ daily value?
 * The code has to do the following:
 * - Get the daily windspeed from the host model
 *   Note: Windspeed is not yet available in DVM-DOS-TEM, but it should be soon.
 * - Calculate the windspeed at ~2m (or some other value?) if the provided wind speed is at another
 *   height.
 * - Possibly: Estimate a sub-daily from a daily mean value.  The fire may occur on timescale where
 *   it occurs at a specific time of day.  In this case it would be good to estimate what the wind
 *   speed was at this time.  This may be difficult as wind can be highly variable and may not be
 *   well predicted by diurnal cycles.
 * - Convert to m/min if needed.
 */
double GetMidflameWindSpeed()//Could pass in the desired height or time of day?
{
	//Temporary stub, return an arbitrary value!!!!!:
	//Summer average wind speeds are ~6 mph in Fairbanks Alaska.
	//6 * ftPerMi / 60 / ftPerM = 160.9344
	return 160.9344;
}

/** Calculate fuel moisture based on recent weather:
*
* This is a fire science dependent function and as such should probably be moved into the FW library.
* Here the emphasis is on determining the available inputs we can pass in.
* An alternative to calculating this at the time of fire would be to make fuel moisture a constantly
* calculated state.  This would moke more sense if litter existed as a distinct stock as well.

The 1 and 10 hour values could probably be estimated from a daily value, though we have to
assume a time of day.  We can superimpose a diurnal cycle.  The longer equilibrating fuels will
need some way of estimating the humidity value for the recent past.
Days since last rain is available in wf->edall.d_atms.

The Nelson and Fosberg models should be the first to consider.
These models calculate moistures for ideal hour classes, which have matching diameters.  We could
better inform the moistures with the actual SAV's diameters.
Note: consider adding diameteres to the fuel model class!!!!!

Inputs:
- Humidity, temp, recent precip?

 * @param tempAir The air temperature (degrees C).
 * @param slope				Units??????


 * @returns M_f_ij The fuel moisture for all fuel classes.  This is not returned in the fuel model
 * passed in because we don't know if curing is being applied.
 */
// FuelMoisture CalculateDeadFuelMoisture()
// {
//   FuelMoisture dfmc;//FMC = dead fuel moisture content.
//   
//   //...
//   
//   return fmc;
// }
std::vector <double> CalculateFuelMoisture(FuelModel& fm, double tempAir, double slope)
{
  //thisCohort = the cohort we are in!!!!!
  
  std::vector <double> M_f_ij(fm.numClasses, 0);//Return value.

  //Dead fuel moisture:
  //

  int hourOfDay = 15;//

  int monthOfYearIndex = X;//0 based			Get!!!!!!
  int monthOfYear = monthOfYearIndex + 1;
  int dayOfMonthIndex = Y;//0 based...		Get!!!!!!

  //Calculate relative humidity:
  p_hPa = Z;//Get the atmospheric pressure...

  //Partial pressure of water vapor...
  P = thisCohort->climate.vapo_d[dayOfMonthIndex];//My reading is that vapo_d is a vector of daily values for the whole year.
  
  double P_s = SaturationVaporPressureBuck(tempAir, p_hPa);//Or Climate.cpp calculate_saturated_vapor_pressure()? 
  //Appears be available as thisCohort->climate.svp_d[]?

  double rh = RHfromVP(P, P_s);//Calculate relative humidity.

  //Get the aspect...

  double oneHrFM = FosbergNWCG_1HrFM(std::string tableA_Path, std::string tableB_Path, std::string tableC_Path,
                         std::string tableD_Path,
                         tempAir, rh, 
                         monthOfYear, hourOfDay,
                         double slopePct, char aspectCardinal, bool shaded = false,
                         char elevation = 'L', UnitsType units = US);

  double tenHrFM = NWCG_10hrFM(double oneHrFM);
  double hundredHrFM = NWCG_100hrFM(double oneHrFM);

  //Live fuel moisture:
  
  //Minimum and maximum daily temperature are not currently available but will be soon.  We use a
  //stub for now.
  double tempCMin = GetDailyMinimumTemperature(tempAir);

  //Calculate VPD from available mean daily weather values:
  //double p_hPa = ????;
  double vpdPa = VPDfromRHBuck(tempAir, rh, p_hPa);// Change units!!!!!
  //Or use Climate.cpp calculate_vpd() and calculate_saturated_vapor_pressure().
  //Or may be available as thisCohort->climate.vpd_d[]?

  //The day length can be obtained from the length_of_day() routine in TEMUtilityFunctions.h, which
  //is in hours.
  int dayOfYear = temutil::day_of_year(monthOfYearIndex, dayOfMonthIndex);
  //double dayLength = length_of_day();
  double dayLengthSec = temutil::length_of_day(cohort->lat, dayOfYear) * 60 * 60;//Need the cohort!!!!!
  
  double gsi = GrowingSeasonIndex(tempCMin, vpdPa, dayLengthSec);
  double HerbaceousLiveFuelMoisture(gsi);
  double WoodyLiveFuelMoisture(gsi);

  //Combine the live and dead moisture:
  
  fm.M_f_ij[0] = oneHrFM;
  fm.M_f_ij[1] =
  fm.M_f_ij[2] =
  fm.M_f_ij[3] =
  fm.M_f_ij[4] =

  return M_f_ij;
}

/** A stub to use until the minimum daily temperature is available as an input:
 */
double GetDailyMinimumTemperature(double meanDailyAirTempC)
{
	return meanDailyAirTempC - 10;
}

/** Determine fuel loadings based on the cohorts states:
 *
 * This process represts a theory what represent fuels in DVM-DOS-TEM.  Once the mapping is defined
 * the fuel loadings are know, since model states are clearly defined.  Since the mapping itself is
 * a theory it represets a place where assumptions could change.
 *
//The fuels must be estimated from:
// - Live fuels: graminoid, herbaceous, shrub PFTs, and the moss layer.
// - Litter: The first soil layer dead vegetation component (rawc) and woody debris (wdebrisc),
//which currently is only present after a previous fire.
//The mass is upscaled from the carbon stocks.  The litter is distributed into the SAV bins based
//on guesstimates informed the PFT inputs and by some ideas of breakdown rates.
//Previous or current cfall could be used to help estimate this but some assumptions still need
//to be made.
 *
 * In order to be able to update the stocks after a fire occurs it will be a lot simpler if we keep
 * fuels of different sources as separate fuel types, even it they have the same SAV.  This is
 * especially relevant to moss, which could potentially be combined with fine dead material.
 *
 * Either return w_o_ij or update it in the fuel model.  We might as well do the later if we pass
 * the model in.  we could just pass SAV_ij in?  Passing in the fuel model allows it to be updated.
 */
//void CohortStatesToLoading(const CohortData *cd, FuelModel& fm)
void CohortStatesToLoading(const Cohort& cd, FuelModel& fm)
{
  const c2b = 2.0//0.5;//The carbon to biomass ratio for vegetation on a dry basis.
  //Can this be fixed for all vegetation or does it need vary?

  //Dead fuels:
  //Most of the dead surface fuel comes from the rawc component of the first non-moss soil layer.
  //I don't think the soil is is accessible for via the WildFire object or its CohortData.  The
  //layers are linked from the Cohort itself.
  Layer* topFibric = theCohort.soilbgc.ground.fstshlwl;
  //Cast to SoilLayer?  No, rawc is part of Layer.

  double rawC = topFibric->rawc;

  //Distribute the rawc to the dead size classes:
  int numDead = std::count(fm.liveDead.begin(), fm.liveDead.end(), Dead);
  
  std::vector <double> SAV_DeadFM(fm.SAV_ij.begin(), fm.SAV_ij.begin() + numDead);
  
  //For now we assume fm is a standard fuel model and assume the default loadings represent a
  //typical size distribution.  This is probably not the case and a better 
  std::vector <double> w_o_DeadFM(fm.w_o_ij.begin(), fm.w_o_ij.begin() + numDead);
  //Or something like:
  //std::vector <double> sizeWts = GetSizeDistribution();

  //DistributeDeadFuelsD2() is the current R draft of this function.  It needed to be finalized and
  //ported to C++.  THe arguments are SAVsIn, weights, litterMass, SAVsOut.
  std::vector <double> w_o_Dead = DistributeDeadFuelsD2(SAV_DeadFM, w_o_DeadFM, (rawC * c2), SAV_DeadFM);

  //Update the loadings (or do below):
  for (int i = 0; i < numDead; i++)
  {
    fm.w_o_ij[i] = w_o_Dead[i];
  }

  //Live fuels:
  //The live fuels are generally split into herbaceous and woody (shrubs):

  //We have to look through the PFTs...
  //I can't find where the PFT's carbon identifiers and stocks live.  Pseudo-coding for now!:
  for [EACH PFT]
  //for (int j = 0; j < numPFTs; j++)
  {
    if (IS_GRAMINOID(PFT) || IS_HERBACEOUS(PFT))
    {
    	//Put graminoid and herbaceous PFTs in the herbaceous:
    	liveHerbIndex = FuelClassIndex(fm.liveDead, Live, 0);//Can we guarantee the herbaceous will alway be the first live?????
    	fm.w_o_ij[liveHerbIndex] += (StemC + LeafC) * c2;
    }
    else if (IS_WOODYSHRUB(PFT))
    {
    	//Put in the woody.
    }
    //Skip tree PFTs.
  }


  /*Moss is tricky.  It's has it's own biomass stock so it's amount is simple determine but it is a
  live fuel that can act like a dead fuel.  It is also always present in DVM-DOS-TEM but may not
  be what the autors were thinking about when they created the standard fuel models.  It also has
  unique moisture dynamics.  It shouldn't just be shoved into either the fine dead of live
  herbaceous type.
  I think it makes sense to treat it as a live fuel or potentially as a dynamic fuel.  It should be
  given its own SAV if the live herbaceous is too coarse.  The model has a couple of different moss
  PFTs and an SAV could be given for each.*/
  
}

/** Simulate ground fire returning or updating the burn depth...
 *

Burnup produces energy over time so it may be better to link the calculations?

 */
void SimulateGroundFire()
{
  double burnDepth = 0;
  
  //Calculate if the energy output is sufficiency to ignite the surface layer.
  //If not record that ignition failed and return.
  
  //Otherwise continue to calculate progressive smoldering downward.
  //This is the same problem of dying and heating to combustion as we move down.  However, we can
  //safely assume that the fire will not continue if we reach mineral soil, bedrock, permafrost, or
  //the water table.
  
  //return burnDepth;
}
