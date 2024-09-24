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
  std::vector <double> M_f_ij = CalculateFuelMoisture(theCohort, fm);

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
 * An alternative to calculating this at the time of fire would be to make fuel moisture a
 * continuously calculated state.  This would moke more sense if litter existed as a distinct stock
 * as well.
 *
 * @param thisCohort The cohort object for this site.
 * @param fm The fuel model object for this site.
 *
 * @returns M_f_ij The fuel moisture for all fuel classes.  This is not returned in the fuel model
 * passed in because we don't know if curing is being applied.
 */
//std::vector <double> CalculateFuelMoisture(FuelModel& fm, double tempAir, double slope)
std::vector <double> CalculateFuelMoisture(const Cohort& thisCohort, const FuelModel& fm)
{
  std::vector <double> M_f_ij(fm.numClasses, 0);//Return value.

  //Get the current date:
  //The month is stored in the cohorts CohortData but that may not be the correct place to get it.
  //We may need to get from somewhere else upstream.
  int monthOfYearIndex = thisCohort.cd.month;//0 based
  int monthOfYear = monthOfYearIndex + 1;

  //The day of month is really up to us.  The middle of the month seems resonable.  We'll use the
  //15th for now and add a calculation later?
  int dayOfMonthIndex = 14;//0 based.

  //The daily data provided in Climate is a vector and should be 0 indexed but the code notes that
  //it currently returns 366 days.
  int dayOfYearIndex = temutil::day_of_year(monthOfYearIndex, dayOfMonthIndex);


  //Dead fuel moisture:---------------------------

  //Current air temperature:
  //float tempAir = thisCohort.climate.tair[monthOfYearIndex]//Monthly
  float tempAir = thisCohort.climate.tair_d[dayOfYearIndex];

  //Calculate current relative humidity:

  //Atmospheric pressure is needed for the the Fireweed code method.
  //Sea level barometric pressure will be available soon but is not currently.  Use this along with
  //thisCohort.cd.cell_elevation() to get the pressure at this elevation.
  //p_hPa = Z;

  //Partial pressure of water vapor:
  //This is an input variable.  There has been some discussion recently about the units of this.
  //The docs say it is in hPa but it is used in the code as if it is in in Pa
  //(see Climate.cpp calculate_vpd()).
  float P_Pa = thisCohort.climate.vapo_d[dayOfYearIndex];//My reading is that vapo_d is a vector of daily values for the whole year.

  //Saturated vapor pressure:
  //double P_s_hPa = SaturationVaporPressureBuck(tempAir, p_hPa);//Fireweed method returns hPa.
  //This is a computed climate variable.
  float P_s_Pa = thisCohort.climate.svp_d[dayOfYearIndex];//From the code this must be in Pa.

  //Use the pressures to calculate relative humidity:
  double rhPct = RHfromVP(P_Pa, P_s_Pa);

  //We assume that fires occur in the mid-afternoon.  The exact hour might need to be shared across
  //the entire fire code.  The Fosberg model was originally devise with afternoon (~14:30) values in
  //mind, but the table method used below can take any daylight time.
  int hourOfDay = 15;

  //Slope:
  //Also duplicated in the Wildfire object.
  double slopePct = thisCohort.cd.cell_slope;//Percent

  //Get the aspect:
  //Also duplicated in the Wildfire object.
  double aspect = thisCohort.cd.cell_aspect;//Units are probably degrees?????

  //Shading:
  //Shading should take into consideration both cloudiness and canopy cover but we arr currently
  //only able to address the later.
  bool shaded = false;
  if (thisCohort.climate.cld_d[dayOfYearIndex] > 50.0)
  {
  	shaded = true;
  }

  //Need to add paths to config file!!!!!
  double oneHrFM = FosbergNWCG_1HrFM(std::string tableA_Path, std::string tableB_Path, std::string tableC_Path,
                                     std::string tableD_Path,
                                     tempAir, rhPct, monthOfYear, hourOfDay,
                                     slopePct,
                                     aspectCardinal,//Need to convert from degrees.
                                     shaded,
                                     'L',//The default value for elevation could be omitted if the default for units was changed to metric.
                                     Metric);

  double tenHrFM = NWCG_10hrFM(double oneHrFM);
  double hundredHrFM = NWCG_100hrFM(double oneHrFM);

  //Live fuel moisture:---------------------------
  //The GSI based fuel moistures should be averages of the 21 days up to today:
  //Monthly values /might work but I suspect minimum daily temp could cause big departures at some
  //times of the year.

  double herbLFM = 0;
  double woodyLFM = 0;
  for (int i = dayOfYearIndex - 20; i <= dayOfYearIndex; i++)
  {
    //Minimum and maximum daily temperature are not currently available but will be soon.  We
    //use a hack value for now:
    float tempCMin = thisCohort.climate.tair_d[i] - 10;
    
    //Get the VPD:
    float vpdPa = thisCohort.climate.vpd_d[i];
    //Or calculate VPD from mean daily weather values:
    //double p_hPa = ?????;
    //...
    //double vpd_hPa = VPDfromRHBuck(tempAir, rh, p_hPa);

    //temutil::length_of_day gives the he day length in hours:
    float dayLengthSec = temutil::length_of_day(theCohort->lat, dayOfYearIndex) * 60 * 60;

    double gsi = GrowingSeasonIndex(tempCMin, vpdPa, dayLengthSec);
    herbLFM += HerbaceousLiveFuelMoisture(gsi);
    woodyLFM += WoodyLiveFuelMoisture(gsi);
  }
  //Average:
  herbLFM = herbLFM / 21;
  woodyLFM = woodyLFM / 21;


  //Combine the live and dead moisture:-----------
  //This makes assumption that the order is that of a standard fuel model.
  fm.M_f_ij[0] = oneHrFM;
  fm.M_f_ij[1] = tenHrFM;
  fm.M_f_ij[2] = hundredHrFM;
  fm.M_f_ij[3] = herbLFM;
  fm.M_f_ij[4] = woodyLFM;

  return M_f_ij;
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
