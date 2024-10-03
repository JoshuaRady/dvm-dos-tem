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

//Constants:
const c2b = 2.0//The carbon to biomass multiplier for vegetation on a dry basis.
//This could vary but many models use a single value.


/** Calculate wildfire behavior and effects using the modeled vegetation, fuels, and meteorology,
 *
The function needs the ecosystem state and meteorology...
 * @param thisCohort The cohort object for this site.
 *   The WildFire object doesn't have all the data we need.  The Cohort gives access to stocks and
 *   meteorology we need.
 * @param monthIndex The current month as a zero based index.
 *

Output:
The output should be stored in the cohort object...
The soil burn depth should be returned/recorded to take advantage of the existing code.

ToDo:
- Add switch for dynamic moisture?

 */
//void RevisedFire(WildFire* wf)//The name will definitely change.
void RevisedFire(const Cohort& thisCohort, int monthIndex)
{

  //Determine the fuel model matching the location's CMT:------------------
  
  //Get the CMT for the grid cell:
  //The CMT number is in the CohortData member of the cohort.  The WildFire object maintains a
  //pointer to it's parent cohort data but it is  private.  This code needs to be in the class,
  //a friend or have access to the cohort.
  //int theCMTnumber = wf->cd.cmttype;//Private data access!!!!!
  int theCMTnumber = thisCohort.cd.cmttype;

  //Crosswalk from CMT to fuel model:
  fuelModelNumber = GetMatchingFuelModel(theCMTnumber);

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
  double slope = wf->slope;//This is percent slope.  Need to convert to fractional slope!!!!!
  //Shortwave radiation may be needed for the moisture calculations.


  //Calculate fuel moisture:----------------------
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
 *
 * @param cmt The number of the CMT for this location.
 *
 * @returns The fuel model number matching the CMT input.	STUB!!!!!
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

/** Determine fuel loadings based on the cohorts states and store then in the sites' fuel model:
 *
 * This process is a theory what represent fuels in DVM-DOS-TEM.  Once the mapping is defined
 * the fuel loadings are know, since model states are clearly defined.  Since the mapping itself is
 * a theory it represets a place where assumptions could change.  The current fuel mapping is:
 *
 * Dead fuels:
 * Dead surface fuels include litter and fallen debris.  In DVM-DOS-TEM these are represented by:
 * - Litter: The rawc component of the first non-moss soil layer.  rawc is part of the 'soil' but
 * is where all aboveground litter-fall ends up.  We assume this is not yet buried.
 * This may include some root cfall as well??????
 * We upscaled the total 'litter' carbon stock to biomass and distribute it among the fuel model's
 * dead fuel size classes.
 * - Large woody debris: Only after fire there is a pool of woody debris (wdebrisc) (and standing
 * dead stems?).
 * We could put all of that in the largest class but some of it might belong in the smaller classes
 * or in a 1000 hour class, which will have a limited effect on fire behavior and is missing from most
 * fuel models.
 * - Possibly live moss (see below).
 *
 * Live fuels:
 * Live fuels include graminoids, forbs, shrubs, and moss (see below).  We convert PFT aboveground
 * carbon stocks to biomass and sort into herbaceous (graminoids & forbs) and woody (shrubs) size
 * classes.
 *
 * Moss:
 * Moss is tricky.  It is a live fuel that can act like a dead fuel, with highly dynamic moisture.
 * The standard fuel models were developed with a focus on the lower 48 so the authors were not
 * thinking a lot about moss.  Moss is hoever an important driver of fire behavior in the north and
 * is ubiquitous in the DVM-DOS-TEM CMTs.  Figuring out how to handle moss is ongoing work.
 *
 * Given this uncertainty we are starting with the two simplest options, place it into either the
 * fine dead or live herbaceous type.  There are several other options.
 * I think it may make sense to treat it as a custom live fuel or potentially as a dynamic fuel.  It
 * should be given its own SAV if the live herbaceous is too coarse.  The model has a couple of
 * different moss PFTs and an SAV could be given for each.
 *
 * We don't include dead moss as a surface fuel.  We treat that as duff, part of the ground fuels.
 *
 * @param thisCohort The cohort object for this site.
 * @param fm The fuel model for the site (with default loadings).
 *
 * returns Nothing but the fuel model loadding (w_o_ij) is updated on return.
 *
 * ToDo:
 * - Record the mapping of stocks to fuels in some way record the value before fire so we can update
 * stocks appropriately afterwards.
 * - Deal with wdebrisc.
 */
void CohortStatesToFuelLoading(const Cohort& theCohort, FuelModel& fm)
{
  //Dead fuels:
  Layer* topFibric = theCohort.soilbgc.ground.fstshlwl;//Get the top non-moss layer.
  double rawC = topFibric->rawc;//Get the total 'litter' carbon mass.

  //Distribute the rawc to the dead size classes:
  
  //Get the dead fuel SAVs:
  int numDead = std::count(fm.liveDead.begin(), fm.liveDead.end(), Dead);//The FuelModel class should provide an interface for this!!!!!
  std::vector <double> savDead(fm.SAV_ij.begin(), fm.SAV_ij.begin() + numDead);

  //Get an estimated distribution of fuel size:
  //For now we assume fm is a standard fuel model and assume the default loadings represent a
  //typical size distribution.  This is probably not the case and a better method for estimating
  //the size distribution will be developed in the future using literature or calculations.
  std::vector <double> w_o_Dead(fm.w_o_ij.begin(), fm.w_o_ij.begin() + numDead);
  //Or something like:
  //std::vector <double> sizeWts = GetSizeDistribution();

  //DistributeDeadFuelsD2() is the current R draft of this function.  It needed to be finalized and
  //ported to C++.  THe arguments are SAVsIn, weights, litterMass, SAVsOut.
  std::vector <double> w_o_Dead = DistributeDeadFuelsD2(savDead, w_o_Dead, (rawC * c2b), savDead);

  //Update the loadings (or do below):
  for (int i = 0; i < numDead; i++)
  {
    fm.w_o_ij[i] = w_o_Dead[i];
  }

  //Live fuels:
  //Standard fuel models have herbaceous and woody classes but in some models both are not actually
  //occupied.  Making sure the fuel classes match the CMT's PFTs is something that should be done
  //in advance but we do checking here for safety.

  //Sort the PFT biomass into the appropriate live fuel classes:
  
  //The fuel model may contain default loadings that should be disregarded:
  liveHerbIndex = FuelClassIndex(fm.liveDead, Live, 0);//Can we guarantee the herbaceous will alway be the first live?????
  liveWoodyIndex = XXXXX;
  fm.w_o_ij[liveHerbIndex] = 0;
  fm.w_o_ij[liveWoodyIndex] = 0;
  
  //I can't find where the PFT's carbon identifiers and stocks live.  Pseudo-coding for now!:
  for (int pftNum = 0; pftNum < NUM_PFT; pftNum++)
  {
    if (!theCohort.cd.d_veg.ifwoody(pftNum))//Or m_veg??????
    {
      //Put all herbaceous PFTs in the herbaceous:
      //Note: There is currently tell graminoids and forbs appart.
      if (nonvascular == 0)
      {
        //Check if class is present!!!!

        //Include aboveground parts:
        double leafC = theCohort.bd[pftNum].m_vegs.c[I_leaf];
        double stemC = theCohort.bd[pftNum].m_vegs.c[I_stem];
        theCohort.bd.fm.w_o_ij[liveHerbIndex] += (leafC + stemC) * c2b;//Convert to dry biomass.
      }
      else//Mosses:
      {
        //Add switch!!!!!
      }
    }
    else//Woody PFTs:
    {
      //Put shrubs in the woody:
      if (IsShrub(theCohort.cd.d_veg, pftNum))
      {
        //Check if class is present!!!!

        //Include aboveground parts:
        double leafC = theCohort.bd[pftNum].m_vegs.c[I_leaf];
        double stemC = theCohort.bd[pftNum].m_vegs.c[I_stem];
        theCohort.bd.fm.w_o_ij[liveHerbIndex] += (leafC + stemC) * c2b;//Convert to dry biomass.
      }
      //Ignore trees.
    }
  }
}

/** Is this PFT a shrub?
 *
 * Currently there is no shrub or tree flag so we have to use other information to infer
 * shrubbiness.  This is a work in progress and may be replaced by a flag later.
 *
 * Since plants can have different growth forms and dwarf statures especially in harsh climates we
 * could try to used stature to help infer shrub status.
 *
 * @param vegState The object containing PFT attributes for this site/
 * @param 
 *
 * @returns True if this PFT is a shrub. 
 */
bool IsShrub(const& vegstate_dim vegState, int pftNum)
{
	
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
 * 			@param fm The fuel model object for this site.			Not currently being used!!!!!
 * @param monthIndex The current month as a zero based index.
 #
 * @returns M_f_ij, the fuel moisture for all fuel classes.  This is not returned in the fuel model
 * passed in because we don't know if curing is being applied.
 
 ToDo:
 - Make Fosberg table files paths available.
 
 */
std::vector <double> CalculateFuelMoisture(const Cohort& thisCohort, int monthIndex)//, const FuelModel& fm)
{
  //std::vector <double> M_f_ij(fm.numClasses, 0);//Return value.
  //We get the number of fuel classes here and then assume the number below.  To handle more than
  //the standard five fuel types we need additional tools to undestand what they are.
  std::vector <double> M_f_ij(fm.numClasses, 0);//Return value.

  //Get the current date:
  //There is a month member in the cohorts CohortData but that apparently is not the current month.
  int monthOfYear = monthOfYearIndex + 1;

  //The day of month is really up to us.  The middle of the month seems resonable.  We'll use the
  //15th for now and add a calculation later?
  int dayOfMonthIndex = 14;//0 based.

  //The daily data provided in Climate is a vector and should be 0 indexed but the code notes that
  //it currently returns 366 days.
  int dayOfYearIndex = temutil::day_of_year(monthIndex, dayOfMonthIndex);


  //Dead fuel moisture:---------------------------

  //Current air temperature:
  //float tempAir = thisCohort.climate.tair[monthIndex]//Monthly
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
  double aspect = thisCohort.cd.cell_aspect;//Degrees

  //Shading should take into consideration both cloudiness and canopy cover:
  bool shaded = false;

  //Canopy cover is equivalent to the projected foliage area.  For the Fosberg NWCG method this is
  //estimated visually.  As such it best to consider this the as the overhead tree canopy.  We can't
  //get this exactly but we can approximate it by excluding the herbaceous and moss leaf area:
  double fpcWoody = 0;
  for (int i = 0; i < NUM_PFT; i++)
  {
    if (thisCohort.cd.d_veg.ifwoody[i])
    {
      fpcWoody += thisCohort.cd.d_veg.fpc[i];
    }
  }

  //To combine the cloudiness (expressed as a percentage) and canopy cover (a fraction) we need to
  //invert the the values so the product increases the effect when combined.  We then re-invert:
  double shadeFrac = 1 - ((1 - (thisCohort.climate.cld_d[dayOfYearIndex] / 100)) * (1 - fpcWoody));

  if (shadeFrac > 0.5)
  {
  	shaded = true;
  }

  //Need to add paths for the Fosberg table files to the config file!!!!!
  double oneHrFM = FosbergNWCG_1HrFM(tableA_Path, tableB_Path, tableC_Path, tableD_Path,//!!!!!
                                     tempAir, rhPct, monthOfYear, hourOfDay,
                                     slopePct, aspect, shaded);//Default values for the rest.

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
  //It would be better to inform the numbers using information from the fuel model.  See above.
  M_f_ij[0] = oneHrFM;
  M_f_ij[1] = tenHrFM;
  M_f_ij[2] = hundredHrFM;
  M_f_ij[3] = herbLFM;
  M_f_ij[4] = woodyLFM;

  return M_f_ij;
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
