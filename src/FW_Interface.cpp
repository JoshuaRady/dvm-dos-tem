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
 * 
 * All of the functions requiring direct access to WildFire members have been moved into the class
 * but have left here pending testing so the new code can be more clearly identified.  There are
 * number of attendant functions. most of which can also be moved into the class but are being left
 * as is for now to make dependancies stark, i.e. they do not need access to member data or
 * functions.
 */

#include "../include/FW_Interface.h"
#include "../include/TEMLogger.h"
#include "../include/TEMUtilityFunctions.h"//For length_of_day().
#include "../include/WildFire.h"

#include "FireweedDeadFuelMoistureFosberg.h"
#include "FireweedFuelTools.h"
#include "FireweedLiveFuelMoistureGSI.h"
#include "FireweedMetUtils.h"

#include <cmath>//Temporary for isnan().

extern src::severity_logger< severity_level > glg;

//Constants:
const double c2b = 2.0;//The carbon to biomass multiplier for vegetation on a dry basis.
//This could vary by PFT but many models use a single value.


/** Calculate wildfire behavior and effects using the modeled vegetation, fuels, and meteorology.
 *
 * This is the main entry point for the revised wildfire model.
 * The function needs the ecosystem state, meteorology, and time of year as inputs.
 *
 * JMR_NOTE: I attempted multiple ways to pass this information in but ultimately found it necessary to
 * modify WildFire to provide access to most of the necessary inputs.
 * 
 * @param[in] monthIndex The current month as a zero based index.
 *
 * @returns The soil burn depth from ground fire (meters).  Other output are stored in model objects.
 * 
 * @note As discussed in a note in WildFire::burn() the burn depth is stored in the FirData member so the
 * return value is not strictly needed.  We are keeping this code parallel to getBurnOrgSoilthick()
 * for now.
 */
double WildFire::RevisedFire(const int monthIndex)//Name could change.
{
  BOOST_LOG_SEV(glg, debug) << "Entering WildFire::RevisedFire()...";

  //Gather weather and environmental conditions:---------------------------
  double tempAir = edall->d_atms.ta;//Daily air temp (at surface).
  //Can we use temperature from climate instead?????
  //Current humidity is needed for calculating fuel moisture but that code handles it itself.
  //We may also need it for duff moisture soon.
  
  double windSpeed = GetMidflameWindSpeed();
  
  //The percent slope is stored in the CohortData object and also in the wildfire object:
  double slopeSteepness = SlopePctToSteepness(cd->cell_slope);

  //Shortwave radiation may be needed for more advanced moisture calculations added in the future.


  //Determine the fuel model for the location and set its parameters:------

  //Get the CMT for the grid cell:
  //JMR_NOTE: The CMT number is in the CohortData member of the cohort.  The WildFire object maintains a
  //pointer to it's parent cohort data but it is private.  This code needs to be in the class,
  //a friend or have access to the cohort.
  int theCMTnumber = cd->cmttype;

  //Crosswalk from CMT to the matching fuel model:
  int fuelModelNumber = GetMatchingFuelModel(theCMTnumber);

  //Load the fuel model from the fuel model table file:
  FuelModel fm = GetFuelModelFromCSV(md.fire_fuel_model_file, fuelModelNumber);//Or tab delimited!!!!!
  fm.ConvertUnits(Metric);//Convert to metric units.

  //Determine the surface fuels from the model vegetation and soil states and update the fuel
  //loadings from their default values:
  CohortStatesToFuelLoading(fm, md.fire_moss_as_dead_fuel);

  //Save the fuel loading prior to fire: (will be compared below...)
  std::vector <double> fuelLoadingBefore = fm.w_o_ij;

  //Calculate the fuel bed depth:
  CalculateFuelBedDepth(fm, md.fire_calculate_delta);

  //Calculate fuel moisture:
  //Note: It is better to calculate fuel moisture after calculating fuel loadings since that process
  //might change the fuel sizes.
  std::vector <double> M_f_ij = CalculateFuelMoisture(fm, monthIndex);

  //Add the moisture to the fuel model possibly computing dynamic fuel moisture:
  BOOST_LOG_SEV(glg, debug) << "Apply fuel moisture to fuel model...";
  if (md.fire_dynamic_fuel)
  {
    fm.CalculateDynamicFuelCuring(M_f_ij);
  }
  else
  {
    fm.SetFuelMoisture(M_f_ij);
  }

  //Feed fuels and weather conditions into the surface fire models:--------

  //Dump the fuel model if debugging:  This may be temporary?????
  BOOST_LOG_SEV(glg, debug) << "Dumping the fuel model prior to fire:";
  BOOST_LOG_SEV(glg, debug) << fm;

  //First to Rothermel & Albini spread rate model:
  //It takes a fuel model (and attendant data) and returns the calculation details.
  //M_f_ij does not need to be included since it is added to the fuel model object above.
  BOOST_LOG_SEV(glg, debug) << "Perform surface file spread rate calculations...";
  SpreadCalcs raData = SpreadCalcsRothermelAlbini_Het(fm, windSpeed, slopeSteepness);
                                                      //M_f_ij,//fuel moisture
                                                      //false,//useWindLimit
                                                      //FALSE);//debug

  //Dump the output of the spread rate calculations if debugging:
  //JMR_NOTE: This may be overkill?????
  BOOST_LOG_SEV(glg, debug) << "Dump the spread calculations:";
  BOOST_LOG_SEV(glg, debug) << raData;
//   {
//   	BOOST_LOG_SEV(glg, debug) << "WildFire::RevisedFire() Spread rate calculation R = " << raData.R;
//   }


  //Simulate the combustion of surface fuels:------------------------------
  BurnupSim = SimulateSurfaceCombustion(fm, raData, tempAir, windSpeed);
  BOOST_LOG_SEV(glg, debug) << "Dump the combustion calculations:";
  BOOST_LOG_SEV(glg, debug) << BurnupSim;

  //Update litter, moss, and soil carbon stocks:
  //Some of this could be done directly here but we could also use parts of the 'original' code like
  // WildFire::updateBurntOrgSoil(), which I split out.


  //Add crown fire!!!!!

  //Simulate ground fire:--------------------------------------------------
  //Use the energy flux from the aboveground fire into the soil surface (RA + Burnup + crown) and
  //the fire air temp (Burnup fire environmental temperature?) as input to the ground fire model:
  //...
  //Also pass in the soil profile conditions...

  double burnDepth = SimulateGroundFire();

  //Match messaging of getBurnOrgSoilthick():
  BOOST_LOG_SEV(glg, debug) << "Setting the burn thickness in FirData...";
  fd->fire_soid.burnthick = burnDepth;

  BOOST_LOG_SEV(glg, info) << "Final Calculated Organic Burn Thickness: " << burnDepth;
  return burnDepth;
}

/** Determine the fuel model (number) matching a given community type:
 *
 * Each CMT has [will have] a predetermined fuel model assigned to it via a parameter file.
 * It needs to be determined if this will be the CMT parmameter file or an additional lookup table.
 *
 * @param[in] cmt The number of the CMT for this location.
 *
 * @returns The fuel model number (not the index) matching the CMT input.	STUB!!!!!
 */
int GetMatchingFuelModel(const int cmt)
{
  //Get the number of the fuel model from the crosswalk in the parameter files.
  //This crosswalk needs to be made!!!!!
  int fuelModelNumber = 161;//Temporarily hardwired.
  //We use TU1 = 161 since it is good for testing.  All fuel types are occupied and it is dynamic.
  
  //There is a bare land fuel model by number but it doesn't have parameters.  Do we need to provide
  //a bare land parameter set or can we just signal the calling code that it should skip fire
  //calculations?
  
  //If no match either throw an error or warn and return a default fuel model.

  return fuelModelNumber;
}

/** Determine fuel loadings based on the cohorts states and store then in the sites' fuel model:
 *
 * This process is a theory of what represent fuels in DVM-DOS-TEM.  Once the mapping is defined
 * the fuel loadings are know, since model states are clearly defined.  Since the mapping itself is
 * a theory it represets a place where assumptions could change.  The current fuel mapping is:
 *
 * @par Dead fuels:
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
 * @par Live fuels:
 * Live fuels include graminoids, forbs, shrubs, and moss (see below).  We convert PFT aboveground
 * carbon stocks to biomass and sort into herbaceous (graminoids & forbs) and woody (shrubs) size
 * classes.
 *
 * @par Moss:
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
 * @param[in,out] fm The fuel model for the site (with default loadings).
 * @param[in] treatMossAsDead Should moss be treated as a fine dead fuel (true) or herbaceous live fuel (false)?
 *
 * @returns Nothing but the fuel model loadding (w_o_ij) is updated on return.
 *
 * ToDo:
 * - Record the mapping of stocks to fuels in some way record the value before fire so we can update
 * stocks appropriately afterwards.  We have done this to some extend in the calling code.
 * - Deal with wdebrisc.
 */
void WildFire::CohortStatesToFuelLoading(FuelModel& fm, const bool treatMossAsDead)
{
  BOOST_LOG_SEV(glg, debug) << "Entering WildFire::CohortStatesToFuelLoading()...";

  const double gPerKg = 1000;//Move to FireweedUnits.h?
  
  //Dead fuels:
  double rawC = GetLitterRawC();//Get the total litter carbon.

  //Distribute the rawc among the dead fuel size classes:

  //Get the dead fuel SAVs:
  int numDead = fm.NumDeadClasses();
  std::vector <double> savsDead(fm.SAV_ij.begin(), fm.SAV_ij.begin() + numDead);

  //Get an estimated distribution of dead fuel sizes:
  std::vector <double> distribSAVs, distribWts;//Empty vectors to hold the distribution.
  GetDeadFuelSizeDistribution(fm, distribSAVs, distribWts);

  //Distribute the litter carbon using the distribution:
  //All carbon stocks are in g/m^2 C and need to be converted to kg/m^2 dry biomass.
  std::vector <double> w_o_Dead = DistributeFuel(distribSAVs, distribWts, (rawC * c2b / gPerKg), savsDead);

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
  int liveHerbIndex = fm.LiveHerbaceousIndex();
  int liveWoodyIndex = fm.LiveWoodyIndex();
  fm.w_o_ij[liveHerbIndex] = 0;
  fm.w_o_ij[liveWoodyIndex] = 0;
  
  for (int pftNum = 0; pftNum < NUM_PFT; pftNum++)
  {
    if (!cd->d_veg.ifwoody[pftNum])//Or m_veg??????
    {
      //Put all herbaceous PFTs in the herbaceous:
      //Note: There is currently no way to tell graminoids and forbs appart.
      if (cd->d_veg.nonvascular[pftNum] == 0)
      {
        //If the herbaceous class is not present adding carbon to it will not influence the fire
        //behavior.  If herbaceous fuels are not important in this system we could ignore them.
        //This is not something we can really resolve at run time.  This is a science question that
        //needs to be addressed during fuel model selection.
        if (!fm.LiveHerbaceousPresent())
        {
          BOOST_LOG_SEV(glg, fatal) << "The live herbaceous fuel type is not active in this fuel model.";
        }

        //Include aboveground parts:
        double leafC = bd[pftNum]->m_vegs.c[I_leaf];
        double stemC = bd[pftNum]->m_vegs.c[I_stem];
        fm.w_o_ij[liveHerbIndex] += (leafC + stemC) * c2b / gPerKg;//Convert to dry biomass.
      }
      else//Mosses:
      {
        //My best reading is that moss is all leaf in DVM-DOS-TEM.  However if they had stems and
        //'roots', i.e. rhizoids = root C, they should burn too.  Include just in case for now:
        double leafC = bd[pftNum]->m_vegs.c[I_leaf];
        double stemC = bd[pftNum]->m_vegs.c[I_stem];//Should be 0.
        double rootC = bd[pftNum]->m_vegs.c[I_root];//Should be 0.
        double mossBiomass = (leafC + stemC + rootC) * c2b / gPerKg;//Convert to dry biomass.

        if (treatMossAsDead)
        {
          fm.w_o_ij[1] += mossBiomass;//Assumes fine fuel is first, which is pretty safe.
        }
        else
        {
          //See notes above.
          if (!fm.LiveHerbaceousPresent())
          {
            BOOST_LOG_SEV(glg, fatal) << "The live herbaceous fuel type is not active in this fuel model.";
          }

          fm.w_o_ij[liveHerbIndex] += mossBiomass;
        }
      }
    }
    else//Woody PFTs:
    {
      //Put shrubs in the live woody fuel:
      if (IsShrub(cd->cmttype, pftNum))
      {
        //If the woody class is not present adding carbon to it will not influence the fire
        //behavior.  See notes for herbaceous fules above.
        if (!fm.LiveHerbaceousPresent())
        {
          BOOST_LOG_SEV(glg, fatal) << "The live woody fuel type is not active in this fuel model.";
        }

        //Include aboveground parts:
        double leafC = bd[pftNum]->m_vegs.c[I_leaf];
        double stemC = bd[pftNum]->m_vegs.c[I_stem];
        fm.w_o_ij[liveWoodyIndex] += (leafC + stemC) * c2b / gPerKg;//Convert to dry biomass.
      }
      //Ignore trees.
    }
  }
}

/** FW_MOD: Get the litter carbon for the site.
 *
 * The 'litter' is defined as the rawC componenet of the first non-moss soil layer.  See
 * CohortStatesToFuelLoading() for more information.
 * 
 * @note This started as a hack and could be moved into CohortStatesToFuelLoading().  However, it
 * does make the code a bit more organized.
 */
double WildFire::GetLitterRawC() const
{
  BOOST_LOG_SEV(glg, debug) << "Entering WildFire::GetLitterRawC()...";

  int topFibricIndex;

  BOOST_LOG_SEV(glg, debug) << "Getting litter carbon.";

  //Determine the top non-moss layer (fstshlwl in Ground terminology):
  //The site's Ground object is a private member of the parent cohort's Soil_Bgc object.  If we
  //could access this we wouldn't have to search the soil layers:
  //Layer* topFibric = thisCohort.soilbgc.ground->fstshlwl;//Get the top non-moss layer.
  //double rawC = topFibric->rawc;//Get the total 'litter' carbon mass.

  for (int i = 0; i < cd->m_soil.numsl; i++)//Assumes we are starting at the top layer.
  {
    if (cd->m_soil.type[i] == 1)//Shallow organic / peat ~ I_FIB
    {
      topFibricIndex = i;
      break;
    }
  }

  return bdall->m_sois.rawc[topFibricIndex];
}

/** Get the estimated dead fuel size distribution for the site.
 * The distribution is returned via the vectors passed.
 *
 * This is a temporary stub-ish implementation of the code.  For now we assume fm is a standard fuel
 * model and that the default loadings represent a typical dead fuel size distribution.  This is
 * probably not the case and a better method for estimating the size distribution will be developed
 * in the future using literature or calculations.  Then this function will be updated to some
 * lookup process using the CMT of fuel model number to get the estimate that is calculated or
 * supplied via an input file.
 * 
 * @param[in] fm A fuel model to base the distribution on. [Temporary!!!!!]
 * @param[out] distribSAVs An empty vector to return the SAVs of the distribution in.
 * @param[out] distribWts An empty vector to return the weights of the distribution in.
 * 
 * @returns Nothing.  Parameters are updated instead.
 */ 
void GetDeadFuelSizeDistribution(const FuelModel& fm, std::vector <double>& distribSAVs,
                                 std::vector <double>& distribWts)
{
  BOOST_LOG_SEV(glg, debug) << "Entering GetDeadFuelSizeDistribution()...";

  //Copy the appropriate dead members:
  int numDead = fm.NumDeadClasses();
  distribSAVs.assign(fm.SAV_ij.begin(), fm.SAV_ij.begin() + numDead);
  distribWts.assign(fm.w_o_ij.begin(), fm.w_o_ij.begin() + numDead);
}

/** Is this PFT a shrub?
 *
 * Currently there is no shrub or tree flag so we have to use other information to infer
 * shrubbiness.  This is a work in progress and may be replaced by a flag later.
 *
 * The vegetation parameter files include names for most of the PFTS, but only in comments, which
 * are not imported.  These include DecShrub and EvrShrub.  The order of PFTs in a CMT are also not
 * systematic.  The only approach that remains is to check the PFT numbers agains their CMTs.  This
 * requires knowledge of the CMT file (cmt_bgcvegetation.txt).  Since the CMTs will change in due
 * course this implementation is not robust and should be replaced soon.
 *
 * Since plants can have different growth forms and dwarf statures, especially in harsh climates, we
 * could try to used stature to help infer when trees are shrubby in a way that is relevant to fuel
 * models.  That might be getting a bit fancy.
 *
 * @param[in] cmtNumber The CMT number.
 * @param[in] pftIdx The index of the PFT to check.
 *
 * @returns True if this PFT is a shrub. 
 */
bool IsShrub(const int cmtNumber, const int pftIdx)
{
  BOOST_LOG_SEV(glg, debug) << "Entering IsShrub()...";

  //I believe all the following PFTs are considered shrubs:
  //A few of these CMT cases could be combined but I'm prioritizing readability over compactness.
  switch (cmtNumber) {
    case 0://Bare ground, contains junk values.
      break;

    case 1://Boreal Black Spruce
    case 2://Boreal White Spruce Forest
    case 3://Boreal Deciduous Forest
      if (pftIdx == 1)//DecidShrub
      {
        return true;
      }
      break;

    case 4://Shrub Tundra (Toolik area)
      if (pftIdx >= 0 && pftIdx <= 3)//Salix, Betula, Decid, EGreen
      {
        return true;
      }
      break;

    case 5://Tussock Tundra (2022)
      if (pftIdx >= 0 && pftIdx <= 2)//Betula, Decid, EGreen
      {
        return true;
      }
      break;

    case 6://Wet Sedge Tundra (Toolik area)
      if (pftIdx == 0)//Decid
      {
        return true;
      }
      break;

    case 7://Heath Tundra (2019)
      if (pftIdx >= 0 && pftIdx <= 1)//Decid, EGreen
      {
        return true;
      }
      break;

    case 12://Lowland Boreal Wet Shrubland
      if (pftIdx >= 0 && pftIdx <= 1)//DecShrub, EvrShrub
      {
        return true;
      }
      break;

    case 20://EML Shrub Tundra
      if (pftIdx == 0)//Betnan (Betuala nana)
      {
        return true;
      }
      break;

    case 21://EML Tussock Tundra (2022)
      if (pftIdx >= 0 && pftIdx <= 1)//Decidsh, Egreensh 
      {
        return true;
      }
      break;

    case 31://Boreal Bog
      if (pftIdx >= 2 && pftIdx <= 3)//DecShrub, EvrShrub
      {
        return true;
      }
      break;

    case 44://Shrub Tundra (Kougarok)
      if (pftIdx >= 0 && pftIdx <= 3)//Salix, Betula, Decid, EGreen
      {
        return true;
      }
      break;

    default:
      BOOST_LOG_SEV(glg, fatal) << "IsShrub() does not know this CMT.";
      break;
  }
  return false;
}

/** Calculate the fuel bed depth and update it in the fuel model passed.
 *
 * Standard fuel models include a fuel bed depth among their default values.  We need a way to
 * estimate this from the model along with fuel loadings.  DVM-DOS-TEM doesn't give us much to work
 * with.  Ecosystem structural information would be most useful but we don't have that.  For now we
 * will make some assumptions.  Here are some options:
 * 
 * - The fuel bed depth is constant.  This may be reasonable for gramaniod dominated ecosystems,
 * e.g. grasslands, tussock tundra, etc., where the vegetation has a typical maximum height but may
 * be more or less dense.
 * 
 * - The mix of fuels and it structure is fairly constant but the amount is not.  That implies the
 * bulk density of the fuel bed is constant but the height varies with fuel loading.  Fuel bed depth
 * is then calculable from fuel loading and an estimate of they typical bulk density.
 *
 * The reality is likely more complicated but we will build our approach around these for now.
 *
 * @param[in,out] fm The fuel model for the site (with default fuel bed depth and updated fuel loadings).
 * @param[in] dynamic Temporary: Should the fuel bed depth be treated as changing with fuel amounts?
 *
 * @returns Nothing.  The fuel bed depth is updated in the fuel model passed.
 */
void CalculateFuelBedDepth(FuelModel& fm, const bool dynamic)
{
  BOOST_LOG_SEV(glg, debug) << "Entering CalculateFuelBedDepth()...";

  //Determine which approach to use based on the fuel model or CMT would be ideal but it will
  //require some research.  For now we use a switch.  If the depth is constant there is nothing to
  //do.  Otherwise calculate the depth based on constant bulk density:
  if (dynamic)
  {
    double totalLoading = std::accumulate(fm.w_o_ij.begin(), fm.w_o_ij.end(), 0.0);
    //This assumes fm.bulkDensity is the standard fuel model default value:
    fm.delta = totalLoading / fm.bulkDensity;//kg/m^2 / kg/m^3 = m
  }
}

/** Get the wind speed at midflame height in m/min.
 *
 * The Rothermel Albini spread model takes midflame wind speed, an ill defined quantity representing
 * the mean wind speed in the flaming zone.  An estimate of near surface wind speed is what we need.
 * 2 meter wind speed should be fine.  The input may not match this so we need to compute our best
 * estimate.
 *
 * We need a sub-daily value ~ daily value?
 * The code has to do the following:
 * - Get the daily windspeed from the host model.
 *   Note: Windspeed is not yet available in DVM-DOS-TEM, but it should be soon.
 * - Calculate the windspeed at ~2m if the provided wind speed is at another height.
 * - Possibly: Estimate a sub-daily from a daily mean value.  The fire may occur on timescale where
 *   it occurs at a specific time of day.  In this case it would be good to estimate what the wind
 *   speed was at this time.  This may be difficult as wind can be highly variable and may not be
 *   well predicted by diurnal cycles.
 * - Convert to m/min if needed.
 *
 * @returns The wind speed at 2 meters height (m/min).
 *
 * @note: We could also pass in the desired height.  We currently need ~2m and will have 2m wind
 * input so height adjustment is not necessary.  Passing in the time of day would allow us to adjust
 * daily to sub-daily wind (see above).
 */
double WildFire::GetMidflameWindSpeed() const
{
  double windSpeed;//Return value.
  
  //Draft:
  //If this is daily how do we know what day of the month it is?
  //The wind speed will be provided as directional components in m/s at 2 meters.
  //double u = edall.d_atms.EasternWindSpeed;//Zonal component U.		Confirm units!!!!!
  //double v = edall.d_atms.NorthernWindSpeed;//Meridional component V.
  //windSpeed = std::sqrt(std::pow(u, 2) + std:pow(v, 2));//Get the vector wind speed.
  //windSpeed *= 60;//Convert from m/s to m/min.
  //This a daily average value.  An afternoon value would probably be better.

  //Temporary stub, return an arbitrary value!!!!!:
  //Summer average wind speeds are ~6 mph in Fairbanks Alaska.
  //6 * ftPerMi / 60 / ftPerM = 160.9344
  windSpeed = 160.9344;
  
  return windSpeed;
}

/** Calculate fuel moisture based on recent weather:
 *
 * An alternative to calculating this at the time of fire would be to make fuel moisture a
 * continuously calculated state.  This would moke more sense if litter existed as a distinct stock
 * as well.
 *
 * @param[in] fm The fuel model object for this site.
 * @param[in] monthIndex The current month as a zero based index.
 #
 * @returns M_f_ij, the fuel moisture for all fuel classes.  This is not returned in the fuel model
 * passed in because we don't know if curing is being applied.
 */
std::vector <double> WildFire::CalculateFuelMoisture(const FuelModel& fm, const int monthIndex) const
{
  BOOST_LOG_SEV(glg, debug) << "Entering WildFire::CalculateFuelMoisture()...";

  //We get the number of fuel classes here and then assume the number below.  To handle more than
  //the standard five fuel types we need additional tools to undestand what they are.
  std::vector <double> M_f_ij(fm.numClasses, 0);//Return value.

  //Get the current date:
  //There is a month member in the cohorts CohortData but that apparently is not the current month.
  int monthOfYear = monthIndex + 1;

  //The day of month is really up to us.  The middle of the month seems resonable.  We'll use the
  //15th for now and add a calculation later?
  int dayOfMonthIndex = 14;//0 based.

  //The daily data provided in Climate is a vector and should be 0 indexed but the code notes that
  //it currently returns 366 days.
  int dayOfYearIndex = temutil::day_of_year(monthIndex, dayOfMonthIndex);


  //Dead fuel moisture:---------------------------
  BOOST_LOG_SEV(glg, debug) << "Calculating dead fuel moisture:";//Temporary?????

  //Current air temperature:
  //float tempAir = climate->tair[monthIndex]//Monthly
  float tempAir = climate->tair_d[dayOfYearIndex];

  //Calculate current relative humidity:

  //Atmospheric pressure is needed for the the Fireweed code method.
  //Sea level barometric pressure will be available soon but is not currently.  Use this along with
  //cd->cell_elevation() to get the pressure at this elevation.
  //p_hPa = Z;

  //Partial pressure of water vapor:
  //This is an input variable.  There has been some discussion recently about the units of this.
  //The docs say it is in hPa but it is used in the code as if it is in Pa
  //(see Climate.cpp calculate_vpd()).
  float P_Pa = climate->vapo_d[dayOfYearIndex];//My reading is that vapo_d is a vector of daily values for the whole year.

  //Saturated vapor pressure:
  //double P_s_hPa = SaturationVaporPressureBuck(tempAir, p_hPa);//Fireweed method returns hPa.
  //This is a computed climate variable.
  float P_s_Pa = climate->svp_d[dayOfYearIndex];//From the code this must be in Pa.

  //Use the pressures to calculate relative humidity:
  double rhPct = RHfromVP(P_Pa, P_s_Pa);

  //We assume that fires occur in the mid-afternoon.  The exact hour might need to be shared across
  //the entire fire code.  The Fosberg model was originally devise with afternoon (~14:30) values in
  //mind, but the table method used below can take any daylight time.
  int hourOfDay = 15;

  //The percent slope is stored in he CohortData object and also in the wildfire object:
  double slopePct = cd->cell_slope;//Percent

  //Get the aspect:
  //Slope is also duplicated in the Wildfire object.
  double aspect = cd->cell_aspect;//Degrees

  //Shading should take into consideration both cloudiness and canopy cover:
  bool shaded = false;

  //Canopy cover is equivalent to the projected foliage area.  For the Fosberg NWCG method this is
  //estimated visually.  As such it best to consider this the as the overhead tree canopy.  We can't
  //get this exactly but we can approximate it by excluding the herbaceous and moss leaf area:
  double fpcWoody = 0;
  for (int i = 0; i < NUM_PFT; i++)
  {
    if (cd->d_veg.ifwoody[i])
    {
      fpcWoody += cd->d_veg.fpc[i];
    }
  }

  //To combine the cloudiness (expressed as a percentage) and canopy cover (a fraction) we need to
  //invert the the values so the product increases the effect when combined.  We then re-invert:
  double shadeFrac = 1 - ((1 - (climate->cld_d[dayOfYearIndex] / 100)) * (1 - fpcWoody));

  if (shadeFrac > 0.5)
  {
  	shaded = true;
  }

  //Perform the 1hr moisture look up and derive the rest from that:
  BOOST_LOG_SEV(glg, debug) << "Calling FosbergNWCG_1HrFM()";//Temporary?????
  double oneHrFM = FosbergNWCG_1HrFM(md.fire_fosberg_a_file, md.fire_fosberg_b_file,
                                     md.fire_fosberg_c_file, md.fire_fosberg_d_file,
                                     tempAir, rhPct, monthOfYear, hourOfDay,
                                     slopePct, aspect, shaded);//Default values for the rest.

  double tenHrFM = NWCG_10hrFM(oneHrFM);
  double hundredHrFM = NWCG_100hrFM(oneHrFM);

  //Live fuel moisture:---------------------------
  //The GSI based fuel moistures should be averages of the 21 days up to today:
  //Monthly values might work but I suspect minimum daily temp could cause big departures at some
  //times of the year.
  BOOST_LOG_SEV(glg, debug) << "Calculating live fuel moisture:";//Temporary?????

  double herbLFM = 0;
  double woodyLFM = 0;
  for (int i = dayOfYearIndex - 20; i <= dayOfYearIndex; i++)
  {
    //Minimum and maximum daily temperature are not currently available but will be soon.  We
    //use a hack value for now:
    float tempCMin = climate->tair_d[i] - 10;
    
    //Get the VPD:
    float vpdPa = climate->vpd_d[i];
    //Or calculate VPD from mean daily weather values:
    //double p_hPa = ?????;
    //...
    //double vpd_hPa = VPDfromRHBuck(tempAir, rh, p_hPa);

    //temutil::length_of_day gives the day length in hours:
    float dayLengthSec = temutil::length_of_day(lat, dayOfYearIndex) * 60 * 60;

    double gsi = GrowingSeasonIndex(tempCMin, vpdPa, dayLengthSec);
    herbLFM += HerbaceousLiveFuelMoisture(gsi);
    woodyLFM += WoodyLiveFuelMoisture(gsi);
  }
  //Average:
  herbLFM = herbLFM / 21;
  woodyLFM = woodyLFM / 21;


  //Combine the live and dead moisture:-----------
  BOOST_LOG_SEV(glg, debug) << "Setting M_f_ij:";//Temporary?????
  //This makes assumption that the order is that of a standard fuel model.
  //It would be better to inform the numbers using information from the fuel model.  See above.
  //The units need to change from percent moisture to moisture fraction.
  M_f_ij[0] = oneHrFM / 100;
  M_f_ij[1] = tenHrFM / 100;
  M_f_ij[2] = hundredHrFM / 100;
  M_f_ij[3] = herbLFM / 100;
  M_f_ij[4] = woodyLFM / 100;

  return M_f_ij;
}

/** Simulate combustion of surface fuels.
 *
 * We simulate surface fuel combustion with the Burnup model of Albini 1995.
 *
 * @param[in] fm The fuel model for the site.
 * @param[in] raData The results of a Rothermel Albini spread model calculation for the site.
 * @param[in] tempAir The ambient air temperature (C).
 * @param[in] windSpeed The wind speed at 2 meters height (m/min).
 *
 * @returns An object containing the simulation results.
 */
BurnupSim SimulateSurfaceCombustion(const FuelModel& fm, const SpreadCalcs raData,
                                    const double tempAir, const double windSpeed)//CalculateSurfaceCombustion?
{
  BOOST_LOG_SEV(glg, debug) << "Entering SimulateSurfaceCombustion()... [In progress]";

  //Burnup takes a number of parameters:
  //Wind and air temperature are needed (passed in):

  //Fuel properties are a superset of those in a standard fire behavior fuel model:
  //Leave all remaining properties at their default values set in BurnupFM().

  //Duff loading and moisture:
  /*It is unclear exactly what should be used here.  The duff concept in Burnup is certainly not
  synonymous with the deep duff of the boreal.  The shallow duff / live moss layer is probably the
  appropriate analog.  We handle that as a surface fuel so we do not include a duff loading to avoid
  double counting.  We treat smoldering ground fire as a separate stage too.  Duff in Burnup treated
  very simply as an additional source of heat rather than a full blown fuel.  Omitting duff seems
  reasonable but we should revisit this assumption in the future.*/
  double duffLoading = 0;
  double duffMoisture = 0;

  //Fire properties are obtained from the spread rate calculations:
  //The model takes fire intensity and residence time.  These can be estimated form the spread rate
  //calculations.

  //Fire intensity of the flaming front must be estimated.  The appropriate calculation to use is
  //the question.  The reaction intensity represents the energy flux generated by the fire.  The
  //extent to which this may need to be reduced to reflect energy lost needs to considered.
  //double fireIntensity = fm.I_R / 60;//Convert kJ/m^2/min -> kW/m^2 kJ/m^2/min
  //kJ/m^2/min = kJ/min/m^2 = kJ/(60 * s) /m^2 =  kJ/s / 60 /m^2 = kW / 60 /m^2 = (kW/m^2) / 60
  //The heat source takes the propagating flux ratio, wind, and slope into consideration:
  double fireIntensity = fm.heatSource / 60;//Convert kJ/m^2/min -> kW/m^2 kJ/m^2/min

  double t_r = ResidenceTime(fm.cSAV, Metric) * 60;//Convert minutes -> seconds.

  //Simulation settings:
  //Start with example settings from FOFEM examples.  We can experiment with these.  Leave other at
  //their default values.
  double dT = 15;
  int nTimeSteps = 3000;

  //Call Burnup:
  //The simulation calculates the both consumption of fuels and the time evolution of the fire
  //behavior.  The output currently only contains the former.  Fire history wil be added.
  BurnupSim output = BurnupFM(fm, duffLoading, duffMoisture, tempAir, windSpeed, fireIntensity,
                              t_r, dT, nTimeSteps);

  return output;
}

/** Simulate ground fire returning the burn depth.
 * 
 * This effectively replaces the functionality of WildFire::getBurnOrgSoilthick() in the original
 * wildfire implementation.
 *
 * This is currently a stub and a place to work out how this simulation phase connects to the other
 * phases of fire.  The ground fire model is currently under development elsewhere.
 *
 * Inputs (proposed):
 * - The structure and state of the soil column.
 * - Energy inputs from aboveground fire components.
 * Burnup produces energy over time so it may be better to link the calculations?
 *
 * @returns The soil burn depth from ground fire (meters).
 */
double SimulateGroundFire()
{
  BOOST_LOG_SEV(glg, debug) << "Entering SimulateGroundFire()... [Stub]";

  double burnDepth = 0;
  
  //Calculate if the energy output is sufficiency to ignite the surface layer.
  //If not record that ignition failed and return.
  
  //Otherwise continue to calculate progressive smoldering downward.
  //This is the same problem of drying and heating to combustion as we move down.  However, we can
  //safely assume that the fire will not continue if we reach mineral soil, bedrock, permafrost, or
  //the water table.
  
  return burnDepth;
}
