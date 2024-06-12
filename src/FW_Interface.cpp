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

/** Calculate wildfire behavior and effects using the modeled vegetation, fuels, and meteorology,
 *
The function takes the ecosystem state (cohort pointer? WildFire pointer) and meteorology...

Output:
The output should be stored in the cohort object...
The soil burn depth should be returned/recorded to take advantage of the existing code.

 */
void RevisedFire(WildFire* wf)//The name will definitely change.
{

  //Determine the fuel model matching the location's CMT:
  //fm = GetMatchingFuelModel(CMTnumber)

  //Determine the surface fuels from the model vegetation and soil states:
  //The fuels must be estimated from:
  // - Live fuels: gramanoid, herbaceous, shrub PFTs, and the moss layer.
  // - Litter: The first soil layer dead vegetation component (rawc) and woody debris (wdebrisc),
  //which currently is only present after a previous fire.
  //The mass is upscaled from the carbon stocks.  The litter is distributed into the SAV bins based
  //on guesstimates informed the PFT inputs and by some ideas of breakdown rates.
  //Previous or current cfall could be used to help estimate this but some assumptions still need
  //to be made.
  
  //Update the fuel masses from their default values:
  //fm.w_o_ij[0] = ...
  
  //Gather weather conditions:
  //Temp, humidity, wind, slope, ...
  double tempAir = wf->edall.d_atms.ta;//Daily air temp (at surface).
  //double humidity = ?????
  //Wind speed (m/min) at midflame height is approximated as 2m wind speed.  We need a sub-daily value ~ daily value?
  //double windSpeed = ?????;//Wind speed does not seem to be available.
  //The wildfire object contains the slope.
  double slope = wf->slope;//What at the units?  May have to convert to fractional slope.
  
  //Fuel moisture must be calculated.
  //The 1 and 10 hour values could probably be estimated from a daily value, though we have to
  //assume a time of day.  We can superimpose a diurnal cycle.  The longer equilibrating fuels will
  //need some way of estimating the humidity value for the recent past.
  //Days since last rain is available in wf->edall.d_atms.
  
  
  //Feed fuels and weather conditions into surface fire models:
  //First to Rothermel & Albini
  
  //Use the fire energy flux into the soil along with the fuel model as input into Burnup:
  
  //Use the energy flux from the aboveground fire into the soil surface (RA + Burnup + crown) and
  //the fire air temp (Burnup fire environmental temperature?) as input to the ground fire model:
  //...
  
  //Also pass in the soil profile conditions...
  
  //SimulateGroundFire();
  //wf->fd->fire_soid.burnthick = burnDepth;
  //Include messaging of original code?????

}

/* Determine the fuel model matching a given community type:*/
FuelModel GetMatchingFuelModel(int cmt)
{
  //Get the number of the fuel model from the crosswalk in the parameter files.
  //This crosswalk needs to be made!!!!!
  int fuelModelNumber = 1;//Temporarily hardwired
  
  //Get the fuel model data:
  //fuelModel = GetFuelModel(fuelModelNumber);
  
  return fuelModel;
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
