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
  FuelModel fm = GetMatchingFuelModel(CMTnumber)

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
  
  
  //De we calculate the fuel bed depth?
  
  
  //Gather weather conditions:
  double tempAir = wf->edall.d_atms.ta;//Daily air temp (at surface).
  //Curten (daily if not hourly) (relative) humidity is needed.  A recent time history of
  //double humidity = ?????
  //Wind speed (m/min) at midflame height is approximated as 2m wind speed.  We need a sub-daily value ~ daily value?
  //double windSpeed = ?????;//Wind speed does not seem to be available.
  //The wildfire object contains the slope.
  double slope = wf->slope;//What at the units?  May have to convert to fractional slope.
  //Shortwave radiation may be needed for the moisture calculations.
  
  
  //Fuel moisture must be calculated.
  FuelMoisture theFuelMoisture = CalculateDeadFuelMoisture()
  //fm.M_f_ij[0] = theFuelMoisture.XXXXX
  
  
  //Feed fuels and weather conditions into surface fire models:
  //First to Rothermel & Albini.  We use the calculations to get some component values.
  //This interface is under development.  It takes a fuel model (and attendant data) and returns the
  //calculation details.
  //FM = SAV_ij, w_o_ij, M_x_1, M_f_ij, h_ij, S_T_ij, S_e_ij, rho_p_ij, liveDead
  RAData raData = SpreadRateRothermelAlbini_Het(fm,
                                     fuelBedDepth?????,//fuelBedDepth
                                     windSpeed,//U
                                     slope,//slopeSteepness
                                     false,//useWindLimit
                                     Metric,//units
                                     FALSE)//debug
  
  //Use the fire energy flux into the soil along with the fuel model as input into Burnup:
  //raData.XXXXX
  
  //Use the energy flux from the aboveground fire into the soil surface (RA + Burnup + crown) and
  //the fire air temp (Burnup fire environmental temperature?) as input to the ground fire model:
  //...
  
  //Also pass in the soil profile conditions...
  
  //SimulateGroundFire();
  //wf->fd->fire_soid.burnthick = burnDepth;
  //Include messaging of original code?????

}

/* Determine the fuel model matching a given community type:

 * Note: FuelModel does not yet exist and is being developed elsewhere!
 */
FuelModel GetMatchingFuelModel(int cmt)
{
  //Get the number of the fuel model from the crosswalk in the parameter files.
  //This crosswalk needs to be made!!!!!
  int fuelModelNumber = 1;//Temporarily hardwired
  
  //Get the fuel model data:
  //fuelModel = GetFuelModel(fuelModelNumber);
  
  return fuelModel;
}

/* Calculate fuel moisture based on recent weather:
*
* This is a fire science dependent function and as such should probably be moved into the FW library.
* Here the emphasis is on determining the available inputs we can pass in.
* An alternative to calculating this at the time od fire would be to make fuel moisture a constantly
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

Output:
A fuel moisture object seems like a compact way to return the moisture values.
Note: FuelMoisture does not currently exist!

What about live fuel moisture?

 */
FuelMoisture CalculateDeadFuelMoisture()
{
  FuelMoisture dfmc;//FMC = dead fuel moisture content.
  
  //...
  
  return fmc;
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
