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
  
  
  //Gather weather conditions:
  //Temp, humidity, wind, slope, ...
  
  
  //Fuel moisture must be calculated.
  //The 1 and 10 hour values could probably be estimated from a daily value, though we have to
  //assume a time of day.  We can superimpose a diurnal cycle.  The longer equilibrating fuels will
  //need some way of estimating the humidity value for the recent past.
  
  
  //Feed fuels and weather conditions into surface fire models:
  //First to Rothermel & Albini
  
  //Use the fire energy flux into the soil along with the fuel model as input into Burnup:
  
  //Use the energy flux into the soil (RA + Burnup) of the surface fire and the fire air temp
  //(Burnup fire environmental temperature?) as input to the ground fire model:
  //...
  

}
