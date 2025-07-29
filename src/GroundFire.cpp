/**
 * @file GroundFire.cpp
 * \author Joshua M. Rady
 * Woodwell Climate Research Center
 * \date 2025
 *
 * @brief This file defines functions to simulate smoldering ground fire using a moderate complexity
 * model where combustion is predicted layer by layer.
 *
 */

#include <algorithm>//For fill().
#include <cmath>
#include <numeric>
#include <string>

#include ""../include/GroundFire.h"//"GroundFire.h"//Note: Difference from source repo.
#include "FireweedMessaging.h"
#include "FireweedUtils.h"

/** Calculate and return the heat of combustion of an organic soil based on its organic content.
 *
 * This is a simple calculation that may not be robust for all organic soils.  It assumes that the
 * inorganic content of the soil has no effect on the energy released during combustion and that the
 * organic material has a fixed heat of combustion.  The former may be somewhat true but the later is
 * certainly not.  The energy stored in organic material differs by source and level of decomposition.
 * This function is biased towards rich organic soils, like peats, where we hope these assmumptions
 * will be less problematic.  The organicHOC parameter can also be changed to account for differences
 * in organic quality if information is available.
 * 
 * Additionally the heat of combustion is generally considered for the case of combustion with excess
 * oxygen.  For flaming combustion this may be appropriate but for smoldering combustion can occur by
 * more than one pathway, with different energies.
 * 
 * The smoldering specific formulation includes the component oxidations:
 * ΔH_c = ΔH peat beta char oxidation + ΔH alpha char oxidation + ΔH beta char oxidation
 * However, we only have those from a few papers so we can only make a first order estimate.  We may
 * be able to improve this function over time with more research.
 * 
 * @param inorganicPct The percent organic content of the soil on the dry basis.
 * @param  organicHOC The (theoretical) heat of combustion for a purely organic soil (MJ/kg).
 * The peat moss in Frandsen 1987 has a heat of combustion of 20.6 MJ/kg at 3.7% inorganic.  We
 * can estimate that of 100% organic as:
 * 20.6 / (100 - 3.7) * 100 = 21.39148 MJ/kg
 *
 * @returns The soil heat of combustion (MJ/kg, ~net fuel heat content).
 */
double SoilHeatOfCombustion(const double inorganicPct, const double organicHOC)
{	
	double soilHOC = organicHOC * (100.0 - inorganicPct) / 100.0;//MJ/kg
	return soilHOC;
}

/** Calculate and return the heat sink for a mass of organic soil.  This is the is the sum of the
 * energy required for drying and to further heat the dry soil to ignition.
 *
 * @param drySoilMass = The mass of dry soil (kg)
 * @param moistureContentPct = Percent moisture content (water mass / dry soil mass * 100%).
 * @param T_int Intial soil temperature (C).
 * @param T_ig Temperature of (auto/self-)ignition or combustion (C).
 * @param c_s Specific heat capacity of the soil (kJ/kg/K). [at 20C?]
 *        This is a function of the soil component heat capacities and the temperature, though for
 *        now  we ignore the temperature dependence as is frequently done.
 *
 * @returns The heat sink (kJ/kg moist soil).
 */
double SoilHeatSink(const double moistureContentPct, const double T_int, const double T_ig, const double c_s)
{
	//Calculate the heat of drying (kJ/kg water):
	//Here we use degrees C but could use K.
	
	//Properties of water (should be moved to constants):
	double DelatH_vap = 2257.0;//ΔH_vap = (Latent) heat of vaporization of water (J/g or kJ/kg)
	double c_w = 4.184;//The specific heat capacity of water (kJ/kg/K at 20C).
	
	//Soil water is in free and bound forms.  The bound water fraction requires additional energy to
	//evaporate.  I'm still trying to figure out how to estimate what the bound fraction is.  There may
	//even be more than one bound fraction as hygroscopic and tightly bound water may be different?
	//As a start we could exclude free water above ~100 - 105% moisture content?
	//These terms are being set to zero until I figure this out.
	double frac_bound = 0.0;//The fraction of water in a bound form.
	double DeltaH_des = 0.0;//ΔH_des = Heat of desorption of water (kJ/kg)
	
	double DeltaH_dry = (100.0 - T_int) * c_w + DelatH_vap + (frac_bound * DeltaH_des);
	
	//Heat of ignition of dry soil (kJ/kg):
	double DeltaH_ig = (T_ig - T_int) * c_s;
	
	//Total heat required, (kJ/kg moist soil):
	//For 1 kg dry soil the soil water (kg) = 1 kg * percent moisture content / 100%.
	double h_sink = (DeltaH_dry * moistureContentPct / 100.0) + DeltaH_ig;
	return h_sink;
}

/** Calculate the distribution of heat over a number of layers, decreasing linearly with depth.
 *
 * The distribution is treated as triangular with the fraction of the distribution in a layer given
 * by the area of a slice of that triangle.  We assume that we are distributing over full layers.
 * The function returns fractional weights adding to one.  The weight order is top to bottom.  The
 * actual heat per layer can be obtained by multiplication.
 *
 * @param numLayers The number of layer for which to compute the distribution.
 *
 * @returns A vector of fractional weights adding to one for each layer from top to bottom.
 */
std::vector <double> HeatDistributionLinear(const int numLayers)
{
	if (numLayers < 1)
	{
		Stop("The number of layer must be at least 1.");
	}
		std::vector <double> layerWeights(numLayers);
		double layerThickFrac = 1.0 / numLayers;//The relative thickness of layers as a fraction of the whole depth.

	for (int i = 0; i < numLayers; i++)
	{
		//The distribution is a right triangle so the width can be defined by a line of form m*x + b,
		//or m*y + b since it is facing down.  We don't know anything about were the layers start, end,
		//or how deep they are.  We treat everything as relative. For convenience we can make the
		//dimensions of the triangle 2 x 1 so it has an area of 1.   It follows that m = 2.  Since we
		//are starting from the top we have to set the intercept b = 2.  Each layer slice has a depth
		//of 1 / numLayers.
		double m = -2.0;
		double b = 2.0;
		
		//The depths of the top and bottom of the slice:
		double d_top = (i - 1.0) * layerThickFrac;
		double d_bottom = i * layerThickFrac;
		
		//The area for a layer from depth d1 to d2 = (d2 - d1) * (w1 + w2)/2, where w is the width of
		//the triangular distribution at depth d.  d2 - d1 always = layerThickFrac:
		double w_top = m * d_top + b;
		double w_bottom = m * d_bottom + b;
		
		layerWeights[i] = layerThickFrac * (w_top + w_bottom) / 2.0;
	}

	if (!FloatCompare(std::accumulate(layerWeights.begin(), layerWeights.end(), 0.0), 1.0))
	{
		Stop("Weights do not sum to 1.");
	}

	return layerWeights;
}

/** Calculate the depth of burn for a generic soil profile.
 * 
 * Ground fire is simulated as smoldering combustion that moves progressively in the downward
 * direction.  The soil column is represented as layers of equal thickness.  The soil column is
 * treated as homogeneous in the vertical direction and vertical spread is not considered.
 *
 * Ground fire initiates when energy input from the aboveground fire into the soil drys and heats
 * the top layer to the temperature of ignition.  Currently the top layer is not treated in a
 * special way.  However, given the importance of ignition we may change this in the future to
 * simulate the surface fire input more explicitly in time and with a full heat transfer
 * representation.
 *
 * The fire extinguishes when the energy to dry and heat the next layer to combustion is less than
 * that transfered from combustion above.  The model uses simplified representations of the process
 * of heat transfer and loss.  The calculation occurs layer by layer without the timing of
 * progression simulated directly, i.e. the model used layer-steps rather than time-steps.
 *
 * This function preserves the model iteration 17 (Proj_11_Exp_17_Analysis.r DominoGroundFire17()).
 *
 * @param soilCol A GFProfile object representing the soil profile for the simulation.
 * @param fireHeatInput Total longwave heat input into the soil from the surface fire (kJ/m^3).
 * @param heatLossFactor The fraction of heat lost from the soil (per cm).  This intended to
 *        represent radiant and  convective losses from the soil surface in the time is take for a
 *        layer to combust. The downward efficiency factor from Benscoter et al. 2011 = downwardEF.
 *        As a heat loss = 0.83.
 * @param d_max The maximum depth of heat transfer for the heat transfer approximation (cm).
 *
 * @returns The burn depth (cm).
 */
double DominoGroundFire(GFProfile& soilCol, const double fireHeatInput,
                        const double heatLossFactor, const double d_max)
{
	//Check for layers of equal depth?

	double burnDepth = 0.0;//Return value.
	bool linearHeatTransfer = false;//Was a parameter.

	//The soil profile contains layer thicknesses for each layer.  That was added with the plan to
	//enable variable depth layers.  We haven't done that yet and the following code assumes that
	//it is constant.
	double layerThickness_cm = soilCol.thickness_cm[0];
	
	//Calculate the number of layers to distribute the heat from the surface fire to:
	//We currently transfer the heat to the top 1 cm.  This is a guesstimate and should probably be
	//changed to model parameter.  The heat is also spread evenly.  This should be probably be
	//updated as well.
	//The number of layers needs to be a whole number.  We round down but it might be better to
	//require that the layer thickness be an even divisor of 1.
	int numSurfLayers = 1.0 / layerThickness_cm;
	double heatPerSurfLayer = fireHeatInput / numSurfLayers;
	for (int i = 0; i < numSurfLayers; i++)
	{
		soilCol.heatSource[i] = heatPerSurfLayer;
	}

	//For each layer see if the heat from above (the fire or previous layer) is enough to ignite it:
	for (int i = 0; i < soilCol.NumLayers(); i++)
	{
		//The heat source for this layer was previously calculated and stored:
		double heatSourceToLayer = soilCol.heatSource[i];

		//The (total) heat sink for a layer is determined exclusively by its physical properties:
		soilCol.heatSink[i] = SoilHeatSink(soilCol.moistureContentPct[i], soilCol.tempC[i],
		                                   soilCol.t_ig, soilCol.c_s[i]);

		if (heatSourceToLayer >= soilCol.heatSink[i])
		{
			soilCol.burnt[i] = true;

			//Energy in excess of that needed to dry and heat this layer has to go somewhere.  Some
			//will be lost and some will heat lower layers.
			double layerExcessHeat = heatSourceToLayer - soilCol.heatSink[i];

			//Calculate the heat released from combustion for this layer:
			//This is approximately the higher heating value?????
			double layerHeatOfCombustion = SoilHeatOfCombustion(soilCol.inorganicPct[i]) *
				soilCol.DrySoilMassKg(i) * 1000.0;//MJ/kg -> kJ/kg

			//The total heat available from this layer:
			double totalLayerHeat = layerHeatOfCombustion + layerExcessHeat;

			//Compute the heat transferred to lower layer(s):
			double heatSourceDown = totalLayerHeat * (1.0 - heatLossFactor);//Or propagating heat?

			//Compute how the heat is distributed:
			//This setting could be expanded to set of values.
			if (!linearHeatTransfer)//(Pure) domino heat transfer scheme:
			{
				//In the simplest representation all the heat is transferred to the next layer:
				if (i != soilCol.NumLayers())//Don't walk off the bottom of the soil column.
				{
					soilCol.heatSource[i + 1] = heatSourceDown;
					//Need to preserve fire inputs if they go to more than one layer!!!!!!  This should be considered a bug?
				}
			}
			else//Simple linear heat transfer scheme:
			{
				if (i == soilCol.NumLayers())//If we have reached the bottom there is nowhere to transfer heat to:
				{
					break;
				}
		
				//Distribute the heat to layers below:
				int distLayers = round(d_max / layerThickness_cm);//The number of layers to distribute heat to.
				//The top and bottom layers to distribute heat to:
				int topLayer = i + 1;
				int bottomLayer = i + distLayers;
				
				std::vector <double> wts = HeatDistributionLinear(distLayers);
		
				for (int j = topLayer; j <= bottomLayer; j++)
				{
					double heatToLayerJ = wts[j - i] * heatSourceDown;
		
					if (heatToLayerJ < 0.0)
					{
						Stop("Negative heat flux for layer" + std::to_string(j));
					}
		
					if (j <= soilCol.NumLayers())//Don't walk off the bottom of the soil column.
					{
						//Accumulate the heat:
						soilCol.heatSource[j] = soilCol.heatSource[j] + heatToLayerJ;
					}
				}
			}
		}
		else
		{
			break;
		}
	}

	//To return the burn depth we need the depth at the bottom of the last layer burned:
	int lowest = -1;
	for (int k = soilCol.NumLayers(); k >= 0; k--)
	{
		if (soilCol.burnt[k])
		{
			lowest = k;
			break;
		}
	}
	
	if (lowest != -1)
	{
		burnDepth = soilCol.layerDepth[lowest] + layerThickness_cm;
	}

	//soilCol.Print(std::cout);//For debugging.

	return burnDepth;
}
