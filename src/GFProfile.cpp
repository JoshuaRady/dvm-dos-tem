/**
 * @file GFProfile.cpp
 * \author Joshua M. Rady
 * Woodwell Climate Research Center
 * \date 2025
 *
 * @brief This file defines an object used to represent a soil column during simulations of
 * smoldering ground fire.
 *
 */

#include <cmath>
#include <string>

#include "GFProfile.h"
#include "FireweedMessaging.h"
#include "FireweedStringUtils.h"
#include "FireweedUtils.h"

//Public Functions:---------------------------------------------------------------------------------

/** Constructor
 *
 * This constructor makes an 'blank' profile of the specified size.  The soil properties are set to
 * invalid values that should be replaced.
 */
GFProfile::GFProfile(const int numLayers)
{
	this->numLayers = numLayers;
	t_ig = ProfileUndef;
	
	//Size the elements and set to invalid values:
	thickness_cm.resize(numLayers, ProfileUndef);
	layerDepth.resize(numLayers, ProfileUndef);
	tempC.resize(numLayers, ProfileUndef);
	bulkDensity.resize(numLayers, ProfileUndef);
	inorganicPct.resize(numLayers, ProfileUndef);
	moistureContentPct.resize(numLayers, ProfileUndef);
	c_s.resize(numLayers, ProfileUndef);
	
	//Outputs:
	burnt.resize(numLayers, false);
	heatSink.resize(numLayers, 0.0);
	heatSource.resize(numLayers, 0.0);
}

/** Get the number of layers.
 * 
 * @returns The number of layers in the soil column.
 */
int GFProfile::NumLayers() const
{
	return numLayers;
}

/** Return the dry soil mass for a layer.
 * 
 * @param The index of the layer.
 * 
 * @returns The dry soil mass (kg/m^2/layer).
 */
double GFProfile::DrySoilMassKg(const int layerIndex) const
{
	if (layerIndex < 0 || layerIndex >= numLayers)
	{
		Stop("Invalid layer index: " + std::to_string(layerIndex));
	}

	return bulkDensity[layerIndex] * (thickness_cm[layerIndex] / 100);//kg/m^3 * (cm / 100 cm/m) / layer = kg/m^2/layer.
}

/** Interpolate the profile.
 *
 * Check the profile and expand it to have layers of equal thickness and interpolate the current
 * values for the new layers, if necessary.
 *
 * @param newLayerThickness The desired layer thickness (cm) to convert to.
 *
 * @returns Nothing.
 
 * @note Our interpolation process begins by nudging the properties of the original layers to the
 * nearest new layer, on the basis of the layer centers.  The profile is then interpolated from
 * there.  This approach preserves the original values exactly in the profile, but their depths
 * may shift a bit.  The alternative would be to interpolate everywhere.  This prioritise the depths
 * over the exact values, since all layers would be interpolated in this case.  It the new layers
 * are thin the difference between the two approaches should be small.
 *
 * @note An issue for our interpolation is that each layer has a thickness.  If we associate our
 * know values with the middle of the layer this causes some issues.  If the top or bottom layer are
 * thick we don't have any information to interpolate on the outsides.  We fill new layers at the
 * top and bottom with the properties of original top and bottom layers.
 */
void GFProfile::Interpolate(const double newLayerThickness)
{
	//Check that newLayerThickness is some small round value?

	//At the end we will have a set of new layers.  The bottom of the deepest layer should be pretty
	//close to that of the current profile but will be a multiple of the new newLayerThickness:
	int initBottomIndex = numLayers - 1;//The position (index) of the bottom layer before interpolation.
	//This should be used carefully since the numbers of layers will change during this function.

	//Depth calculations;
	double columnBottom = layerDepth[numLayers - 1] + thickness_cm[numLayers - 1];//The current bottom of the column.
	double newColumnBottom = round(columnBottom / newLayerThickness) * newLayerThickness;//The adjusted bottom of the column.
	double newBottomLayerDepth = newColumnBottom - newLayerThickness;

	//Check if interpolation is necessary:
	bool needed = false;
	for (int i = 0; i < numLayers; i++)
	{
		if (thickness_cm[i] != newLayerThickness)//We don't give any wiggle room.
		{
			needed = true;
			break;
		}
	}

	if (needed)
	{
		//The first should be at zero, otherwise there are problems.
		if (layerDepth[0] != 0.0)
		{
			Stop("The top layer is not at the surface!");
		}

		//Nudge the existing layers to the nearest equidistant layer postion:
		//The most we will move a value is by newLayerThickness / 2, which will probably not matter
		//if newLayerThickness is small.  This process may mean that the top and bottom layers are
		//no longer at the top and bottom (see notes).
		for (int j = 0; j < numLayers; j++)
		{
			double layerCenter = layerDepth[j] + (thickness_cm[j] / 2);//The center depth of the original layer.
			
			//Nudge the original layer properties to the new layer that the original center falls in:
			int newLayerIndex = layerCenter / newLayerThickness;//Round down / truncate the division to get the layer.
			layerDepth[j] = newLayerIndex * newLayerThickness;
			thickness_cm[j] = newLayerThickness;
		}

		//Start the interpolation from the bottom of the column.  By doing this any layers added
		//will be below the next layer we have to examine.  That means the layer indexes won't
		//change as we are going, which makes things less confusing.

		//Add layers below the bottom layer if needed:
		//If the bottom of the lowest layer after nudging is still where we expect the new bottom it
		//is not necessary.
		if (layerDepth[numLayers - 1] != newBottomLayerDepth)
		{
			int numNewLayers = (newBottomLayerDepth - layerDepth[numLayers - 1]) / newLayerThickness;
			AddLayersToBottom(numNewLayers);
		}

		//Starting from the original bottom layer see if layers need to interpolated between them:
		for (int j = initBottomIndex; j > 0; j--)
		{
			//See if adjacent layers are one new layer thickness apart:
			double deltaZ = layerDepth[j] - layerDepth[j - 1];
			if (deltaZ != newLayerThickness)
			{
				//If not add layers to fill in between them:
				int numNewLayers = (deltaZ / newLayerThickness) - 1.0;
				//The original layers should have been nudged to regular positions and the distance
				//between them should be some exact multiple of the new layer thickness.  Check to
				//be sure:
				if ((numNewLayers + 1) * newLayerThickness != deltaZ)
				{
					Stop("Distance between layers is not a multiple of the new layer thickness.");
				}
			
				//Interpolate values for the new layers being inserted::
				InterpolateBetween(j - 1, j, numNewLayers);
			}
		}

		//Add layers above the top layer if needed:
		if (layerDepth[0] != 0.0)
		{
			int numNewLayers = layerDepth[0] / newLayerThickness;
			AddLayersToTop(numNewLayers);
		}

		//Perform quality control:
		//Check spacing?
		if (!Validate())
		{
			Print(std::cout);
		}
	}
}

/** Interpolate new layers between two soil layers.
 *
 * The function adds new layers between the layers specified and interpolates the soil physical
 * properties evenly across them.  If there is space between the parent layers then the new layers
 * will be given equal thicknesses to fill the gap. Otherwise they will be left blank and a warning
 * will be given. The interpolation but does not change the parent layer thicknesses.
 *
 * @param topLayer The top layer index to interpolate from.
 * @param bottomLayer The bottom layer index to interpolate to.
 * @param numNewLayers The number of layers to interpolate between them.
 *
 * @note Specifying the top and bottom layers for interpolation makes this function call clear to
 * read but specifying the bottom is not really necessary since it only makes sense to interpolate
 * between adjacent layers.  We could allow writing over intervening layers but the utility of that
 * is not clear.
 *
 * @returns Nothing,
 */
void GFProfile::InterpolateBetween(const int topLayer, int bottomLayer, const int numNewLayers)
{
	double newLayerThickness = ProfileUndef;

	if (bottomLayer != (topLayer + 1))
	{
		Stop("InterpolateBetween() can only interpolate between adjacent layers.");
	}

	//Insert new layers for the interpolation.
	InsertBlankLayers(bottomLayer, numNewLayers);
	int firstNewLayer = bottomLayer;//The position of the top new layer.
	bottomLayer = bottomLayer + numNewLayers;//The bottom layer moves down.

	//Check the parent layer spacing:
	//double zBetween = layerDepth[topLayer] + thickness_cm[topLayer] - layerDepth[bottomLayer];
	double zBetween = layerDepth[bottomLayer] - (layerDepth[topLayer] + thickness_cm[topLayer]);
	if (zBetween == 0.0)//The layers are touching.
	{
		Warning("The parent layers are touching, leaving interpolated layer thicknesses and depths blank.");
	}
	else if (zBetween > 0.0)//There is space between layers.
	{
		newLayerThickness = zBetween / numNewLayers;
	}
	else//Negative values indicate an error in the profile.
	{
		Stop("Layer spacing doesn't make sense.");
	}

	//Interpolate the properties:
	//Calculate the interpolation steps as fractional distances between the two values:
	double stepSize = 1.0 / (numNewLayers + 1);
	double step = stepSize;

	for (int i = 0; i < numNewLayers; i++)
	{
		thickness_cm[firstNewLayer + i] = newLayerThickness;
		layerDepth[firstNewLayer + i] = layerDepth[firstNewLayer + i - 1] + thickness_cm[firstNewLayer + i];
		
		tempC[firstNewLayer + i] = std::lerp(tempC[topLayer], tempC[bottomLayer], step);
		bulkDensity[firstNewLayer + i] = std::lerp(bulkDensity[topLayer],
		                                           bulkDensity[bottomLayer], step);
		inorganicPct[firstNewLayer + i] = std::lerp(inorganicPct[topLayer],
		                                            inorganicPct[bottomLayer], step);
		moistureContentPct[firstNewLayer + i] = std::lerp(moistureContentPct[topLayer],
		                                                  moistureContentPct[bottomLayer], step);
		c_s[firstNewLayer + i] = std::lerp(c_s[topLayer], c_s[bottomLayer], step);

		step += stepSize;
	}
	//Leave burnt, heatSink, and heatSource at their defaults.
}

/** Add layers to the top of the profile, cloning the properties of the current top layer.
 *
 * This function is used when the surface layer is being sliced into thinner layers and the code
 * makes assumptions consistant with this.  It would also probably be safe to assume that after this
 * the new top layer will be at the surface, but we don't currently check for that here.
 *
 * @param numNewLayers Number of layers to add.
 *
 * @returns Nothing,
 */
void GFProfile::AddLayersToTop(const int numNewLayers)
{
	if (numNewLayers < 1)
	{
		Stop("AddLayersToTop(): Invalid value for numNewLayers: " + std::to_string(numNewLayers));
	}

	double layerThickness = thickness_cm[0];

	thickness_cm.insert(thickness_cm.begin(), numNewLayers, layerThickness);
	
	for (int i = numNewLayers - 1; i >= 0; i--)
	{
		layerDepth.insert(layerDepth.begin(), layerThickness * i);
	}

	tempC.insert(tempC.begin(), numNewLayers, tempC[0]);
	bulkDensity.insert(bulkDensity.begin(), numNewLayers, bulkDensity[0]);
	inorganicPct.insert(inorganicPct.begin(), numNewLayers, inorganicPct[0]);
	moistureContentPct.insert(moistureContentPct.begin(), numNewLayers, moistureContentPct[0]);
	c_s.insert(c_s.begin(), numNewLayers, c_s[0]);

	//The profile should only be modified before a simulation so the outputs are at initial values:
	burnt.insert(burnt.begin(), numNewLayers, false);
	heatSink.insert(heatSink.begin(), numNewLayers, 0.0);
	heatSource.insert(heatSource.begin(), numNewLayers, 0.0);

	numLayers = numLayers + numNewLayers;
}

/** Add layers to the bottom of the profile, cloning the properties of the current bottom layer.
 *
 * This function is used when the bottom layer is being sliced into thinner layers and the code
 * makes assumptions consistant with this.
 *
 * @param numNewLayers Number of layers to add.
 *
 * @returns Nothing,
*/
void GFProfile::AddLayersToBottom(const int numNewLayers)
{
	if (numNewLayers < 1)
	{
		Stop("AddLayersToBottom(): Invalid value for numNewLayers: " + std::to_string(numNewLayers));
	}

	thickness_cm.insert(thickness_cm.end(), numNewLayers, thickness_cm[numLayers - 1]);
	
	//for (int i = 0; i < numLayers; i++)
	for (int i = 0; i < numNewLayers; i++)
	{
		layerDepth.push_back(layerDepth[numLayers - 1 + i] + thickness_cm[numLayers - 1]);
	}

	tempC.insert(tempC.end(), numNewLayers, tempC[numLayers - 1]);
	bulkDensity.insert(bulkDensity.end(), numNewLayers, bulkDensity[numLayers - 1]);
	inorganicPct.insert(inorganicPct.end(), numNewLayers, inorganicPct[numLayers - 1]);
	moistureContentPct.insert(moistureContentPct.end(), numNewLayers, moistureContentPct[numLayers - 1]);
	c_s.insert(c_s.end(), numNewLayers, c_s[numLayers - 1]);

	//The profile should only be modified before a simulation so the outputs are at initial values:
	burnt.insert(burnt.end(), numNewLayers, false);
	heatSink.insert(heatSink.end(), numNewLayers, 0.0);
	heatSource.insert(heatSource.end(), numNewLayers, 0.0);

	numLayers += numNewLayers;
}

/** Insert new blank layers into the profile. (Properties must be added afterwards.)
 *
 * @param insertAt The index to insert at (i.e. the layer currently at index insertAt will be shifted down).
 * @param numNewLayers Number of layers to insert.
 *
 * @returns Nothing,
 */
void GFProfile::InsertBlankLayers(const int insertAt, const int numNewLayers)
{
	if (numNewLayers < 1)
	{
		Stop("InsertBlankLayers(): Invalid value for numNewLayers: " + std::to_string(numNewLayers));
	}

	thickness_cm.insert(thickness_cm.begin() + insertAt, numNewLayers, ProfileUndef);//Or add values?
	layerDepth.insert(layerDepth.begin() + insertAt, numNewLayers, ProfileUndef);//Or add values?

	//Insert and set to invalid values:
	tempC.insert(tempC.begin() + insertAt, numNewLayers, ProfileUndef);
	bulkDensity.insert(bulkDensity.begin() + insertAt, numNewLayers, ProfileUndef);
	inorganicPct.insert(inorganicPct.begin() + insertAt, numNewLayers, ProfileUndef);
	moistureContentPct.insert(moistureContentPct.begin() + insertAt, numNewLayers, ProfileUndef);
	c_s.insert(c_s.begin() + insertAt, numNewLayers, ProfileUndef);

	//The profile should only be modified before a simulation so the outputs are at initial values:
	burnt.insert(burnt.begin() + insertAt, numNewLayers, false);
	heatSink.insert(heatSink.begin() + insertAt, numNewLayers, 0.0);
	heatSource.insert(heatSource.begin() + insertAt, numNewLayers, 0.0);

	numLayers += numNewLayers;
}

/** Perform some checks on the profile to see if looks valid.
 *
 * @returns Whether the profile passed the checks.
 */
bool GFProfile::Validate() const
{
	bool valid = true;

	if (numLayers < 1)
	{
		Warning("Invalid value of numLayer.");
		valid = false;
	}

	//Make sure the number of property elements are consistant:
	if (numLayers != thickness_cm.size() ||
	    !SameLengths(thickness_cm, layerDepth, tempC, bulkDensity, inorganicPct, moistureContentPct,
	                 c_s, burnt, heatSink, heatSource))
	{
		Warning("Some profile elements are not the same length.");
		valid = false;
	}

	//Confirm the depth profile is continuous:
	for (int i = 0; i < numLayers; i++)
	{
		if (thickness_cm[i] <= 0)
		{
			Warning("Invalid thickness for layer " + std::to_string(i));
			valid = false;
		}

		if (i == 0)
		{
			if (layerDepth[i] != 0.0)
			{
				Warning("The top layer is not at the surface!");
				valid = false;
			}
		}
		else
		{
			if (layerDepth[i] != (layerDepth[i - 1] + thickness_cm[i - 1]))
			{
				Warning("Layer dimensions not consistant for layer " + std::to_string(i));
				valid = false;
			}
		}
	}

	//Check for invalid values:
	if (!InRange(tempC, -60.0, 50.0))//A broad guess range for non-fire temperatures.
	{
		Warning("Unlikely temperature value(s).");
		valid = false;
	}

	//Histisols can have low bulk densities approaching zero for some peats.  Compacted glacial till
	//on the other hand approaches that of concrete to we need to give a pretty wide range.  Values
	//pulled from figure 4.44 in:
	//Brady, N.C. and Weil, R.R., 2016.
	//The nature and properties of soils (Vol. 15, pg. 163). Columbus, OH: Pearson.
	if (!InRange(bulkDensity, 0.0, 2400.0))
	{
		Warning("Unlikely bulk density value(s).");
	}

	if (!InRange(inorganicPct, 0.0, 100.0))
	{
		Warning("Invalid percentage(s) for inorganic content");
		valid = false;
	}

	if (!InRange(moistureContentPct, 0.0, 1000.0))//Don't have a great idea of the upper limit.
	{
		Warning("Invalid soil moisture content(s).");
		valid = false;
	}

	//We don't check the calculated properties.

	return valid;
}

/** Print the soil profile data to an output stream.
 *
 * @param output The output stream to print to.
 *
 * @returns The ostream so it can be conatinated to.
 */
std::ostream& GFProfile::Print(std::ostream& output) const
{
	output << "Ground fire soil profile:" << std::endl;
	output << "The column has " << numLayers << " layers." << std::endl;
	output << "The soil ignition temperature is " << t_ig << " C." << std::endl;

	//output << "Layer thickness (cm): " << thickness_cm << std::endl;
	output << "Layer thickness (cm): ";
	PrintVector(output, thickness_cm);
	output << "Layer depth (cm to top): ";
	PrintVector(output, layerDepth);
	output << "Layer temperature (C): ";
	PrintVector(output, tempC);
	output << "Layer bulk density (kg/m^3): ";
	PrintVector(output, bulkDensity);
	
	output << "Layer dry mass (kg/m^2/layer): ";
	for (int i = 0; i < numLayers - 1; i++)
	{
		output << DrySoilMassKg(i) << ", ";
	}
	output << DrySoilMassKg(numLayers - 1) << std::endl;

	output << "Layer % inorganic: ";
	PrintVector(output, inorganicPct);
	output << "Layer % moisture content: ";
	PrintVector(output, moistureContentPct);
	output << "Layer heat capacity (kJ/kg/K, AKA c_s): ";
	PrintVector(output, c_s);

	output << "Layer burnt: ";
	output << std::boolalpha;//Print as burnt as strings.
	//PrintVector(output, burnt);
	for (int i = 0; i < burnt.size() - 1; i++)
	{
		output << burnt[i] << ", ";
	}
	output << burnt[burnt.size() - 1] << std::endl;
	output << std::noboolalpha;
	
	output << "Layer heat sink: ";
	PrintVector(output, heatSink);
	output << "Layer heat source: ";
	PrintVector(output, heatSource);

	return output;
}

/** Print the soil profile to an output stream as a set of delimited data rows for each layer..
 *
 * @param output The output stream to print to.
 * @param delim The delimiter character.  Defaults to the tab character.
 *
 * @returns The ostream so it can be conatinated to.
 */
std::ostream& GFProfile::PrintDelimited(std::ostream& output, const char delim) const
{
	//Print the header:
	output << "LayerThickness" << delim << "LayerDepth" << delim << "TempC" << delim << "t_ig" <<
	          delim << "BulkDensity" << delim << "InorganicPct" << delim << "MoistureContentPct" <<
	          delim << "c_s" << delim << "Burnt" << delim << "HeatSink" << delim << "HeatSource" <<
	          std::endl;//Or newline?
	//Omit the soil mass?

	//Print the layer data in rows:
	output << std::boolalpha;//Print as burnt as strings.
	for (int i = 0; i < numLayers; i++)
	{
		//Consider rounding values...
		output << thickness_cm[i] << delim << layerDepth[i] << delim << tempC[i] << delim << t_ig <<
		          delim << bulkDensity[i] << delim << inorganicPct[i] << delim <<
		          moistureContentPct[i] << delim << c_s[i] << delim << burnt[i] << delim <<
		          heatSink[i] << delim << heatSource[i] << delim << std::endl;//Or newline?
	}
	output << std::noboolalpha;

	return output;
}

//External functions:-------------------------------------------------------------------------------

/* Overloaded stream print operator for GFProfile.
 *
 */
std::ostream& operator<<(std::ostream& output, const GFProfile& gfProfile)
{
	gfProfile.Print(output);
	return output;
}
