# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # DEM coregistration and glacier geodetic mass balance with xdem
#
# ## Advanced exercise
#
# <div class="alert alert-info" style="font-size:110%">
#  <h3>Objectives</h3>
#  
# This notebook showcases how to use xdem to:
# - coregister two DEMs
# - visualize biases and uncertainties in the elevation change map
# - calculate the glacier geodetic mass balance.
#   
# For this, we rely on two DEMs derived from Pleiades stereo images over Mera glacier, Nepal.
# </div>

# ### Import the necessary modules

# +
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib widget

import geoutils as gu
import xdem

# To avoid interpolation in plt.imshow
plt.rcParams['image.interpolation'] = 'none'
# -


# ## 1 - Load input data
#
# <div class="alert alert-info" style="font-size:110%">
# All the data is taken from the Mera mass balance geodetic tutorial located in the "mb_Mera" folder. We will use
# <ul>
#     <li> 2012 Pleiades DEM: Mera_Pleiades_2012-11-25_DEM_4m.tif
#     <li> 2018 Pleiades DEM: Mera_Pleiades_2018-10-28_DEM_4m.tif
#     <li> Copernicus (reference) DEM: Mera_COP30_DEM_UTM45_data.tif
#     <li> RGI outlines of all glaciers in the region: Glacier_inventory_around_Mera.shp
#     <li> Mera outline of 2012: Mera_outline_2012_realigned.shp
#     <li> Mera outlines of 2018: Mera_outline_2018_realigned.shp
# </ul>
# </div>
#

# #### Load all DEMs
#
# All DEMs must be cropped to a common extent and reprojected onto the 2018 DEM grid. \
# Note: you may use the function `gu.raster.load_multiple_rasters`.




# #### Load glacier outlines: RGI outlines, Mera 2012 and Mera 2018 outlines



# ### Calculate and plot elevation change

# #### Calculate the 2012-2018 elevation change




# #### Plot original dh map with glacier contours



# There are non-zeros elevation changes outside of glaciers => These DEMs need to be coregistered.

# ## 2 - DEM coregistration

# #### Prepare inputs for coregistration
# First we create a mask, i.e. a raster of same shape as our dh map, to mask pixels on glaciers. `gl_mask` is `True` on glaciers, `False` elsewhere. \
# Since the RGI outlines can be slightly inaccurate, it is recommended to use a 100 m buffer.



# Then we mask pixels in steep slopes (> 40 degrees) and gross blunders (abs(dh) > 50).



# We plot the final mask of pixels to be used for coregistration.



# ### Run coregistration
# #### We use the Nuth & Kaab (2011) algorithm to estimate a horizontal offset, then we remove a possible vertical bias by removing the median dh value in stable terrain.

# First, we need to define the coregistration function. \
# Then we estimate the coregistration needed between our two Pleiades DEMs. \
# Finally, we apply that coregistration to the 2012 DEM.



# #### Print the output of the coregistration: offset in east and north direction, vertical offset.



# ### <span style='color:red '> **Question:** </span> How many iterations of the Nuth & Kaab algorithm were run?
#

# #### Answer: ...

# ### Calculate new elevation change and plot






# #### We see that elevation changes outside glaciers are close to zero (yellow color).

# ### Calculate statistics of before/after coregistration
#
# Print the number of inlier pixels, mean, median and NMAD of elevation change before and after coregistration.
#
# Because `dh.data` is a masked array, we use `compressed()` to output only unmasked values.



# ### Plot statistics of median elevation change as a function of slope, aspect and reference elevation

# First, we calculate slope, aspect and reference elevation for pixels in inlier mask.





# We calculate elevation change before and after coregistration in inlier mask.




# We use the function `xdem.spatialstats.nd_binning` to calculate statistics (default is count, median, NMAD) for bins of each category, for elevation change before coregistration.




# The results are pandas dataframes containing the statistics for each bin



# We do the same for elevation change after coregistration.




# And finally we make plots summarizing it for before/after coregistration.






# ### BONUS: Calculate statistics as a function of curvature and along-track coordinates



# ## 3 - Geodetic mass balance of Mera glacier
#
# ### Calculate mean dh over glaciers, per elevation bins
# We use the function `xdem.volume.hypsometric_binning` (see brief documentation [here](https://xdem.readthedocs.io/en/stable/gen_modules/xdem.volume.html)) to calculate the mean (or median) value of dh within elevation bands. \
# The calculation requires a (raster) mask of the glacier, the dh values, and a reference DEM used to get the elevation of the pixels. All must be on the same grid. \
# The default elevation bins have a height of 50 m, from the DEM minimum to maximum. Here we set the bin height to 25 m with the argument `bins`.



# The output is a panda's DataFrame containing the elevation bin edges, the mean dh ('value') and observations (valid pixels) count. Print it to screen.



# ### Calculate the glacier area within each elevation bin
# We use the function `xdem.volume.calculate_hypsometry_area`. \
# This is particularly needed for data with gaps, as the values in `ddem_bins['count']` are the number of pixels with observations, not total pixel count. \
# The result is in m$^2$.



# ### <span style='color:red '> **Questions:** </span> 
# - what is the median dh at 5400 m ?
# - at which elevation is the highest thinning?
# - what is the area of glaciers in the elevation band at 6000 m?
#

# ### Plot the results
# 1. Plot glacier area vs elevation
# 2. Plot the elevation change vs elevation

# The glacier area is plotted with blue bars (x axis at the bottom) and elevation change in dark line (x axis at the top) versus elevation (y axis).



# ### Calculate total volume and mass change
# - volume change is in m$^3$
# - volume is converted to mass (kg) assuming a volume to mass conversion factor of 850 kg/m$^3$ (Huss, 2013, TC)



# ### Calculate specific mass balance
# In meter water equivalent (m w.e.)
# We divide by the average area between the two dates.



# ### Print results



# ### BONUS - Use geopandas' explore and/or folium to plot the results


# <div class="alert alert-info" style="font-size:110%">
#
# - create an interactive map of the RGI outlines with `rgi_outlines.ds.explore`. Use `style_kwds` to plot empty polygons with black edges. You can set a different background with argument `tiles` (list of options available [here](https://leaflet-extras.github.io/leaflet-providers/preview/)).
# - convert the dh Raster into a numpy array with Nan, projected into web mercator (ESRI:53004). I suggest downsampling to 20 m for speed.
# - add a Folium image layer with `folium.raster_layers.ImageOverlay`. You need to set a colormap.
# - optionally, try displaying the colorbar using information from [this link](https://leafmap.org/notebooks/62_folium_colorbar/)
# - save the output as an html file for later use
# - Play with the different options to customize your plot
#
# </div>


