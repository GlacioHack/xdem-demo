# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Manipulating raster/vector data with geoutils - Exercise & solutions

# In this notebook, we apply the tools learnt in the `intro_geoutils` notebook to a different data set.
#
# We use the data from the Mera geodetic mass balance tutorials to be presented in the following days. They are available in the data folder "04_mb_Mera".\
# There are 2 Pleiades DEMs from November 2012 and October 2018.\
# There are manual outlines of Mera glacier for the same dates and a sample of RGI outlines for the study area.


# ## Import the necessary modules

# +
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib widget

import geoutils as gu
import xdem
# -

# ## Loading the data ###

# ### Load the two Pleiades DEMs from 2012 and 2018
#
# Files are named "Mera_Pleiades_2012-11-25_DEM_4m.tif" and "Mera_Pleiades_2018-10-28_DEM_4m.tif".

dem_2012 = gu.Raster("../data/04_mb_Mera/rasters/Mera_Pleiades_2012-11-25_DEM_4m.tif")
dem_2018 = gu.Raster("../data/04_mb_Mera/rasters/Mera_Pleiades_2018-10-28_DEM_4m.tif")

# ### Load the glacier outlines: Mera outlines from 2012 and 2018, RGI outlines
#
# Files are named: "Mera_outline_2012_realigned.shp", "Mera_outline_2018_realigned.shp" and "Glacier_inventory_around_Mera.shp".

outlines_2012 = gu.Vector("../data/04_mb_Mera/glacier_outlines/Mera_outline_2012_realigned.shp")
outlines_2018 = gu.Vector("../data/04_mb_Mera/glacier_outlines/Mera_outline_2018_realigned.shp")
outlines_rgi = gu.Vector("../data/04_mb_Mera/RGI_shapefiles/Glacier_inventory_around_Mera.shp")

# #### **Questions:** 
# - What is the spatial resolution of the DEMs? Is it the same for both?
# - What is the coordinate reference system of the DEMs and outlines? Is it the same for all?
# - Are the DEM extents the same?
#
# Use the code cell below to type your commands and find the answer to the questions.

print(dem_2012.info())

print(dem_2018.info())

print(dem_2012.crs.to_wkt())

print(outlines_2012.crs.to_wkt())

print(outlines_rgi.crs.to_wkt())

# ### Answers:
# => Both DEMs have the same spatial resolution: 4 m  
# => All files have the same CRS: EPSG 32645 = UTM zone 45N, except for RGI which is in lat/lon  
# => The DEM extents differ and hence the raster sizes !

# # Terrain attributes

# #### Calculate the hillshade of 2012 and 2018 DEMs

hillshade_2012 = xdem.terrain.hillshade(dem_2012)
hillshade_2018 = xdem.terrain.hillshade(dem_2018)

# #### Plot both hillshades and overlay the outlines of Mera glacier for the corresponding year

# +
fig = plt.figure(figsize=(10,6))

ax1 = plt.subplot(121)
hillshade_2012.show(cmap='gray', add_cbar=False)
outlines_2012.show(ax=ax1, facecolor='none', edgecolor='k')
ax1.set_title('Mera glacier and 2012 DEM')

ax2 = plt.subplot(122)
hillshade_2018.show(cmap='gray', add_cbar=False)
outlines_2018.show(ax=ax2, facecolor='none', edgecolor='k')
ax2.set_title('Mera glacier and 2012 DEM')

plt.tight_layout()
plt.show()
# -

# ### Comment what you see:
# - what are the white (transparent) areas on the figure?
# - which DEM has a larger extent?

# #### Calculate and plot slope and aspect for the 2018 DEM

slope = xdem.terrain.slope(dem_2018)
plt.figure(figsize=(10,10))
slope.show(cbar_title="Slope (degrees)")
plt.show()

# #### Calculate and plot aspect for the 2018 DEM

aspect = xdem.terrain.aspect(dem_2018)
plt.figure(figsize=(10,10))
aspect.show(cbar_title="Aspect (degrees)", cmap="twilight")
plt.show()

# # Calculating the difference between the two DEMs

ddem = dem_2018 - dem_2012

# #### An error occurs ! Why is that?

# # Handling rasters projection and georeferencing

# ### Reproject the two DEMs on the same grid
# Check afterwards that both DEM have same shape and georeferences.

dem_2012_proj = dem_2012.reproject(dst_ref=dem_2018)
print(dem_2012_proj.info())
print(dem_2018.info())

# ### Now calculate the elevation change again

ddem = dem_2018 - dem_2012_proj

# ## Plot the elevation change map

# #### Make a plot with the following features: 
# 1) Plot the elevation difference with a color scale of +50/-50 m
# 2) Plot Mera glacier outlines in black
# 3) Include a color scale label
# 4) Include a figure title

plt.figure(figsize=(10, 10))
ax = plt.subplot(111)
outlines_2012.show(ax=ax, facecolor='none', edgecolor='b', zorder=2)
outlines_2018.show(ax=ax, facecolor='none', edgecolor='r', zorder=3)
outlines_rgi.show(ax=ax, facecolor='none', edgecolor='k', zorder=4)
ddem.show(ax=ax, cmap='RdYlBu', vmin=-50, vmax=50, cbar_title='Elevation change 2012 - 2018 (m)', zorder=1)
ax.set_title('Thinning glaciers near Mera peak')
plt.tight_layout()
plt.show()

# ## Save the results in a file named "tmp_Pleiades_2012_2018_dh.tif"

ddem.save("../data/tmp_Pleiades_2012_2018_dh.tif")

# # Calculate zonal statistics
# Here we want to calculate the mean elevation change on and off glaciers.

# ## Rasterize the RGI glacier outlines on the same grid as ddem

glaciers_mask = outlines_rgi.create_mask(ddem)
plt.figure(figsize=(8,8))
glaciers_mask.show(add_cbar=False)
plt.show()

# ### Calculate mean dh over glaciers or stable terrain

# Over glaciers:

print(np.mean(ddem[glaciers_mask]))

# Over stable terrain

print(np.mean(ddem[~glaciers_mask]))

# ### Can you explain what is the issue?


