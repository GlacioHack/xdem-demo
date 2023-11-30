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

# # Manipulating raster/vector data with geoutils - Exercise

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

dem_2012 = # use the command to load a Raster to open the 2012 DEM
dem_2018 = # use the command to load a Raster to open the 2018 DEM

# ### Load the glacier outlines: Mera outlines from 2012 and 2018, RGI outlines
#
# Files are named: "Mera_outline_2012_realigned.shp", "Mera_outline_2018_realigned.shp" and "Glacier_inventory_around_Mera.shp".

outlines_2012 = # use the command to load a Vector to open the 2012 outlines
outlines_2018 = # use the command to load a Vector to open the 2018 outlines
outlines_rgi = # use the command to load a Vector to open the RGI outlines

# #### **Questions:** 
# - What is the spatial resolution of the DEMs? Is it the same for both?
# - What is the coordinate reference system of the DEMs and outlines? Is it the same for all?
# - Are the DEM extents the same?
#
# Use the code cell below to type your commands and find the answer to the questions.

print(dem_2012.info())

# +
# type in other commands here to find the answer to the questions
# -

# ### Answers:
#

# # Terrain attributes

# #### Calculate the hillshade of 2012 and 2018 DEMs

hillshade_2012 = # to be filled in
hillshade_2018 = # to be filled in

# #### Plot both hillshades and overlay the outlines of Mera glacier for the corresponding year

# +
fig = plt.figure(figsize=(10,6))

ax1 = plt.subplot(121)
# to be filled in

ax2 = plt.subplot(122)
# to be filled in

plt.tight_layout()
plt.show()
# -

# ### Comment what you see:
# - what are the white (transparent) areas on the figure?
# - which DEM has a larger extent?

# #### Calculate and plot slope and aspect for the 2018 DEM

slope = # to be filled in
plt.figure(figsize=(10,10))
# to be filled in
plt.show()

# #### Calculate and plot aspect for the 2018 DEM

aspect = # to be filled in
plt.figure(figsize=(10,10))
# to be filled in
plt.show()

# # Calculating the difference between the two DEMs

ddem = # to be filled in

# #### An error occurs ! Why is that?

# # Handling rasters projection and georeferencing

# ### Reproject the two DEMs on the same grid
# Check afterwards that both DEM have same shape and georeferences.

dem_2012_proj = # to be filled in
print(dem_2012_proj.info())
print(dem_2018.info())

# ### Now calculate the elevation change again

ddem = # to be filled in

# ## Plot the elevation change map

# #### Make a plot with the following features: 
# 1) Plot the elevation difference with a color scale of +50/-50 m
# 2) Plot Mera glacier outlines in black
# 3) Include a color scale label
# 4) Include a figure title

plt.figure(figsize=(10, 10))
ax = plt.subplot(111)
# to be filled in
plt.tight_layout()
plt.show()

# ## Save the results in a file named "Pleiades_2012_2018_dh.tif"

# +
# to be filled in
# -

# # Calculate zonal statistics
# Here we want to calculate the mean elevation change on and off glaciers.

# ## Rasterize the RGI glacier outlines on the same grid as ddem

glaciers_mask = # to be filled in
plt.figure(figsize=(8,8))
# to be filled in
plt.show()

# ### Calculate mean dh over glaciers or stable terrain

# Over glaciers:

# +
# to be filled in
# -

# Over stable terrain

# +
# to be filled in
# -

# ### Can you explain what is the issue?


