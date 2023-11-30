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

# # Manipulating raster/vector data with geoutils

# **geoutils** was developed by glaciologists (R. Hugonnet, E. Mannerfelt, A. Dehecq), to make it easier to handle raster and vector data in Python.\
# **xDEM** uses all functionnalities of geoutils and add additional tools specific to working with DEMs. \
# The full documentation can be found at: \
# https://geoutils.readthedocs.io. \
# https://xdem.readthedocs.io.
#
# This tutorial highlights some basic functionalities of geoutils/xdem:
# - read, write, reproject and display rasters and vectors
# - calculate terrain attributes from DEMs (slope, aspect etc.)
# - calculate difference of DEMs and basic statistics.
#
# For more examples, check the examples gallery on both geoutils and xdem's documentation.

# ## Import the necessary modules

# +
import matplotlib.pyplot as plt
import numpy as np

import geoutils as gu
import xdem
# -

# ## Download the sample data set (if not done already) - Should take a few seconds ###
# The code comes with some data for showcasing the functionality, here two DEMs over Longyearbyen, Norway and glacier outlines over Svalbard. \
# The files are already on disk, we only nned to find their location.

xdem.examples.download_longyearbyen_examples(overwrite=False)
print(xdem.examples.get_path("longyearbyen_ref_dem"))
print(xdem.examples.get_path("longyearbyen_tba_dem"))
print(xdem.examples.get_path("longyearbyen_glacier_outlines"))

# ### **Side work:** Find these files on your computer and open them in QGIS.

# ## Read the two DEMs and glacier outlines for the region ###

# Raster files (e.g. GeoTiff) can be loaded in one line with `gu.Raster(path_to_file)`.

dem_2009 = gu.Raster(xdem.examples.get_path("longyearbyen_ref_dem"))
dem_1990 = gu.Raster(xdem.examples.get_path("longyearbyen_tba_dem"))

# ### **Side work:** Run the same commands, but copy/paste the file path. Uncomment the lines below (remove the # symbol) and replace ... with the file path.

# +
# dem_2009 = gu.Raster(...)
# dem_1990 = gu.Raster(...)
# -

# Vector files (e.g. ESRI shapefiles) can be loaded in one line with `gu.Vector(path_to_file)`.

outlines_1990 = gu.Vector(xdem.examples.get_path("longyearbyen_glacier_outlines"))

# ## Quickly visualize a raster
# Since a Raster object comes with all atributes, it can be quickly plotted with its georeferencing information.

dem_2009.show()

# It is easier to visualize as a hillshade

dem_2009_hs = xdem.terrain.hillshade(dem_2009)
dem_2009_hs.show(cmap='gray')

# ## Quickly visualize vector data

outlines_1990.show()

# ## Notes on the Raster and Vector classes
#
# ### The `geoutils.Raster` class is based upon **[rasterio](https://rasterio.readthedocs.io)** and inherits most of its metadata. 
# These objects contain the raster metadata, with the same convention as rasterio. 
# To georeference an object, one needs to know the coordinate reference system (called `crs` in rasterio/geopandas) and any 3 of these four information:
# - the raster width and height
# - the pixel resolution
# - the extent of the raster (called `bounds` in rasterio/geopandas)
# - the position of one pixel, traditionally, the upper-left \
# These variables are inter-dependent, e.g. if one knows the raster's extent and width and height, the pixel resolution is fixed.
# All these variables are stored in the `xdem.DEM` instance with the following attributes:

print(f"As PROJ4 string: {dem_2009.crs.to_proj4()}\n")
print(f"As WKT string: \n{dem_2009.crs.to_wkt()}")

dem_2009.width

dem_2009.height

dem_2009.res

dem_2009.bounds

# The resolution and position of the upper left pixel are traditionally stored in a so-called [transform](https://rasterio.readthedocs.io/en/latest/topics/transforms.html):

dem_2009.transform

# These information, and more, can all be obtained **at once** with the command

print(dem_2009.info())

# Along with these metadata, the `xdem.DEM` object contains the data, stored as a numpy masked array in the `self.data` attribute:

dem_2009.data

# ### `gu.Vector` instances are based upon geopandas. 
# The class contains several useful methods (`self.create_mask` is showcased below), and the GeoDataFrame can accessed via:

outlines_1990.ds

# # Terrain attributes
#
# **xDEM** enables calculating many terrain attributes: slope, aspect, hillshade, curbature etc. The full list is available here: https://xdem.readthedocs.io/en/stable/terrain.html
#
# Here we will demonstrate a few of them.

# ### Slope

slope = xdem.terrain.slope(dem_1990)
slope.show(cbar_title="Slope (degrees)")

# ### Aspect

aspect = xdem.terrain.aspect(dem_1990)
aspect.show(cbar_title="Aspect (degrees)", cmap="twilight")

# #### **Side work:** 
# 1) Plot the terrain rugosity
# 2) Calculate the maximum slope (using `np.max`)

# +
# rugosity = xdem.terrain.???(dem_1990)
# rugosity.show(cbar_title="...", cmap="...")

# +
# print("Maximum slope is:", np.max(...))
# -

# # Rasters operations

# ## Reproject the two DEMs on the same grid
# Here the function returns a warning because the two DEMs are already on the same grid, so nothing is done.

dem_1990 = dem_1990.reproject(dst_ref=dem_2009)
print(dem_1990.info())
print(dem_2009.info())

# ## Reproject to a given resolution, bounds, or CRS

# #### Change pixel resolution

dem_test = dem_1990.reproject(dst_res=60)
print(dem_test.info())

# #### **Question:** What is the new raster height?

# #### Change extent/bounds

dem_test = dem_1990.reproject(dst_bounds={"left":502810, "top":8674000, "right":529430, "bottom": 8654290})
print(dem_test.info())

# #### **Question:** What is the new raster height?

# #### Change Coordinate Reference System (CRS) i.e. projection

dem_test = dem_1990.reproject(dst_crs='epsg:4326')
print(dem_test.info())

# #### **Question:** What is the new pixel resolution? What are the units?

# ## Reproject the outlines in the same coordinate system as DEMs

outlines_proj = outlines_1990.reproject(dst_crs=dem_2009.crs)

# # Calculating the difference between two DEMs

ddem = dem_2009 - dem_1990

# ## Plot the elevation change map
# #### Note: 
# - `ax` is used here to share the same subplot between the raster and outlines (the default is to create a new figure)
# - ddem is plotted last, to preserve the extent, as glacier outlines cover all of Svalbard.
# - zorder is used to plot in the right sequence (outlines on top)

# +
vmax = max(abs(np.max(ddem.data)), abs(np.min(ddem.data)))

ax = plt.subplot(111)
outlines_proj.show(ax=ax, facecolor='none', edgecolor='k', zorder=2)
ddem.show(ax=ax, cmap='RdYlBu', vmin=-vmax, vmax=vmax, cbar_title='Elevation change 2009 - 1990 (m)', zorder=1)
ax.set_title('Thinning glaciers near Longyearbyen')
plt.tight_layout()
plt.show()
# -

# #### **Side work:** Make changes to the plot: 
# 1) Plot the glacier outlines in red.
# 2) Use a different min/max value for colorscale.
# 3) Use a different colormap (choose among the list [here](https://matplotlib.org/stable/gallery/color/colormap_reference.html)).
# 4) Change the figure title.
#
# (Uncomment the lines below and replace "..." with a different value)

# +
# fig = plt.figure()
# ax = plt.subplot(111)
# outlines_proj.show(ax=ax, facecolor='none', edgecolor='...', zorder=2)
# ddem.show(ax=ax, cmap='coolwarm', vmin=..., vmax=..., cbar_title='Elevation change 2009 - 1990 (m)', zorder=1)
# ax.set_title('...')
# plt.tight_layout()
# plt.show()
# -

# ## Saving the results

ddem.save("temp_ddem.tif")

# # Calculate zonal statistics
# Here we want to calculate the mean elevation change on and off glaciers.

# ## Rasterize the glacier outlines on the same grid as ddem
# `glacier_mask` is `True` on glaciers, `False` elsewhere.

glacier_mask = outlines_1990.create_mask(ddem)
glacier_mask.show(add_cbar=False)

# ### Calculate mean dh over glaciers or stable terrain

# Over glaciers:

print(np.mean(ddem[glacier_mask]))

# Over stable terrain

print(np.mean(ddem[~glacier_mask]))

# ### Something is wrong, mean dh over stable terrain should be ~0 => we need to coregister the DEMs. This is some proper work for our next example !


