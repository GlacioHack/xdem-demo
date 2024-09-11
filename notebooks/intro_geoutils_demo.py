# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Manipulating raster/vector data in Python

# <div class="alert alert-info" style="font-size:110%">
# <h3>Objectives</h3>
#     
# In this tutorial, we will learn some basic operations on raster/vector data in Python, including DEMs, such as:
# - read, write, reproject and display rasters and vectors
# - calculate terrain attributes from DEMs (slope, aspect etc.)
# - calculate difference of DEMs and basic statistics.
#
# We will make use of the **geoutils** and **xdem** libraries. They are not "standard" Python libraries but were developed by glaciologists (R. Hugonnet, E. Mannerfelt, A. Dehecq), to make it easier to work with such data in Python. They are based upon the more known **rasterio/GDAL** and **geopandas** libraries.
# - **geoutils** provides functionalities for working with raster/vector data.
# - **xDEM** uses all functionnalities of geoutils and add additional tools specific to working with DEMs.
#
# The full documentation can be found at: \
# https://geoutils.readthedocs.io. \
# https://xdem.readthedocs.io.
#
# For more features, check the examples gallery on both geoutils and xdem's documentation.
#
#  </div>

# ## Import the necessary modules

# +
import matplotlib.pyplot as plt
import numpy as np

import geoutils as gu
import xdem

# for interactive plots only, in case of issues with plots, comment this line and restart the kernel
# %matplotlib widget

# To prevent interpolation on plots
plt.rcParams["image.interpolation"] = "none"
# -

# # 1) Quick-start: loading and visualizing georeferenced data

# ## Download the sample data set (if not done already) - Should take a few seconds ###
# The code comes with some data for showcasing the functionality, here two DEMs over Longyearbyen, Norway and glacier outlines over Svalbard. \
# The files are already on disk, we only nned to find their location.

xdem.examples.download_longyearbyen_examples(overwrite=False)
print(xdem.examples.get_path("longyearbyen_ref_dem"))
print(xdem.examples.get_path("longyearbyen_tba_dem"))
print(xdem.examples.get_path("longyearbyen_glacier_outlines"))

# ### <span style='color:red '> **TO DO:** </span> Find these files on your computer.
#
#

# ## Read the two DEMs and glacier outlines for the region ###

# Raster files (e.g. GeoTiff) can be loaded in one line with `gu.Raster(path_to_file)`. The command load the metadata, but data is loaded only if/when needed.

dem_2009 = gu.Raster(xdem.examples.get_path("longyearbyen_ref_dem"))
dem_1990 = gu.Raster(xdem.examples.get_path("longyearbyen_tba_dem"))

# ### <span style='color:red '> **TO DO:** </span> Run the same commands, but copy/paste the file path. Uncomment the lines below (remove the # symbol) and replace ... with the file path.

# +
# dem_2009 = gu.Raster(...)
# dem_1990 = gu.Raster(...)
# -

# Vector files (e.g. ESRI shapefiles) can be loaded in one line with `gu.Vector(path_to_file)`.

outlines_1990 = gu.Vector(xdem.examples.get_path("longyearbyen_glacier_outlines"))

# ## Quickly visualize a raster
# Since a Raster object comes with all atributes, it can be quickly plotted with its georeferencing information.

plt.figure()
dem_2009.plot()
plt.show()

# It is easier to visualize as a hillshade

# +
dem_2009_hs = xdem.terrain.hillshade(dem_2009)

plt.figure()
dem_2009_hs.plot(cmap='gray')
plt.show()
# -

# ## Quickly visualize vector data

fig, ax = plt.subplots(figsize=(8, 10))
outlines_1990.plot(ax=ax)
plt.tight_layout()
plt.show()

# # 2) More details about the Raster class
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

# **Coordinate reference system (CRS)**

print(f"EPSG code: {dem_2009.crs.to_epsg()}\n")
print(f"Printed as WKT string: \n{dem_2009.crs.to_wkt()}")

# **Raster width, height, resolution and bounds**

dem_2009.width

dem_2009.height

dem_2009.res

dem_2009.bounds

# The resolution and position of the upper left pixel are traditionally stored in a so-called [transform](https://rasterio.readthedocs.io/en/latest/topics/transforms.html):

dem_2009.transform

# **These information, and more, can all be obtained at once with the command**

dem_2009.info()

# Along with these metadata, the `geoutils.Raster` object contains the data, stored as a numpy masked array in the `self.data` attribute:

dem_2009.data

# ### Note on nodata
#
# Nodata values are read from file and nodata values are masked using np.ma.masked_arrays.
# The nodata value is stored in `self.nodata` and the mask in `self.data.mask`.

print(dem_2009.nodata)

# The best way to modify the mask is through `self.set_mask`.
# All pixels where mask is set to True or > 0 will be masked (in addition to previously masked pixels).

# Let's mask all pixels with elevation above 500 m

dem_2009_masked = dem_2009.copy()
dem_2009_masked.set_mask(dem_2009 > 500)

fig, axes = plt.subplots(1, 2, figsize=(11, 5))
dem_2009.plot(ax=axes[0], title="Before masking")
dem_2009_masked.plot(ax=axes[1], title="After masking")
plt.tight_layout()
plt.show()

# #### Plot the mask

plt.figure()
plt.imshow(dem_2009_masked.data.mask)
plt.show()

# ## Raster reprojection

# ### Reproject the two DEMs on the same grid
# NB: Here the function returns a warning because the two DEMs are already on the same grid, so nothing is actually done. In most situation, your two input DEMs will be on different grids.

dem_1990 = dem_1990.reproject(ref=dem_2009)

# Check that both DEMs indeed have the same grid

dem_1990.info()
dem_2009.info()

# ### Reproject to a given resolution, bounds, or CRS

# #### Change pixel resolution

dem_test = dem_1990.reproject(res=60)
dem_test.info()

# #### **Question:** What is the new raster height?

# #### Change extent/bounds

dem_test = dem_1990.reproject(bounds={"left":502810, "top":8674000, "right":529430, "bottom": 8654290})
dem_test.info()

# #### **Question:** What is the new raster height?

# #### Change Coordinate Reference System (CRS) i.e. projection

dem_test = dem_1990.reproject(crs='epsg:4326')
dem_test.info()

# #### **Question:** What is the new pixel resolution? What are the units?


# ## Arithmetic operations
#
# A Raster can be applied any pythonic arithmetic operation (+, -, /, //, *, **, %) with another Raster, ndarray or number. It will output one or two Rasters. NumPy coercion rules apply for data type.

calc_rast = (dem_2009 + 300) / 2

plt.figure()
calc_rast.plot()
plt.show()

# ## Numpy array interface
#
# A Raster can be applied any NumPy universal functions and most mathematical, logical or masked-array functions with another Raster, ndarray or number.
# No data values are ignored during calculations

print(np.max(dem_2009))

print(np.mean(dem_2009))


# ## More functionalities
# The Raster class of geoutils has many more features such as: crop, polygonize, interpolate, proximity or export to other formats. For more details, check this page: https://geoutils.readthedocs.io/en/stable/raster_class.html

# # 3) More details about the Vector class
#
# ### The `gu.Vector` instances are based upon **[GeoPandas](https://geopandas.org)**. 
#
# The underlying GeoPandas' DataFrame can be accessed via:

outlines_1990.ds

# #### Additionally, it has similar attributes to those of the `Raster` class, such as:

outlines_1990.crs

outlines_1990.bounds

# ## Reproject the outlines in the same coordinate system as our DEM

outlines_proj = outlines_1990.reproject(crs=dem_2009.crs)

# ## More functionalities
# More features are available such as reproject, crop, rasterize etc. \
# A full list can be found at this link: https://geoutils.readthedocs.io/en/stable/vector_class.html.
#

# # 4) Introduction to working with DEMs

# The **xDEM** modules adds functionalities on top of the Raster class that are specific to Digital Elevation Models or DEMs.
# More examples will be provided in the next tutorial, but here we showcase a few basic functionalities.

# ## Terrain attributes
#
#  One feature of **xDEM** is the calculation of terrain attributes: slope, aspect, hillshade, curvature etc. The full list is available here: https://xdem.readthedocs.io/en/stable/terrain.html.
#
# Here we will demonstrate a few of them.

# ### Slope

slope = xdem.terrain.slope(dem_1990)
fig, ax = plt.subplots(figsize=(8, 6))
slope.plot(cbar_title="Slope (degrees)")
plt.show()

# ### Aspect

aspect = xdem.terrain.aspect(dem_1990)
fig, ax = plt.subplots(figsize=(8, 6))
aspect.plot(cbar_title="Aspect (degrees)", cmap="twilight")
plt.show()

# #### <span style='color:red '> **TO DO:** </span> 
# 1) Plot the terrain rugosity
# 2) Calculate the maximum slope (using `np.max`)  
# Uncomment (remove the #) below.

# +
# rugosity = xdem.terrain.???(dem_1990)
# rugosity.plot(cbar_title="...", cmap="...")

# +
# print("Maximum slope is:", np.max(...))
# -



# ## Calculating and plotting the difference between two DEMs

ddem = dem_2009 - dem_1990

# ### Plot the elevation change map
# #### Note: 
# - `ax` is used here to share the same subplot between the raster and outlines (the default is to create a new figure)
# - ddem is plotted last, to preserve the extent, as glacier outlines cover all of Svalbard.
# - zorder is used to plot in the right sequence (outlines on top)

# +
vmax = max(abs(np.max(ddem.data)), abs(np.min(ddem.data)))

fig, ax = plt.subplots(figsize=(10, 8))
ax = plt.subplot(111)
outlines_proj.plot(ax=ax, facecolor='none', edgecolor='k', zorder=2)
ddem.plot(ax=ax, cmap='coolwarm_r', vmin=-vmax, vmax=vmax, cbar_title='Elevation change 2009 - 1990 (m)', zorder=1)
ax.set_title('Thinning glaciers near Longyearbyen')
plt.tight_layout()
plt.show()
# -

# #### <span style='color:red '> **TO DO:** </span> Make changes to the plot: 
# 1) Plot the glacier outlines in red.
# 2) Use a different min/max value for colorscale.
# 3) Use a different colormap (choose among the list [here](https://matplotlib.org/stable/gallery/color/colormap_reference.html)).
# 4) Change the figure title.
#
# (Uncomment the lines below and replace "..." with a different value)

# +
# fig, ax = plt.subplots(figsize=(10, 8))
# ax = plt.subplot(111)
# outlines_proj.plot(ax=ax, facecolor='none', edgecolor='...', zorder=2)
# ddem.plot(ax=ax, cmap='coolwarm', vmin=..., vmax=..., cbar_title='Elevation change 2009 - 1990 (m)', zorder=1)
# ax.set_title('...')
# plt.tight_layout()
# plt.show()
# -

# ### Saving the results

ddem.save("temp_ddem.tif")

# ## Calculate zonal statistics
# Here we want to calculate the mean elevation change on and off glaciers.

# ### Rasterize the glacier outlines on the same grid as ddem
# `glacier_mask` is `True` on glaciers, `False` elsewhere.

# +
glacier_mask = outlines_1990.create_mask(ddem)

plt.figure()
glacier_mask.plot(add_cbar=False)
plt.show()
# -

# ### Calculate mean dh over glaciers or stable terrain

# Over glaciers:

print(np.mean(ddem[glacier_mask]))

# Over stable terrain

print(np.mean(ddem[~glacier_mask]))

# ### Something is incorrect... mean dh over stable terrain should be ~0 => we need to coregister the DEMs. This is some proper work for our next example !
