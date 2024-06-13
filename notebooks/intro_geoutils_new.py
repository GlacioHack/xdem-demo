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

# # Manipulating raster/vector data in Python

# In this tutorial, we will learn some basic operations on raster/vector data in Python such as:
# - read, write, reproject and display rasters and vectors
# - running arithmetic operations on rasters
# - rasterizing vectors, polygonizing rasters etc
#
# We will make use of the **geoutils** library. This is not a "standard" Python library but was developed by glaciologists (R. Hugonnet, E. Mannerfelt, A. Dehecq, A. Tedstone), to make it easier to work with geospatial data (raster, vector) in Python. They are based upon the more known **rasterio/GDAL** and **geopandas** libraries.
#
# The full documentation can be found at: \
# https://geoutils.readthedocs.io. \
#
# For more features, check the examples gallery on the documentation.

# ## Import the necessary modules

# +
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib widget

import geoutils as gu

plt.rcParams["image.interpolation"] = "none"
# -

# ## Sample data ###
# The code comes with some data for showcasing the functionality. If not done already, the data will be downloaded when needed (takes a few seconds). \
# Here we use two Landsat scenes over Mount Everest (one single-band, one RGB image), an ASTER DEM over Patagonia and glacier outlines from RGI. \
# All available samples can be shown with:

print(gu.examples.available)

# The files are already on disk, we only need to find their location, e.g.

print(gu.examples.get_path("everest_landsat_b4"))

# ### <span style='color:red '> **TO DO:** </span> Find these files on your computer and open them in QGIS.
#
#

# ## Read the three rasters and glacier outlines ###

# Raster files (e.g. GeoTiff) can be loaded in one line with `gu.Raster(path_to_file)`.

landsat_b4 = gu.Raster(gu.examples.get_path("everest_landsat_b4"))
landsat_rgb = gu.Raster(gu.examples.get_path("everest_landsat_rgb"))
aster_dem = gu.Raster(gu.examples.get_path("exploradores_aster_dem"))


# ### <span style='color:red '> **TO DO:** </span> Run the same commands, but copy/paste the file path. Uncomment the lines below (remove the # symbol) and replace ... with the file path.

# +
# dem_2009 = gu.Raster(...)
# dem_1990 = gu.Raster(...)
# -

# Vector files (e.g. ESRI shapefiles) can be loaded in one line with `gu.Vector(path_to_file)`.

outlines = gu.Vector(gu.examples.get_path("everest_rgi_outlines"))

# ## Quickly visualize a raster
# Since a Raster object comes with all atributes, it can be quickly plotted with its georeferencing information.

plt.figure()
landsat_b4.plot()
plt.show()

# The default colorbar is viridis, one can set it grays to display in B&W

plt.figure()
landsat_b4.plot(cmap="gray", add_cbar=None)
plt.show()

# For 3-band images, `Raster.plot` displays as RGB by default

plt.figure()
landsat_rgb.plot()
plt.show()

# ## Quickly visualize vector data

plt.figure()
outlines.plot()
plt.show()

# ## Superimpose layers or play with axes
#
# One can pass a matplotlib Axes instance as argument to the plot function.


# +
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

aster_dem.plot(ax=axes[0], title="ASTER DEM")

landsat_rgb.plot(ax=axes[1], title="Landsat RGB + outlines")
outlines.reproject(landsat_rgb).plot(ax=axes[1], ec="k", fc="none")

plt.tight_layout()
plt.show()
# -

# ## Attributes of the Raster and Vector classes
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

print(f"EPSG code: {landsat_rgb.crs.to_epsg()}\n")
print(f"Printed as WKT string: \n{landsat_rgb.crs.to_wkt()}")

# **Raster width, height, resolution and bounds**

landsat_rgb.width

landsat_rgb.height

landsat_rgb.res

landsat_rgb.bounds

# The resolution and position of the upper left pixel are traditionally stored in a so-called [transform](https://rasterio.readthedocs.io/en/latest/topics/transforms.html):

landsat_rgb.transform

# **These information, and more, can all be obtained at once with the command**

landsat_rgb.info()

# Along with these metadata, the `geoutils.Raster` object contains the data, stored as a numpy masked array in the `self.data` attribute:

landsat_rgb.data

# ### The `gu.Vector` instances are based upon **[GeoPandas](https://geopandas.org)**. 
#
# You can see a full list of features at [this link](https://geoutils.readthedocs.io/en/stable/vector_class.html). \
# The underlying GeoPandas' DataFrame can be accessed via:

outlines.ds


# #### Additionally, it has similar attributes to those of the `Raster` class, such as:

outlines.crs

outlines.bounds

# # Raster operations

# ### Note on nodata
#
# Nodata values are read from file and nodata values are masked using np.ma.masked_arrays.
# The nodata value is stored in `self.nodata` and the mask in `self.data.mask`.

print(aster_dem.nodata)

plt.figure()
plt.imshow(aster_dem.data.mask)
plt.show()

# The best way to modify the mask is through `self.set_mask`.
# All pixels where mask is set to True or > 0 will be masked (in addition to previously masked pixels).

# Let's mask all pixels with elevation above 2000 m

aster_dem_masked = aster_dem.copy()
aster_dem_masked.set_mask(aster_dem > 2000)

fig, axes = plt.subplots(1, 2, figsize=(11, 5))
aster_dem.plot(ax=axes[0], title="Before masking")
aster_dem_masked.plot(ax=axes[1], title="After masking")
plt.tight_layout()
plt.show()

# ### Reproject the two DEMs on the same grid
# NB: Here the function returns a warning because the two DEMs are already on the same grid, so nothing is actually done.

landsat_b4 = landsat_b4.reproject(ref=landsat_rgb)

# Check that both DEMs indeed have the same grid

landsat_b4.info()
landsat_rgb.info()

# ### Reproject to a given resolution, bounds, or CRS

# #### Change pixel resolution

rast_test = landsat_b4.reproject(res=60)
rast_test.info()

# #### **Question:** What is the new raster height?

# #### Change extent/bounds

rast_test = landsat_b4.reproject(bounds={"left":483430.0, "bottom":3093260.0, "right":498190.0, "top":3102710.0})
rast_test.info()

# #### **Question:** What is the new raster height?

# #### Change Coordinate Reference System (CRS) i.e. projection

rast_test = landsat_b4.reproject(crs='epsg:4326')
rast_test.info()

# #### **Question:** What is the new pixel resolution? What are the units?

# ### Reproject the outlines in the same coordinate system as DEMs

outlines_proj = outlines.reproject(crs=landsat_b4.crs)

# ### Cropping
# Cropping differ from reprojecting in the sense that it does not perform any resampling, i.e. if the bounds do not exactly match the raster grid, they will be rounded

rast_test = landsat_b4.crop(crop_geom=(483530.0, 3093260.0, 498190.0, 3102730.0))
rast_test.info()

# #### **Question:** By how much does the left bound differ from the requested value?

# ### Arithmetic operations
#
# A Raster can be applied any pythonic arithmetic operation (+, -, /, //, *, **, %) with another Raster, ndarray or number. It will output one or two Rasters. NumPy coercion rules apply for dtype.

calc_rast = (landsat_b4 + 1)/2

plt.figure()
calc_rast.plot()
plt.show()

# #### **Question:** What happens with white areas?

# A solution is to change the data type before running the operation.

calc_rast = (landsat_b4.astype(np.float32) + 1)/2

plt.figure()
calc_rast.plot()
plt.show()

# A Raster can also be applied any pythonic logical comparison operation (==, !=, >=, >, <=, <) with another Raster, ndarray or number. It will cast to a Mask.

calc_mask = landsat_b4 < 100

plt.figure()
calc_mask.plot()
plt.show()


# ### Array interface
#
# A Raster can be applied any NumPy universal functions and most mathematical, logical or masked-array functions with another Raster, ndarray or number.

calc_rast = np.sqrt(landsat_b4)

np.max(calc_rast)


# ### Polygonize

# Get a Vector instance containing polygons of all pixels with landsat B4 values > 200 (~snow)

snow = (landsat_b4 > 200).polygonize()

fig, ax = plt.subplots(1, 1)
landsat_b4.plot(cmap="gray", ax=ax, add_cbar=False)
snow.plot(ax=ax, alpha=0.5)
plt.show()

# ### Proximity

snow_proximity = snow.proximity()

plt.figure()
snow_proximity.plot()
plt.show()

# ### Saving and exporting

# Saving a raster/vector to file

snow_proximity.save("my_raster.tif")
snow.save("my_vector.gpkg")

# Exporting raster to numpy array with nans

calc_rast.get_nanarray()

# Exporting a raster to xarray or rasterio

aster_dem.to_xarray()

aster_dem.to_rio_dataset()

# Exporting a raster to vector

landsat_pc = landsat_rgb.reproject(res=200).to_pointcloud()
landsat_pc.plot(column="b1", cmap="gray")
plt.show()

# # Vector operations
#
# ### All `Raster` methods that have an equivalent for vectors also exist for `Vector` instances, this includes: reproject (can only update CRS), crop and proximity. Other methods are illustrated below.

# ### Rasterizing 
#
# -> returns a Mask

glacier_mask = outlines.create_mask(landsat_b4)
glacier_mask

plt.figure()
glacier_mask.plot(add_cbar=False)
plt.show()

# ### Calculating zonal statistics
#
# We can use the mask to calculate e.g. the mean value on and off glacier

# Over glaciers:

print(np.mean(landsat_b4[glacier_mask]))

# Over stable terrain

print(np.mean(landsat_b4[~glacier_mask]))

# ### Buffer
#
# We grow the outlines by 500 m, which returns a new Vector instance.

# +
outlines_metric = outlines.reproject(crs=landsat_b4.crs)
buffered_outlines = outlines_metric.buffer(500)

plt.figure()
buffered_outlines.plot()
plt.show()
