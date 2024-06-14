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

# # Vegetation in Paris

# <div class="alert alert-info" style="font-size:110%">
#  <p>
#  <h3>Activity</h3>
#  </p>
#  
#  Locate vegetation cover in Paris, thanks to Sentinel-2 satellite images.
#  
# </div>

# +
import numpy as np
import geopandas as gpd

import matplotlib
import matplotlib.pyplot as plt

from rasterio.plot import show_hist
import geoutils as gu
# -

# ## Discover the structure of an image file
#
# An image is a set of matrices with associated metadata
#
# ![image.png](img/image.png)

# <div class="alert alert-info" style="font-size:110%">
#  
#  <ol>
#     <li>Open the satellite image of Paris (file : T31UDQ_20190705T105031_TCI_10m.tif), in true colors, taken in July 2019.</li>
#     <li>Display its metadata :</li>
#     <ul>
#         <li>Filename</li>
#         <li>Projection</li>
#         <li>Spatial extent</li>
#         <li>Number of bands</li>
#      </ul>
#  </ol>
#  
# All data are available in folder <b>data/Paris/</b>
# </div>

rgb_rst = gu.Raster("../../data/PARIS/images/T31UDQ_20190705T105031_TCI_10m.tif")
print(rgb_rst.info())

# <div class="alert alert-info" style="font-size:110%">
#
#   <ul>
#     <li>Read and display the number of bands</li>
#     <li>Display band types and shape</li>
#   </ul>
#
# </div>

rgb_rst.count

rgb_rst.dtype, rgb_rst.shape

rgb_rst.data

# <div class="alert alert-info" style="font-size:110%">
#
#   <ul>
#     <li>Display the image in true color</li>
#   </ul>
#
# </div>

rgb_rst.plot()


# ## Making a false color composition
#
# False color compositions with infrared are common to highlight elements in the infrared, in particular vegetation.
#
# ![Infrarouge.jpg](img/Infrarouge.jpg)  
# *Infrared image, vegetation particularly stands out.*

# #### Exploring the NIR band
# <div class="alert alert-info" style="font-size:110%">
#
# - Open the near infrared (NIR) band, acquired at the same time (file : T31UDQ_20190705T105031_B08_10m.tif).
# - Compare the values of the pixels from the true color image and the one from the infrared band. You can use the function `show_hist` to plot histograms.
#
# </div>

nir_rst = gu.Raster("../../data/PARIS/images/T31UDQ_20190705T105031_B08_10m.tif")
nir_rst.plot()

fig, axes = plt.subplots(ncols=2, figsize=(12, 6))
show_hist(rgb_rst.data, bins=50, ax=axes[0], label=["Red", "Green", "Blue"], alpha=0.3, histtype="stepfilled")
show_hist(nir_rst.data, bins=50, ax=axes[1], label="PIR", alpha=0.3, histtype="stepfilled")


# #### Normalizing the NIR band
#
# <div class="alert alert-info" style="font-size:110%">
#
# - Create a function `normalize` to normalize the data in the NIR band between 0-255 and apply it to the opened band
# </div>

def normalize(rst):
    band_min, band_max = rst.data.min(), rst.data.max()
    normalized = ((rst - band_min)/(band_max - band_min)) * 255
    return normalized.astype("int32")


nnir_rst = normalize(nir_rst)
show_hist(nnir_rst.data, bins=50, label="PIR", alpha=0.3, histtype="stepfilled")

# #### Creating the false color composite
#
# <div class="alert alert-info" style="font-size:110%">
#
#   <ul>
#     <li>Display a false color composition (NIR, Red, Green), using the normalized NIR band.</li>
#   </ul>
#
# Note: you can use function `Raster.split_bands` to first create single layer Rasters, then `Raster.from_array` to generate the composite.
# </div>

red, green, blue = rgb_rst.split_bands()

nnir_rv = gu.Raster.from_array(
    data=np.stack((nnir_rst.data, red.data, green.data)), 
    crs=rgb_rst.crs, 
    transform=rgb_rst.transform,
    nodata=rgb_rst.nodata
)
nnir_rv.plot()

# ## Index calculation - NDVI
# Calculating an index from raster is a common practice for highlighting certain features. The best-known of these indices is NDVI: *Normalized Difference Vegetation Index*, which highlights vegetation zones.
#
# It is based on pixel values in the R and NIR bands:
#
# $NDVI = \frac{PIR-R}{PIR+R}$
#
# This index ranges from -1 to 1, and its values are generally interpreted according to the following scale:
# * Negative values correspond to surfaces composed essentially of water (snow, rivers, etc.).
# * Pixels with values around zero describe bare or artificial ground.
# * Pixels with values between 0.2 and 0.3 (*approximately*) represent areas of low vegetation (shrubs, meadows, etc.).
# * Higher values ((*approximately*) 0.6-0.8) indicate highly vegetated areas (forests).

# <div class="alert alert-info" style="font-size:110%">
#
#   <ul>
#     <li>Calculate the NDVI on 5 July 2019.</li>
#     <li>Save the results in an image file.</li>
#     <li>Display the NDVI, using the “RdYlGn” color palette.</li>
#   </ul>
#
# Remember that an image is a set of matrix(es) and metadata. Matrices cannot be saved on their own without metadata.
#
# </div>

ndvi = (nnir_rst - red) / (nnir_rst + red)
ndvi.data[(nnir_rst + red).data == 0] = 0
ndvi.data.mask[(nnir_rst + red).data == 0] = False
ndvi.plot(cmap="RdYlGn")

# ### Cropping the image

# <div class="alert alert-info" style="font-size:110%">
#
#   <ul>
#     <li>Crop the NDVI according to Paris boundaries.</li>
#     <li>Display the result with the same color palette as above.</li>
#   </ul>
#
# A file with Paris boundaries is located in **data/PARIS/limites_paris.geojson**.
#
# </div>

paris_shape = gu.Vector("../../data/PARIS/limites_paris.geojson")
paris_shape.plot()

ndvi_crop = ndvi.crop(paris_shape)
ndvi_crop.plot(cmap="RdYlGn")

paris_mask = paris_shape.create_mask(ndvi)
ndvi_clip = ndvi_crop.copy()
ndvi_clip.set_mask(~paris_mask)
ndvi_clip.plot(cmap="RdYlGn")

# ### Extracting vegetation zones

# <div class="alert alert-info" style="font-size:110%">
#
#   <ul>
#     <li>Apply a threshold to the NDVI to extract vegetated areas (two cases below).</li>
#     <li>Plot the results.</li>  
#     <li>Save the extracted areas as vector files.</li>
#   </ul>
#
# </div>

threshold_veg = 0.2
threshold_water = -0.9

# #### Case with 2 categories (True = vegetation, False = no vegetation)

is_vegetation = (ndvi > threshold_veg)
print(type(is_vegetation))
is_vegetation.plot(add_cbar=False)

# ### Case with 4 categories
# 0 = water \
# 1 = bare ground \
# 2 = vegetation \
# 10 = out of bounds pixels

paris_classif = np.ones(ndvi_clip.shape, dtype="uint8")
paris_classif[ndvi_clip.data.mask] = 10
paris_classif[(ndvi_clip.data > threshold_veg) & (paris_classif < 10)] = 2
paris_classif[(ndvi_clip.data <= threshold_water)] = 0
classif_rst = nir_rst.copy(new_array=paris_classif)

classif_rst.plot(interpolation="none", vmax=3)

classif_vect = classif_rst.polygonize(target_values=[0, 1, 2])
type(classif_vect)

classif_vect.ds

# Rename column "raster_value" by "category"
classif_vect.ds = classif_vect.ds.rename(columns={"raster_value": "category"})
classif_vect.ds

classif_vect.plot(column="category")

classif_vect.save("my_vegetation.gpkg")

m = classif_vect.ds.explore("category", name="vegetation", legend=True)
m

# ### BONUS - Use folium to plot the results

# <div class="alert alert-info" style="font-size:110%">
#
# - convert the NDVI Raster into a numpy array with Nan, projected into web mercator (ESRI:53004)
# - create the interactive map with `geopandas.ds.explore` as above (easier)
# - add a Folium image layer with `folium.raster_layers.ImageOverlay`. You need to set a colormap.
# - display the colorbar using information from [this link](https://leafmap.org/notebooks/62_folium_colorbar/)
# - save the output as an html file for later use
# - Play with the different options to customize your plot
#
# </div>

ndvi_reproj = ndvi_clip.reproject(crs="ESRI:53004")
bounds = ndvi_reproj.get_bounds_projected("EPSG:4326")
ndvi_nan = ndvi_reproj.get_nanarray()

plt.imshow(ndvi_nan)
plt.colorbar()

cmap = plt.get_cmap("viridis")
cmap(0)

# +
import folium

m = classif_vect.ds.explore("category", name="vegetation", legend=True)
img = folium.raster_layers.ImageOverlay(
        name="NDVI",
        image=ndvi_nan,
        bounds=[[bounds.bottom, bounds.left], [bounds.top, bounds.right]],
        colormap=lambda x: cmap(x),
        opacity=0.6,
    )
m.add_child(img)
folium.LayerControl().add_to(m)

m
# -
m.save("test.html")



