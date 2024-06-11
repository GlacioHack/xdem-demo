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



# <div class="alert alert-info" style="font-size:110%">
#
#   <ul>
#     <li>Read and display the number of bands</li>
#     <li>Display band types and shape</li>
#   </ul>
#
# </div>



# <div class="alert alert-info" style="font-size:110%">
#
#   <ul>
#     <li>Display the image in true color</li>
#   </ul>
#
# </div>



# ### Making a false color composition
#
# False color compositions with infrared are common to highlight elements in the infrared, in particular vegetation.
#
# ![Infrarouge.jpg](img/Infrarouge.jpg)  
# *Infrared image, vegetation particularly stands out.*

# <div class="alert alert-info" style="font-size:110%">
#
#   <ul>
#     <li>Open the near infrared (NIR) band, acquired at the same time (file : T31UDQ_20190705T105031_B08_10m.tif).</li>
#     <li>Compare the values of the pixels from the true color image and the one from the infrared band, using histograms.</li>
#     <li>Display a false color composition (NIR, Red, Green), using a normalization of the NIR band.</li>
#   </ul>
#
# </div>

# #### Exploring the NIR band



# #### Normalizing the NIR band

def normalize(rst):
    ...
    return normalized



# #### Creating the false color composite
#
# Note: you can use function `Raster.split_bands` to first create single layer Rasters, then `Raster.from_array` to generate the composite.



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



# ### Extracting vegetation zones

# <div class="alert alert-info" style="font-size:110%">
#
#   <ul>
#     <li>Apply a threshold to the NDVI to extract vegetated areas.</li>
#     <li>Save extracted areas as vector files.</li>
#   </ul>
#
# </div>

threshold_veg = 0.2
threshold_water = -0.9

# #### Case with 2 categories (True = vegetation, False = no vegetation)



# ### Case with 4 categories
# 0 = water \
# 1 = bare ground \
# 2 = vegetation \
# 10 = out of bounds pixels



# ### BONUS - Use folium to plot results


