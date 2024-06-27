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

fn_dem_2012 = "../../data/mb_Mera/rasters/Mera_Pleiades_2012-11-25_DEM_4m.tif"
fn_dem_2018 = "../../data/mb_Mera/rasters/Mera_Pleiades_2018-10-28_DEM_4m.tif"
fn_ref_dem = "../../data/mb_Mera/rasters/Mera_COP30_DEM_UTM45_data.tif"
dem_2012, dem_2018, ref_dem = gu.raster.load_multiple_rasters([fn_dem_2012, fn_dem_2018, fn_ref_dem], crop=True, ref_grid=1)


# #### Load glacier outlines: RGI outlines, Mera 2012 and Mera 2018 outlines

rgi_shpfile = "../../data/mb_Mera/RGI_shapefiles/Glacier_inventory_around_Mera.shp"
mera_shpfile_2012 = "../../data/mb_Mera/glacier_outlines/Mera_outline_2012_realigned.shp"
mera_shpfile_2018 = "../../data/mb_Mera/glacier_outlines/Mera_outline_2018_realigned.shp"
rgi_outlines = gu.Vector(rgi_shpfile)
mera_outlines_2012 = gu.Vector(mera_shpfile_2012)
mera_outlines_2018 = gu.Vector(mera_shpfile_2018)

# ### Calculate and plot elevation change

# #### Calculate the 2012-2018 elevation change

dh = dem_2018 - dem_2012


# #### Plot original dh map with glacier contours

# First we need to reproject the RGI outlines in the same CRS as the DEMs
rgi_outlines = rgi_outlines.reproject(dem_2012)

vmax=30
plt.figure(figsize=(10, 10))
ax = plt.subplot(111)
rgi_outlines.plot(ax=ax, facecolor='none', edgecolor='k', lw=0.8, zorder=2)
mera_outlines_2012.plot(ax=ax, facecolor='none', edgecolor='turquoise', zorder=3)
dh.plot(ax=ax, cmap='RdYlBu', vmin=-vmax, vmax=vmax, cbar_title='Elevation change 2012 - 2018 (m)', zorder=1)
ax.set_title('Mera glacier and surroundings')
plt.tight_layout()
plt.show()

# There are non-zeros elevation changes outside of glaciers => These DEMs need to be coregistered.

# ## 2 - DEM coregistration

# #### Prepare inputs for coregistration
# First we create a mask, i.e. a raster of same shape as our dh map, to mask pixels on glaciers. `gl_mask` is `True` on glaciers, `False` elsewhere. \
# Since the RGI outlines can be slightly inaccurate, it is recommended to use a 100 m buffer.

# gl_mask = rgi_outlines.create_mask(dh)
gl_mask = rgi_outlines.buffer(100).create_mask(dh)

# Then we mask pixels in steep slopes (> 40 degrees) and gross blunders (abs(dh) > 50).

slope = xdem.terrain.slope(ref_dem)
slope_mask = (slope < 40)
outlier_mask = (np.abs(dh) < 50)

# We plot the final mask of pixels used for coregistration.

inlier_mask = ~gl_mask & slope_mask & outlier_mask

plt.figure(figsize=(8, 8))
inlier_mask.plot()
plt.show()

# To avoid issues with some numpy functions later (np.median), convert to numpy array with nodata set to False

inlier_mask = inlier_mask.data.filled(False)

# Free memory (needed when running on binder)

del slope, slope_mask, outlier_mask

# ### Run coregistration
# #### We use the Nuth & Kaab (2011) algorithm to estimate a horizontal offset, then we remove a possible vertical bias by removing the median dh value in stable terrain.

# First, we need to define the coregistration function. \
# Then we estimate the coregistration needed between our two Pleiades DEMs. \
# Finally, we apply that coregistration to the 2012 DEM.

coreg = xdem.coreg.NuthKaab() + xdem.coreg.VerticalShift(vshift_reduc_func=np.median)
# coreg = xdem.coreg.NuthKaab() + xdem.coreg.Deramp()
coreg.fit(dem_2018, dem_2012, inlier_mask, verbose=True)
dem_2012_coreg = coreg.apply(dem_2012)

# #### Print the output of the coregistration: offset in east and north direction, vertical offset.

print(coreg.pipeline[0].meta)

# ### <span style='color:red '> **Question:** </span> How many iterations of the Nuth & Kaab algorithm were run?
#

# #### Answer: ...

# ### Calculate new elevation change and plot

dh_coreg = dem_2018 - dem_2012_coreg

# +
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
mera_outlines_2012.plot(ax=axes[0], facecolor='none', edgecolor='k', zorder=3)
dh.plot(ax=axes[0], cmap='RdYlBu', vmin=-vmax, vmax=vmax, cbar_title='Elevation change 2012 - 2018 (m)', 
        title="Before coregistration", zorder=1)

mera_outlines_2012.plot(ax=axes[1], facecolor='none', edgecolor='k', zorder=3)
dh_coreg.plot(ax=axes[1], cmap='RdYlBu', vmin=-vmax, vmax=vmax, cbar_title='Elevation change 2012 - 2018 (m)', 
        title="After coregistration", zorder=1)

plt.tight_layout()
plt.show()

fig, ax = plt.subplots(1, 1)
rgi_outlines.plot(ax=ax, facecolor='none', edgecolor='k', lw=0.5, zorder=2)
mera_outlines_2012.plot(ax=ax, facecolor='none', edgecolor='k', zorder=3)
dh_coreg.plot(ax=ax, cmap='RdYlBu', vmin=-4, vmax=4, cbar_title='Elevation change 2012 - 2018 (m)', 
        title="After coreg - detailed", zorder=1)
plt.tight_layout()
plt.show()
# -


# #### We see that elevation changes outside glaciers are close to zero (yellow color).

# ### Calculate statistics of before/after coregistration
#
# Print the number of inlier pixels, mean, median and NMAD of elevation change before and after coregistration.
#
# Because `dh.data` is a masked array, we use `compressed()` to output only unmasked values.

# +
inlier_orig = dh[inlier_mask]
nstable_orig, mean_orig = len(inlier_orig), np.mean(inlier_orig)
med_orig, nmad_orig = np.median(inlier_orig), xdem.spatialstats.nmad(inlier_orig)
print(f"Number of stable pixels: {nstable_orig}")
print(f"Before coregistration:\
      \n\tMean dh: {mean_orig:.2f}\
      \n\tMedian dh: {med_orig:.2f}\
      \n\tNMAD dh: {nmad_orig:.2f}")

inlier_coreg = dh_coreg[inlier_mask]
nstable_coreg, mean_coreg = len(inlier_coreg), np.mean(inlier_coreg)
med_coreg, nmad_coreg = np.median(inlier_coreg), xdem.spatialstats.nmad(inlier_coreg)
print(f"After coregistration:\
      \n\tMean dh: {mean_coreg:.2f}\
      \n\tMedian dh: {med_coreg:.2f}\
      \n\tNMAD dh: {nmad_coreg:.2f}")
# -

# ### Plot statistics of median elevation change as a function of slope, aspect and reference elevation

# First, we calculate slope, aspect and reference elevation for pixels in inlier mask.

slope, aspect = xdem.terrain.get_terrain_attribute(ref_dem, ["slope", "aspect"])

slope_inliers = slope[inlier_mask]
aspect_inliers = aspect[inlier_mask]
z_inliers = ref_dem[inlier_mask]

# We calculate elevation change before and after coregistration in inlier mask.


dh_inliers = dh[inlier_mask]
dh_coreg_inliers = dh_coreg[inlier_mask]

# We use the function `xdem.spatialstats.nd_binning` to calculate statistics (default is count, median, NMAD) for bins of each category, for elevation change before coregistration.

dh_stats = {}
dh_stats["slope"] = xdem.spatialstats.nd_binning(dh_inliers, [slope_inliers,], list_var_names=["slope",])
dh_stats["aspect"] = xdem.spatialstats.nd_binning(dh_inliers, [aspect_inliers,], list_var_names=["aspect",])
dh_stats["elevation"] = xdem.spatialstats.nd_binning(dh_inliers, [z_inliers,], list_var_names=["elevation",])


# The results are pandas dataframes containing the statistics for each bin

dh_stats["slope"]

# We do the same for elevation change after coregistration.

dh_coreg_stats = {}
dh_coreg_stats["slope"] = xdem.spatialstats.nd_binning(dh_coreg_inliers, [slope_inliers,], list_var_names=["slope",])
dh_coreg_stats["aspect"] = xdem.spatialstats.nd_binning(dh_coreg_inliers, [aspect_inliers,], list_var_names=["aspect",])
dh_coreg_stats["elevation"] = xdem.spatialstats.nd_binning(dh_coreg_inliers, [z_inliers,], list_var_names=["elevation",])


# And finally we make plots summarizing it for before/after coregistration.

fig, axes = plt.subplots(1, 3, figsize=(11, 4))
xdem.spatialstats.plot_1d_binning(dh_stats["slope"], "slope", "nanmedian", ax=axes[0])
xdem.spatialstats.plot_1d_binning(dh_stats["aspect"], "aspect", "nanmedian", ax=axes[1])
xdem.spatialstats.plot_1d_binning(dh_stats["elevation"], "elevation", "nanmedian", ax=axes[2])
plt.suptitle("Before correction")
plt.tight_layout()
plt.show()

fig, axes = plt.subplots(1, 3, figsize=(11, 4))
xdem.spatialstats.plot_1d_binning(dh_coreg_stats["slope"], "slope", "nanmedian", ax=axes[0])
xdem.spatialstats.plot_1d_binning(dh_coreg_stats["aspect"], "aspect", "nanmedian", ax=axes[1])
xdem.spatialstats.plot_1d_binning(dh_coreg_stats["elevation"], "elevation", "nanmedian", ax=axes[2])
plt.suptitle("After correction")
plt.tight_layout()
plt.show()


# ### BONUS: Calculate statistics as a function of curvature and along-track coordinates

# +
# Test used to roughly estimate the along-track angle, which maximizes the amplitude of along-track error
# for angle in np.arange(-26, -20, 1):
#    xx, yy = gu.raster.get_xy_rotated(raster=dh, along_track_angle=angle)
#    xx_inliers = xx[inlier_mask]
#    stats = xdem.spatialstats.nd_binning(dh_coreg_inliers, [xx_inliers,], list_var_names=["along-track",])
#    print(angle, stats["nanmedian"].max() - stats["nanmedian"].min())
#    #xdem.spatialstats.plot_1d_binning(stats, "along-track", "nanmedian")
#    #plt.show()
# -


curvature = xdem.terrain.get_terrain_attribute(ref_dem, ["curvature"])
curvature_inliers = curvature[inlier_mask]
xx, yy = gu.raster.get_xy_rotated(raster=dh, along_track_angle=-23)
yy_inliers = yy[inlier_mask]

dh_stats["curvature"] = xdem.spatialstats.nd_binning(dh_inliers, [curvature_inliers,], list_var_names=["curvature",], 
                                                     list_var_bins=[np.arange(-2, 2, 0.2)])
dh_stats["alongtrack"] = xdem.spatialstats.nd_binning(dh_inliers, [yy_inliers,], list_var_names=["alongtrack",])

dh_coreg_stats["curvature"] = xdem.spatialstats.nd_binning(dh_coreg_inliers, [curvature_inliers,], list_var_names=["curvature",], 
                                                           list_var_bins=[np.arange(-2, 2, 0.2)])
dh_coreg_stats["alongtrack"] = xdem.spatialstats.nd_binning(dh_coreg_inliers, [yy_inliers,], list_var_names=["alongtrack",])

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
xdem.spatialstats.plot_1d_binning(dh_stats["curvature"], "curvature", "nanmedian", ax=axes[0])
xdem.spatialstats.plot_1d_binning(dh_stats["alongtrack"], "alongtrack", "nanmedian", ax=axes[1])
plt.suptitle("Before correction")
plt.tight_layout()
plt.show()

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
xdem.spatialstats.plot_1d_binning(dh_coreg_stats["curvature"], "curvature", "nanmedian", ax=axes[0])
xdem.spatialstats.plot_1d_binning(dh_coreg_stats["alongtrack"], "alongtrack", "nanmedian", ax=axes[1])
plt.suptitle("After correction")
plt.tight_layout()
plt.show()

# ## 3 - Geodetic mass balance of Mera glacier
#
# ### Calculate mean dh over glaciers, per elevation bins
# We use the function `xdem.volume.hypsometric_binning` (see brief documentation [here](https://xdem.readthedocs.io/en/stable/gen_modules/xdem.volume.html)) to calculate the mean (or median) value of dh within elevation bands. \
# The calculation requires a (raster) mask of the glacier, the dh values, and a reference DEM used to get the elevation of the pixels. All must be on the same grid. \
# The default elevation bins have a height of 50 m, from the DEM minimum to maximum. Here we set the bin height to 25 m with the argument `bins`.

mera_mask = mera_outlines_2012.create_mask(dh_coreg)
ddem_bins = xdem.volume.hypsometric_binning(dh_coreg[mera_mask], ref_dem[mera_mask], bins=25, aggregation_function=np.median)

# The output is a panda's DataFrame containing the elevation bin edges, the mean dh ('value') and observations (valid pixels) count.

print(ddem_bins)

# ### Calculate the glacier area within each elevation bin
# We use the function `xdem.volume.calculate_hypsometry_area`. \
# This is particularly needed for data with gaps, as the values in `ddem_bins['count']` are the number of pixels with observations, not total pixel count. \
# The result is in m$^2$.

bins_area = xdem.volume.calculate_hypsometry_area(ddem_bins, ref_dem[mera_mask], pixel_size=dh_coreg.res)
print(bins_area)

# ### <span style='color:red '> **Questions:** </span> 
# - what is the median dh at 5400 m ?
# - at which elevation is the highest thinning?
# - what is the area of glaciers in the elevation band at 6000 m?
#

# ### Plot the results
# 1. Plot glacier area vs elevation
# 2. Plot the elevation change vs elevation

# +
plt.rcParams['font.size'] = 14
plt.figure(figsize=(8, 8))

plt.barh(
    y=ddem_bins.index.mid,
    width=bins_area/1e6,
    left=0,
    height=(ddem_bins.index.left - ddem_bins.index.right) * 1,
    alpha=0.8,
)
plt.xlabel("Glacier area per elevation bins (km\u00b2)")
plt.ylabel("Elevation (m)")

plt.twiny()
plt.plot(ddem_bins["value"].values, ddem_bins.index.mid.values, "k")
plt.xlabel("Elevation change (m/a)")

plt.tight_layout()
plt.show()
# -

# The glacier area is plotted with blue bars (x axis at the bottom) and elevation change in dark line (x axis at the top) versus elevation (y axis).

# ### Calculate total volume and mass change
# - volume change is in m$^3$
# - volume is converted to mass (kg) assuming a volume to mass conversion factor of 850 kg/m$^3$ (Huss, 2013, TC)

dV = np.sum(ddem_bins['value'] * bins_area)  # in m^3
dM = 0.85 * dV  # in kg

# ### Calculate specific mass balance
# In meter water equivalent (m w.e.)
# We divide by the average area between the two dates.

mean_area = float((mera_outlines_2012.ds.Area + mera_outlines_2018.ds.Area) / 2)
dh_mwe = dM / mean_area

# ### Print results

print(f"Total volume change: {dV:.2f} m\u00b3")
print(f"Total mass change: {dM/1e3:.2f} t")
print(f"Specific mass balance: {dh_mwe:.1f} m w.e.")

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

# Prepare data for plotting: reproject, convert to nan array and get WGS84 bounds
# Downsample dh map to 20 m for speed
dh_reproj = dh_coreg.reproject(crs="ESRI:53004", res=20)
bounds = dh_reproj.get_bounds_projected("EPSG:4326")
dh_nan = dh_reproj.get_nanarray()

# Create a colormap to be used
cmap = plt.get_cmap("RdYlBu")

# +
import folium

# Initialize map with ds.explore (easier) and plot glacier contours
m = rgi_outlines.ds.explore(style_kwds={"fill":False, "color": "black"}, name="RGI outlines", tiles="Esri.WorldImagery") #"area_km2", legend=True)

# Add raster layer with Folium
# Add normalization for the color scale
img = folium.raster_layers.ImageOverlay(
        name="dh",
        image=dh_nan,
        bounds=[[bounds.bottom, bounds.left], [bounds.top, bounds.right]],
        colormap=lambda x: cmap(x / 30. / 2 + 0.5),
        opacity=1,
    )
m.add_child(img)
folium.LayerControl().add_to(m)

# Display
m
# -

m.save("Mera_map.html")
