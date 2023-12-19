# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # DEM coregistration and glacier geodetic mass balance with xdem
#
# ### This notebook showcases how to use xdem to coregister two DEMs and calculate the glacier geodetic mass balance of Mera glacier, Nepal, using two DEMs derived from Pleiades stereo images.
#
#

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
# All the data is taken from the Mera mass balance geodetic tutorial located in the "mb_Mera" folder.

# #### Get input data path ###
# ### <span style='color:red '> **TO DO:** </span> If needed, update the path to the files below.
#

# #### Load all DEMs at once, cropping to the common extent, and reprojecting onto the 2018 DEM grid

fn_dem_2012 = "../data/04_mb_Mera/rasters/Mera_Pleiades_2012-11-25_DEM_4m.tif"
fn_dem_2018 = "../data/04_mb_Mera/rasters/Mera_Pleiades_2018-10-28_DEM_4m.tif"
fn_ref_dem = "../data/04_mb_Mera/rasters/Mera_COP30_DEM_UTM45_data.tif"
dem_2012, dem_2018, ref_dem = gu.raster.load_multiple_rasters([fn_dem_2012, fn_dem_2018, fn_ref_dem], crop=True, ref_grid=1)


# #### Load glacier outlines

rgi_shpfile = "../data/04_mb_Mera/RGI_shapefiles/Glacier_inventory_around_Mera.shp"
mera_shpfile_2012 = "../data/04_mb_Mera/glacier_outlines/Mera_outline_2012_realigned.shp"
mera_shpfile_2018 = "../data/04_mb_Mera/glacier_outlines/Mera_outline_2018_realigned.shp"
rgi_outlines = gu.Vector(rgi_shpfile)
mera_outlines_2012 = gu.Vector(mera_shpfile_2012)
mera_outlines_2018 = gu.Vector(mera_shpfile_2018)

# ### Calculate and plot elevation change

dh = dem_2018 - dem_2012


# #### Plot original dh map with glacier contours

vmax=30
plt.figure(figsize=(10, 10))
ax = plt.subplot(111)
rgi_outlines.show(ax=ax, facecolor='none', edgecolor='k', lw=0.5, zorder=2)
mera_outlines_2012.show(ax=ax, facecolor='none', edgecolor='k', zorder=3)
dh.show(ax=ax, cmap='RdYlBu', vmin=-vmax, vmax=vmax, cbar_title='Elevation change 2012 - 2018 (m)', zorder=1)
ax.set_title('Mera glacier and surroundings')
plt.tight_layout()
plt.show()

# There are non-zeros elevation changes outside of glaciers => These DEMs need to be coregistered.

# ## 2 - DEM coregistration

# #### Prepare inputs for coregistration
# First we create a mask, i.e. a raster of same shape as our dh map, to mask pixels on glaciers. `gl_mask` is `True` on glaciers, `False` elsewhere.

gl_mask = rgi_outlines.create_mask(dh)

# Then we mask pixels in steep slopes and gross blunders

slope = xdem.terrain.slope(ref_dem)
slope_mask = (slope.data < 40).filled(False)
outlier_mask = (np.abs(dh.data) < 50).filled(False)

# We plot the final mask of pixels used for coregistration

inlier_mask = ~gl_mask.data.filled(False) & slope_mask & outlier_mask
plt.figure(figsize=(8, 8))
plt.imshow(inlier_mask.squeeze())
plt.show()

# Free memory (needed when running on binder)

del slope, slope_mask, outlier_mask

# ### Run coregistration
# #### We use the Nuth & Kaab (2011) algorithm to estimate a horizontal offset, then we remove a possible vertical bias by removing the median dh value in stable terrain.

# First, we need to define the coregistration function. \
# Then we estimate the coregistration needed between our two Pleiades DEMs. \
# Finally, we apply that coregistration to the 2012 DEM.

coreg = xdem.coreg.NuthKaab() + xdem.coreg.VerticalShift(vshift_func=np.median)
coreg.fit(dem_2018, dem_2012, inlier_mask, verbose=True)
dem_2012_coreg = coreg.apply(dem_2012)

# The results of the coregistration can be output like this

print(coreg.pipeline[0]._meta)

# ### <span style='color:red '> **Question:** </span> How many iterations of the Nuth & Kaab algorithm were run?
#

# #### Answer: ...

# ### Calculate new elevation change and plot

dh_coreg = dem_2018 - dem_2012_coreg

# +
plt.figure(figsize=(10, 6))
ax1 = plt.subplot(121)
mera_outlines_2012.show(ax=ax1, facecolor='none', edgecolor='k', zorder=3)
dh.show(ax=ax1, cmap='RdYlBu', vmin=-vmax, vmax=vmax, cbar_title='Elevation change 2012 - 2018 (m)', zorder=1)
ax1.set_title('Before coregistration')

ax2 = plt.subplot(122)
mera_outlines_2012.show(ax=ax2, facecolor='none', edgecolor='k', zorder=3)
dh_coreg.show(ax=ax2, cmap='RdYlBu', vmin=-vmax, vmax=vmax, cbar_title='Elevation change 2012 - 2018 (m)', zorder=1)
ax2.set_title('After coregistration')

plt.tight_layout()
plt.show()
# -

# #### We see that elevation changes outside glaciers are close to zero (yellow color).

# ### Calculate statistics of before/after coregistration
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
# This is particularly needed for data with gaps, as the values in `ddem_bins['count']` are the number of pixels with observations, not total pixel count.
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


# ## It's your turn !

# Test running the code with other parameters, here are some ideas:
# 1. Calculate the mass balance using the mean value of each bin instead of the median.
# 2. Try with different elevation bin sizes
# 3. Try removing pixels whose value differ by more than five NMAD from the mean\
# How sensitive are your results to these different parameters?\
# Can your reproduce the results of Wagnon et al. (2021), JoG. Check the methods to see what bin size and filtering method they used. Your results should match their +/- 0.05 m w.e..
# 4. Try calculating the mass balance for another glacier.\
# Hint 1: To select one glacier from the RGI outlines, use `rgi_outlines.ds[rgi_outlines.ds["rgi_id"] == "YOUR_RGI_ID"]` where YOUR_RGI_ID is replaced by the ID of your chosen glacier, within the Pleiades DEM bounds.

# 1. Using the mean

ddem_bins2 = xdem.volume.hypsometric_binning(dh_coreg[mera_mask], ref_dem[mera_mask], bins=25, aggregation_function=np.mean)
dV2 = np.sum(ddem_bins2['value'] * bins_area)  # in m^3
dM2 = 0.85 * dV2  # in kg
dh_mwe2 = dM2 / mean_area
print(dh_mwe2)

# 2. Using different bin sizes

ddem_bins4 = xdem.volume.hypsometric_binning(dh_coreg[mera_mask], ref_dem[mera_mask], bins=10, aggregation_function=np.mean)
bins_area4 = xdem.volume.calculate_hypsometry_area(ddem_bins4, ref_dem[mera_mask], pixel_size=dh_coreg.res)
dV4 = np.sum(ddem_bins4['value'] * bins_area4)  # in m^3
dM4 = 0.85 * dV4  # in kg
dh_mwe4 = dM4 / mean_area
print(dh_mwe4)

# 3. Using a 5 NMAD filter

dh_mera = dh_coreg[mera_mask]
dh_mera_filt = dh_mera.copy()
dh_mera_filt[np.abs(dh_mera - np.mean(dh_mera)) > 5 * xdem.spatialstats.nmad(dh_mera)] = np.ma.masked
ddem_bins3 = xdem.volume.hypsometric_binning(dh_mera_filt, ref_dem[mera_mask], bins=25, aggregation_function=np.mean)
dV3 = np.sum(ddem_bins3['value'] * bins_area)  # in m^3
dM3 = 0.85 * dV3  # in kg
dh_mwe3 = dM3 / mean_area
print(dh_mwe3)

# Matching Wagnon et al. (2021)

ddem_bins5 = xdem.volume.hypsometric_binning(dh_mera_filt, ref_dem[mera_mask], bins=10, aggregation_function=np.mean)
bins_area5 = xdem.volume.calculate_hypsometry_area(ddem_bins5, ref_dem[mera_mask], pixel_size=dh_coreg.res)
dV5 = np.sum(ddem_bins5['value'] * bins_area5)  # in m^3
dM5 = 0.85 * dV5  # in kg
dh_mwe5 = dM5 / mean_area
print(dh_mwe5)

# #### The results of Wagnon et al. (2021) are −2.55 ± 0.34 m w.e. for the period 2012-2018. We get a slightly different result, probably because of different implementations in the coregistration and/or binning. 

# Try calculating the mass balance for another glacier
# Hint: To select one glacier from the RGI outlines, use `rgi_outlines.ds[rgi_outlines.ds["rgi_id"] == "RGI2000-v7.0-C-15-04940"]` 
unknown_gl_outline = gu.Vector(rgi_outlines.ds[rgi_outlines.ds["rgi_id"] == "RGI2000-v7.0-C-15-04940"])
unknown_gl_outline = unknown_gl_outline.reproject(dh)
unknown_gl_mask = unknown_gl_outline.create_mask(dh)

# Plot dh map and outline

fig, ax = plt.subplots()
dh_coreg.show(ax=ax, cmap='RdYlBu', vmin=-vmax, vmax=vmax)
unknown_gl_outline.show(ax=ax, facecolor='none', edgecolor='k')
plt.show()

# Calculate hypsometric bins

ddem_bins = xdem.volume.hypsometric_binning(dh_coreg[unknown_gl_mask], ref_dem[unknown_gl_mask], bins=25)
bins_area = xdem.volume.calculate_hypsometry_area(ddem_bins, ref_dem[unknown_gl_mask], pixel_size=dh_coreg.res)

# ### Plot the results

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

# ### Calculate outputs

dV = np.sum(ddem_bins['value'] * bins_area)  # in m^3
dM = 0.85 * dV  # in kg
mean_area = float(unknown_gl_outline.ds.area_km2.iloc[0])
dh_mwe = dM / mean_area / 1e6

print(f"Total volume change: {dV:.2f} m\u00b3")
print(f"Total mass change: {dM/1e3:.2f} t")
print(f"Specific mass balance: {dh_mwe:.1f} m w.e.")

# ### NB: The example shown here is very simple and does not take into account outliers or gaps that could be present in the data. A more thorough analysis of the result should always be carried out !

