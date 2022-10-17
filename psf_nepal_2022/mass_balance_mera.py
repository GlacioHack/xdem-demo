# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Calculating a glacier geodetic mass balance with xdem
#
# ### This notebook showcases how to use xdem to calculate the glacier geodetic mass balance of Mera glacier, Nepal, using two DEMs derived from Pleiades stereo images.
#
#

# ### Import the necessary modules

# +
import matplotlib.pyplot as plt
import numpy as np

import geoutils as gu
import xdem

# To avoid interpolation in plt.imshow
plt.rcParams['image.interpolation'] = 'none'
# -


# ## 1 - Load input data ###

# #### Get input data path ###

fn_dem_2012 = "./data/mb_Mera/2012-11-25_DEM_4m.tif"
fn_dem_2018 = "./data/mb_Mera/2018-10-28_DEM_4m.tif"
fn_ref_dem = "./data/mb_Mera/Mera_Co30_DEM_UTM45.tif"
rgi_shpfile = "./data/mb_Mera/shapefile_RGI/Glacier_inventory_around_Mera.shp"
mera_shpfile_2012 = "./data/mb_Mera/shapefiles/Mera_outline_2012_realigned.shp"
mera_shpfile_2018 = "./data/mb_Mera/shapefiles/Mera_outline_2018_realigned.shp"

# #### Load all DEMs

dem_2012 = xdem.DEM(fn_dem_2012)
dem_2018 = xdem.DEM(fn_dem_2018)
ref_dem = xdem.DEM(fn_ref_dem)

# #### Crop and reproject DEMs on a same projection/grid

bounds = gu.projtools.merge_bounds([dem_2012, dem_2018, ref_dem], merging_algorithm="intersection")
dem_2018.crop(bounds)
dem_2012 = dem_2012.reproject(dst_ref=dem_2018)
ref_dem = ref_dem.reproject(dst_ref=dem_2018)

# #### Load glacier outlines

rgi_outlines = gu.Vector(rgi_shpfile)
mera_outlines_2012 = gu.Vector(mera_shpfile_2012)
mera_outlines_2018 = gu.Vector(mera_shpfile_2018)

# ### Calculate and plot elevation change

dh = dem_2018 - dem_2012


# #### Plot original dh map with glacier contours

vmax=30
ax = plt.subplot(111)
rgi_outlines.ds.plot(ax=ax, facecolor='none', edgecolor='k', lw=0.5, zorder=2)
mera_outlines_2012.ds.plot(ax=ax, facecolor='none', edgecolor='k', zorder=3)
dh.show(ax=ax, cmap='RdYlBu', vmin=-vmax, vmax=vmax, cb_title='Elevation change 2012 - 2018 (m)', zorder=1)
ax.set_title('Mera glacier and surroundings')
plt.tight_layout()

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

inlier_mask = ~gl_mask & slope_mask & outlier_mask
plt.imshow(inlier_mask.squeeze())

# #### Run coregistration
# We use the Nuth & Kaab (2011) algorithm to estimate a horizontal offset, then we remove a possible vertical bias by removing the median dh value in stable terrain.

# First, we need to define the coregistration function. \
# Then we estimate the coregistration needed between our two Pleiades DEMs. \
# Finally, we apply that coregistration to the 2012 DEM.

coreg = xdem.coreg.NuthKaab() + xdem.coreg.BiasCorr(bias_func=np.median)
coreg.fit(dem_2018, dem_2012, inlier_mask, verbose=True)
dem_2012_coreg = coreg.apply(dem_2012)

# The results of the coregistration can be output like this

print(coreg.pipeline[0]._meta)

# #### Calculate new elevation change and plot

dh_coreg = dem_2018 - dem_2012_coreg

# +
plt.figure(figsize=(12, 6))
ax1 = plt.subplot(121)
mera_outlines_2012.ds.plot(ax=ax1, facecolor='none', edgecolor='k', zorder=3)
dh.show(ax=ax1, cmap='RdYlBu', vmin=-vmax, vmax=vmax, cb_title='Elevation change 2012 - 2018 (m)', zorder=1)
ax1.set_title('Before coregistration')

ax2 = plt.subplot(122)
mera_outlines_2012.ds.plot(ax=ax2, facecolor='none', edgecolor='k', zorder=3)
dh_coreg.show(ax=ax2, cmap='RdYlBu', vmin=-vmax, vmax=vmax, cb_title='Elevation change 2012 - 2018 (m)', zorder=1)
ax2.set_title('After coregistration')

plt.tight_layout()
# -

# We see that elevation changes outside glaciers are close to zero (yellow color).

# #### Calculate statistics of before/after coregistration
# Because `dh.data` is a masked array, we use `compressed()` to output only unmaked values.

# +
inlier_orig = dh.data[inlier_mask].compressed()
nstable_orig, mean_orig = len(inlier_orig), np.mean(inlier_orig)
med_orig, nmad_orig = np.median(inlier_orig), xdem.spatialstats.nmad(inlier_orig)
print(f"Number of stable pixels: {nstable_orig}")
print(f"Before coregistration:\
      \n\tMean dh: {mean_orig:.2f}\
      \n\tMedian dh: {med_orig:.2f}\
      \n\tNMAD dh: {nmad_orig:.2f}")

inlier_coreg = dh_coreg.data[inlier_mask].compressed()
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
# The calculation requires a (raster) mask of the glacier, the dh values, and a reference DEM used to get the elevation of the pixels. All must be on the same grid. \
# The default elevation bins have a height of 50 m, from the DEM minimum to maximum. Here we set the bin height to 25 m with the argument `bins`.

mera_mask = mera_outlines_2012.create_mask(dh_coreg)
ddem_bins = xdem.volume.hypsometric_binning(dh_coreg.data[mera_mask], ref_dem.data[mera_mask], bins=25, aggregation_function=np.median)

# The output is a panda's DataFrame containing the elevation bin edges, the mean dh ('value') and observations (valid pixels) count.

print(ddem_bins)

# ### Calculate the glacier area within each elevation bin
# This is particularly needed for data with gaps, as the values in `ddem_bins['count']` are the number of pixels with observations, not total pixel count.
# The result is in m$^2$.

bins_area = xdem.volume.calculate_hypsometry_area(ddem_bins, ref_dem.data[mera_mask], pixel_size=dh_coreg.res)
print(bins_area)

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
plt.plot(ddem_bins["value"], ddem_bins.index.mid, "k")
plt.xlabel("Elevation change (m/a)")

plt.tight_layout()
# -

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
# Can your reproduce the results of Wagnon et al. (2021), JoG. Check the methods to see what bin size and filtering method they used. Your results should match their +/- 0.03 m w.e..
# 4. Try calculating the mass balance for another glacier.\
# Hint 1: To select one glacier from the RGI outlines, use `rgi_outlines.ds[rgi_outlines.ds.RGIId == "YOUR_RGI_ID"]` where YOUR_RGI_ID is replaced by the ID of your chosen glacier, within the Pleiades DEM bounds.
