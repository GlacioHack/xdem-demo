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

fn_dem_2012 = "../data/mb_Mera/rasters/Mera_Pleiades_2012-11-25_DEM_4m.tif"
fn_dem_2018 = "../data/mb_Mera/rasters/Mera_Pleiades_2018-10-28_DEM_4m.tif"
fn_ref_dem = "../data/mb_Mera/rasters/Mera_COP30_DEM_UTM45_data.tif"
rgi_shpfile = "../data/mb_Mera/RGI_shapefiles/Glacier_inventory_around_Mera.shp"
mera_shpfile_2012 = "../data/mb_Mera/glacier_outlines/Mera_outline_2012_realigned.shp"
mera_shpfile_2018 = "../data/mb_Mera/glacier_outlines/Mera_outline_2018_realigned.shp"

# #### Load all DEMs at once, cropping to the common extent, and reprojecting onto the 2018 DEM grid

dem_2012, dem_2018, ref_dem = gu.raster.load_multiple_rasters([fn_dem_2012, fn_dem_2018, fn_ref_dem], crop=True, ref_grid=1)


# #### Load glacier outlines

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

inlier_mask = ~gl_mask.data.data & slope_mask & outlier_mask
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


# ### What about uncertainties?
#
# No systematic uncertainties, as the co-registration artificially set the median elevation difference to zero on the stable terrain. One can assume random uncertainties to have three independent sources: elevation change uncertainty ($\sigma _{\Delta z}$), delineation uncertainty ($\sigma _A$) and volume-to-mass density conversion uncertainty ($\sigma _\rho$).
#
# Here, we will assume that $\sigma _A$, the delineation error, is the glacier perimeter multiplied by two times the resolution of the satellite images used for the delineation (10 m). $\sigma _\rho$ is the different density uncertainties, i.e. $60 kg. m^{−3}$ for multi-year mass balances (Huss, 2013 - 10.5194/tc-7-877-2013).
#
# For ($\sigma _{\Delta z}$) we make use of the stable terrain to quantify the elevation change's uncertainty.
#
# As a first (crude) approximation, and assuming that the glacier is larger than the squared decorrelation length of the random errors (Brun et al., 2016 - 10.1038%2FNGEO2999), one could simply take the standard deviation of the elevation change on the off-glacier terrain:

# first approximation of uncertainty of elevation change
sigma_dz = nmad_coreg # https://xdem.readthedocs.io/en/stable/intro_robuststats.html
print(f"After coregistration NMAD:{sigma_dz:.2f}")

# The uncertainty on the volume estimate can then be obtained from (Berthier et al., 2014 - 10.5194/tc-8-2275-2014):
#
# $$\sigma _{\Delta V} = \sqrt((\sigma _{\Delta z}(p+5*(1-p))*A_{mean})^2+(\sigma _A \Delta z)^2) $$
#
# Where p is the proportion of non-voided area (voids are penalized by a factor 5), $A_{mean}$ is the mean area between the 2 DEM dates and $\Delta z$ is the mean on-glacier elevation change.

# +
# first approximation of uncertainty of volume change
p = 1-np.sum(dh_coreg[mera_mask].mask)/dh_coreg[mera_mask].size

print(f"proportion of non-voided area p:{p:.2f}")

mera_mask_2012 = mera_outlines_2012.create_mask(dh_coreg)
mera_mask_2018 = mera_outlines_2018.create_mask(dh_coreg)

A_mean = (np.sum(mera_mask_2012)+np.sum(mera_mask_2018))*dh_coreg.res[0]**2/2 # in m2

print(f"mean area between the 2 DEM dates:{A_mean:.2f} m\u00b2")

sigma_A = 2*22000*10 # in m2 - ~22km perimeter times resolution of sensor used for delineation (10 m). Other studies (Brun et al., 2016) take 10% of the area

print(f"delineation error:{sigma_A:.2f} m\u00b2")

dz = np.nanmean(dh_coreg[mera_mask])

print(f"mean on-glacier elevation change:{dz:.2f} m")

sigma_dv = np.sqrt((sigma_dz*(p+(1-p))*A_mean)**2+(sigma_A*dz)**2)

print(f"Total volume change:{dV:.2f} +/- {sigma_dv:.2f} m\u00b3")
# -

# Then the uncertainty on the mass change would be given by:
#
# $$\sigma _{\Delta M} = \sqrt((\sigma _{\Delta V}*\rho)^2+(\sigma _\rho \Delta V)^2) $$

# +
# uncertainty in mass change

sigma_dM = np.sqrt((sigma_dv*0.85)**2+(0.06*dV)**2)

print(f"Total mass change: {dM/1e3:.2f} +/- {sigma_dM/1e3:.2f} t")

# in m w.e.

sigma_dh_mwe = sigma_dM/mean_area

print(f"Specific mass balance: {dh_mwe:.1f} +/- {sigma_dh_mwe:.1f} m w.e.")
# -

# ### Actually, there are some strong assumptions behing what we did above...
#
# And while it gives a reasonable first guess, one would want to account for the heteroscedasticity and spatial correlation of the elevation errors... some general ideas below, but if you want to explore this in more details, best is to refer to the very complete xdem documentation (https://xdem.readthedocs.io/), xdem advanced tutorials and the study by Hugonnet et al. (2022): https://doi.org/10.1109/jstars.2022.3188922
#
# Propagating elevation errors spatially accounting for heteroscedasticity and spatial correlation is complex. It requires computing the pairwise correlations between all points of an area of interest (be it for a sum, mean, or other operation), which is computationally intensive. Here, we rely on published formulations to perform computationally-efficient spatial propagation for the mean of elevation (or elevation differences) in an area.
#
# References: Hugonnet et al. (2022), Figure S16, Equations 17–19 and Rolstad et al. (2009), Equation 8 (http://dx.doi.org/10.3189/002214309789470950).
#
# We load the same data, and perform the same calculations on heteroscedasticity and spatial correlations of errors as in the Elevation error map and Spatial correlation of errors examples.

slope, maximum_curvature = xdem.terrain.get_terrain_attribute(dem_2018, attribute=["slope", "maximum_curvature"])
errors, df_binning, error_function = xdem.spatialstats.infer_heteroscedasticity_from_stable(
    dvalues=dh_coreg, list_var=[slope, maximum_curvature], list_var_names=["slope", "maxc"], unstable_mask=rgi_outlines
)

# We use the error map to standardize the elevation differences before variogram estimation, following Equation 12 of Hugonnet et al. (2022), which is more robust as it removes the variance variability due to heteroscedasticity.

zscores = dh / errors
emp_variogram, params_variogram_model, spatial_corr_function = xdem.spatialstats.infer_spatial_correlation_from_stable(
    dvalues=zscores, list_models=["Gaussian", "Spherical"], unstable_mask=rgi_outlines, random_state=42
)

# With our estimated heteroscedasticity and spatial correlation, we can now perform the spatial propagation of errors. The best estimation of the standard error for Mera Glacier is done by directly providing the shapefile, which relies on Equation 18 of Hugonnet et al. (2022).

areas=[gu.Vector(mera_shpfile_2012).ds[gu.Vector(mera_shpfile_2012).ds["Area"]>0],]
areas

# +
areas=[gu.Vector(mera_shpfile_2012).ds[gu.Vector(mera_shpfile_2012).ds["Area"]>0],]
stderr_gla = xdem.spatialstats.spatial_error_propagation(
    areas=areas, errors=errors, params_variogram_model=params_variogram_model
)

print(f"The error (1-sigma) in mean elevation change for Mera is {stderr_gla:.2f} meters.")
# -

# When passing a numerical area value, we compute an approximation with disk shape from Equation 8 of Rolstad et al. (2009). This approximation is practical to visualize changes in elevation error when averaging over different area sizes, but is less accurate to estimate the standard error of a certain area shape.

areas = 10 ** np.linspace(1, 12)
stderrs = xdem.spatialstats.spatial_error_propagation(
    areas=areas, errors=errors, params_variogram_model=params_variogram_model
)
plt.plot(areas / 10**6, stderrs)
plt.xlabel("Averaging area (km²)")
plt.ylabel("Standard error (m)")
plt.vlines(
    x=np.pi * params_variogram_model["range"].values[0] ** 2 / 10**6,
    ymin=np.min(stderrs),
    ymax=np.max(stderrs),
    colors="red",
    linestyles="dashed",
    label="Disk area with radius the\n1st correlation range of {:,.0f} meters".format(
        params_variogram_model["range"].values[0]
    ),
)
plt.vlines(
    x=np.pi * params_variogram_model["range"].values[1] ** 2 / 10**6,
    ymin=np.min(stderrs),
    ymax=np.max(stderrs),
    colors="blue",
    linestyles="dashed",
    label="Disk area with radius the\n2nd correlation range of {:,.0f} meters".format(
        params_variogram_model["range"].values[1]
    ),
)
plt.xscale("log")
plt.legend()
plt.show()

# ## It's your turn !

# Test running the code with other parameters, here are some ideas:
# 1. Calculate the mass balance using the mean value of each bin instead of the median.
# 2. Try with different elevation bin sizes
# 3. Try removing pixels whose value differ by more than five NMAD from the mean\
# How sensitive are your results to these different parameters?\
# Can your reproduce the results of Wagnon et al. (2021), JoG. Check the methods to see what bin size and filtering method they used. Your results should match their +/- 0.03 m w.e..
# 4. Try calculating the mass balance for another glacier.\
# Hint 1: To select one glacier from the RGI outlines, use `rgi_outlines.ds[rgi_outlines.ds.RGIId == "YOUR_RGI_ID"]` where YOUR_RGI_ID is replaced by the ID of your chosen glacier, within the Pleiades DEM bounds.


