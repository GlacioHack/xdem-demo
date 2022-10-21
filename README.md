# xdem demo - PSF Capable, Nepal, October 2022

This repository contains a set of notebooks to showcase how to handle georeferenced (raster, vector) files and run DEM analysis in Python using [geoutils](https://github.com/GlacioHack/GeoUtils/) and [xdem](https://github.com/GlacioHack/xdem/).

**Running the examples hassle-free, on a web browser**

You can run the examples directly in your web browser without installing Python, using [Binder](https://jupyter.org/binder). To do so, simply click on this button: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GlacioHack/xdem-demo/psf_nepal_2022) (valid at least until 01/11/2022)

Then navigate into the "psf_nepal_2022" folder and open any of the examples. The python scripts can be converted in Jupyter notebooks with right-click -> "Open with..." -> "Notebook" (thanks to [jupytext](https://github.com/mwouts/jupytext)).

The repository contains two scripts/notebooks:
- `intro_geoutils.py` - shows how to open, manipulate and save raster and vector files with xdem and geoutils.
- `mera_mass_balance.py` - explains how to calculate the mass balance of a glacier (here Mera glacier in Nepal) from two DEMs. The data needed to run this notebook is automatically downloaded in the folder "data". It can also be downloaded from [this link](https://filesender.renater.fr/download.php?token=5aab3c13-cdfb-41f5-a930-8da6368e2659&archive_format=undefined&files_ids=18517548) (valid until 01/11/2022). 

**Running the examples from your computer**

The examples can also be downloaded and run on your computer by installing the necessary packages with conda. Simply download the content of this repository ("code" button on the top-right or in command-line `git clone -b psf_nepal_2022 https://github.com/GlacioHack/xdem-demo.git`), then run: 
```
conda env create -f environment.yml
conda activate xdem-demo
jupyter notebook
```
