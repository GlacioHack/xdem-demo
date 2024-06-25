# xdem demo - S3 workshop, Grenoble, June 2024

This repository contains a set of notebooks to showcase how to handle georeferenced (raster, vector) files and run DEM analysis in Python using [geoutils](https://github.com/GlacioHack/GeoUtils/) and [xdem](https://github.com/GlacioHack/xdem/).
They are split into beginner/advanced depending on whether or not you have previous experience with geoutils/xdem.

**Running the examples from your computer**

- Download the code:
The examples can be downloaded and run on your computer by installing the necessary packages with conda. Simply download the content of this repository ("code" button on the top-right or in command-line `git clone -b 2024_ige https://github.com/GlacioHack/xdem-demo.git`), then run:
```
mamba env create -f environment.yml
mamba activate xdem-demo
jupyter notebook
```
Note: you can replace the `mamba` commands with `conda` but mamba is faster.

- Download the data:
The data can be downloaded with the script `download_data.sh` or directly at [this link](https://filesender.renater.fr/?s=download&token=bc73932e-4999-4237-8554-2e1a8862fcd2) (you then need to unzip the file and copy it in the same folder as the one containing the codes).

**Running the examples from a Binder hub**
You can run the examples directly in your web browser without installing Python, using [Binder](https://jupyter.org/binder). There are two options (valid at least until 29/06/2024):
- myBinder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GlacioHack/xdem-demo.git/2024_ige) 
- UGA Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://binderhub.univ-grenoble-alpes.fr/v2/gh/GlacioHack/xdem-demo.git/2024_ige) 

**Description of the notebooks**
The folder `notebooks` contains several examples. The examples are sorted for beginners or advanced users:
- `intro_geoutils.ipynb` showcase the basic functionalities of geoutils, using sample data that comes with the code (the data is already downloaded on your computer). The solutions are provided in `intro_geoutils_solutions.ipynb` if appropriate.
- `mass_balance_mera.ipynb` is an example on how to coregister two DEMs and a calculate a glacier geodetic mass balance using xdem. The solutions to the final exercise are provided in `mass_balance_mera_solutions.ipynb`.
- a third notebook `makina/vegetation_paris.ipyng` further showcases raster analysis not related to DEMs.
