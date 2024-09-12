# ESA training course, Innsbruck, September 2024

Ths repository contains the material for the tutorial on glaciers geodetic mass balance for ESA's 13th ADVANCED TRAINING COURSE ON LAND REMOTE SENSING (Snow and Glaciers), held in Innsbruck from 16 to 20 September.

It contains a set of notebooks to showcase how to handle georeferenced (raster, vector) files and run DEM analysis in Python using [geoutils](https://github.com/GlacioHack/GeoUtils/) and [xdem](https://github.com/GlacioHack/xdem/).

**Running the examples from your computer**

- Download the code:
The examples can be downloaded and run on your computer by installing the necessary packages with conda. Simply download the content of this repository ("code" button on the top-right or in command-line `git clone -b 2024_esa https://github.com/GlacioHack/xdem-demo.git`), then run:
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
The folder `notebooks` contains 3 notebooks:
- `intro_geoutils_demo.ipynb` showcases the basic functionalities of geoutils and xdem, using sample data that comes with the code (the data is already downloaded on your computer). The solutions are provided in `intro_geoutils_demo_solutions.ipynb`.
- `intro_geoutils_exercise.ipnb` is another example similar to `intro_geoutils_demo.ipynb` with a different dataset to let you get experience with geoutils/xdem.
- `mass_balance_mera.ipynb` is an example on how to coregister two DEMs and a calculate a glacier geodetic mass balance using xdem. The solutions to the final exercise are provided in `mass_balance_mera_solutions.ipynb`.
