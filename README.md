# ESA training course, Innsbruck, September 2024

Ths repository contains the material for the tutorial on glaciers geodetic mass balance for ESA's 13th ADVANCED TRAINING COURSE ON LAND REMOTE SENSING (Snow and Glaciers), held in Innsbruck from 16 to 20 September.

It contains a set of Python notebooks to showcase how to handle georeferenced (raster, vector) files and run DEM analysis using [geoutils](https://github.com/GlacioHack/GeoUtils/) and [xdem](https://github.com/GlacioHack/xdem/).

**Running the examples from your computer**

- Installing mamba:
If you don't already have conda and/or mamba (both are Python package managers) installed on your computer, you first need to install them. The instructions can be found at [this link](https://github.com/conda-forge/miniforge).

- Download the code:
The examples can be downloaded and run on your computer by installing the necessary packages with conda. Simply download the content of this repository at [this link](https://github.com/GlacioHack/xdem-demo/archive/refs/heads/2024_esa.zip) ("code" button on the top-right of this page) or in command-line `git clone -b 2024_esa https://github.com/GlacioHack/xdem-demo.git`, then run:
```
mamba env create -f environment.yml
mamba activate xdem-demo
jupyter notebook
```
Note: you can replace the `mamba` commands with `conda` but mamba is faster.

- Download the data:
For Linux/MacOS users, the data can be downloaded with the script `download_data.sh` found with the code. Alternatively, you can download it manually at [this link](https://filesender.renater.fr/download.php?token=00e59d13-98dc-4b7f-9ffd-92bb28e06987&files_ids=43430563). You then need to unzip the file and copy the extracted `data` folder in the same folder as the one containing the codes.

**Running the examples from a Binder hub**
You can run the examples directly in your web browser without installing Python, using [Binder](https://jupyter.org/binder). There are two options (valid at least until 12/10/2024):
- myBinder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GlacioHack/xdem-demo.git/2024_esa) 
- UGA Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://binderhub.univ-grenoble-alpes.fr/v2/gh/GlacioHack/xdem-demo.git/2024_esa)

**Description of the notebooks**
The folder `notebooks` contains 3 notebooks:
- `intro_geoutils_demo.ipynb` showcases the basic functionalities of geoutils and xdem, using sample data that comes with the code (the data is already downloaded on your computer). The solutions are provided in `intro_geoutils_demo_solutions.ipynb`.
- `intro_geoutils_exercise.ipnb` is another example similar to `intro_geoutils_demo.ipynb` with a different dataset to let you get experience with geoutils/xdem.
- `mass_balance_mera.ipynb` is an example on how to coregister two DEMs and a calculate a glacier geodetic mass balance using xdem. The solutions to the final exercise are provided in `mass_balance_mera_solutions.ipynb`.
