# ESA training course, Innsbruck, September 2024

Ths repository contains the material for the tutorial on glaciers geodetic mass balance for ESA's 13th ADVANCED TRAINING COURSE ON LAND REMOTE SENSING (Snow and Glaciers), held in Innsbruck from 16 to 20 September.

It contains a set of Python notebooks to showcase how to handle georeferenced (raster, vector) files and run DEM analysis using [geoutils](https://github.com/GlacioHack/GeoUtils/) and [xdem](https://github.com/GlacioHack/xdem/).

**Running the examples from your computer**

- Installing mamba or anaconda (it should already be installed on your machine):
If you don't already have conda and/or mamba (both are Python package managers) installed on your computer, you first need to install them. The instructions can be found at [this link](https://github.com/conda-forge/miniforge).

- Download the code:
The examples can be downloaded and run on your computer by installing the necessary packages with conda. Simply download the content of this repository at [this link](https://github.com/GlacioHack/xdem-demo/archive/refs/heads/2024_esa.zip) ("code" button on the top-right of this page) or in command-line `git clone -b 2024_esa https://github.com/GlacioHack/xdem-demo.git`. Using the conda console, change your directory to the newly downloaded folder, at the location of the file 'environment.yml' then run:
```
conda env create -f environment.yml
conda activate xdem-demo
jupyter notebook
```
- Data:
For the second and third notebooks you will need the data that is stored in the 'data' folder from the 'GeodeticMB_Kneib' folder that you downloaded at the start of the workshop. Copy this folder and paste it in the same location as the 'notebooks' folder.

**Running the examples from a Binder hub**
You can run the examples directly in your web browser without installing Python, using [Binder](https://jupyter.org/binder). Just use the following link :

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GlacioHack/xdem-demo.git/2024_esa) (valid at least until 12/10/2024)

**Description of the notebooks**
The folder `notebooks` contains 3 notebooks:
- `intro_geoutils_demo.ipynb` showcases the basic functionalities of geoutils and xdem, using sample data that comes with the code (the data is already downloaded on your computer). 
- `intro_geoutils_exercise.ipynb` is another example similar to `intro_geoutils_demo.ipynb` with a different dataset to let you get experience with geoutils/xdem. The solutions are provided in `intro_geoutils_solutions.ipynb`.
- `mass_balance_mera_extended.ipynb` is an example on how to coregister two DEMs and a calculate a glacier geodetic mass balance using xdem. 
