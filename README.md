# xdem demo - Winter school, Nepal, December 2023

This repository contains a set of notebooks to showcase how to handle georeferenced (raster, vector) files and run DEM analysis in Python using [geoutils](https://github.com/GlacioHack/GeoUtils/) and [xdem](https://github.com/GlacioHack/xdem/).

**Running the examples from your computer**

The examples can be downloaded and run on your computer by installing the necessary packages with conda. Simply download the content of this repository ("code" button on the top-right or in command-line `git clone -b psf_nepal_2022 https://github.com/GlacioHack/xdem-demo.git`), then run: 
```
conda env create -f environment.yml
conda activate xdem-demo
jupyter notebook
```
The folder `notebooks` contains several examples:
- `intro_geoutils.ipynb` showcase the basic functionalities of geoutils and xdem, using sample data that comes with the code (the data is already downloaded on your computer).
- `exercice1.ipynb` aims at applying the same tools to a different data set, provided during the winter school. The solutions are provided in `exercice1_solution.ipynb`.
- `mass_balance_mera.ipynb` is an example on how to coregister two DEMs and a calculate a glacier geodetic mass balance using xdem. The solutions to the final exercise are provided in `mass_balance_mera_solutions.ipynb`.
