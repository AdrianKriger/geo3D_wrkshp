# Libraries and Packages

While *geo3D* uses the essential python libraries and packages typically used in all Spatial Data Analysis; such as `numpy`, `matplotlib`, `scipy`, `shapely`, `geopandas` and `gdal`; the following are not:

i. `Triangle `

A LoD1 3D Model is based on a triangulated irregular network (TIN) data-structure. A TIN is a collection of connected 3D triangles that form a continuous closed surface. osm_LoD1_3DCityModel creates a TIN with the python implementation of Jonathan Richard Shewchuk’s Triangle and uses a special algorithm called a Constrained Delaunay triangulation {cite}`Dzhelil2020,shewchuk96b`.

ii. `pydeck`

The interactive pseudo-3D visualisation is via pydeck. Pydeck is a set of Python bindings for making spatial visualisations with deck.gl, a WebGL-powered framework developed by Uber for visual exploratory data analysis of large datasets, optimised for a Jupyter environment {cite}`Visgl2020`.

iii. `Plus Codes`

A Plus Code is a geocode based on a system of regular grids for identifying an area anywhere on the Earth. Plus Codes are designed to be used like street addresses and may be especially useful in places where there is no formal system to identify buildings, such as street names, house numbers, and post codes. Plus Codes are therefore an alternative to traditional addresses and are extremely helpful to under severed communities who cannot access basic services such as receive postal, emergency, and social services {cite}`Goo2015`.

iv. `cjio`

A python command line interface (CLI) to process and manipulate CityJSON files {cite}`cjio2022`.
