# Setup: how to run the tutorial

```{contents}
:local:
```

## a. Execute with a Binder deployment

Binder is a cloud service that allows users to share reproducible interactive computing environments from code repositories. It is generally used to enable other users to easily run your own code through Jupyter notebooks. It is a really cool service offered for free by several organisations (MyBinder through Jupyter, etc.).

Binder is probably the easiest way to execute notebooks in this repository, as you only have to do one click to arrive in a Jupyterlab with all the necessary libraries. However, its public resources are limited, supporting only some *'tens'* (not hundreds) of instances at a time with restricted processing power. The primary purpose of this functionality here; is to demonstrate the plug-and-play nature of *geo3D*. 

Each Notebook in the village and suburbs sections ([interactiveOnly](../notebooks/interactiveOnly), [LoD1_3D_CityModel](../notebooks/osm_LoD1_3DCityModel) and [CityJSONSpatialDataScience](../notebooks/CityJSONspatialDataScience)) have a Binder icon, from which to launch the **that particular** notebook on the Binder service.

Alternatively, to launch the **entire** `geo3D_wrkshp` GitHub repository; click on the Binder below button:

[![Binder](_static/launch-binder.svg)](https://mybinder.org/v2/gh/AdrianKriger/geo3D_wrkshp/HEAD)

Thereafter; navigate and choose a notebook, using the file browser, on the left side of the Jupyterlab screen.

## b. Execute on your own computer

Almost all parts of this tutorial were designed to run with limited computer resources, so it is possible to run on your laptop. It is not an easy one-click-push-button solution; as you will have to install the software environment yourself. 

Steps to run this tutorial on your own computer are listed below and demonstrated through Linux commands only:
```python
#- git clone the geo3D_wrkshp repository
git clone git clone https://github.com/adriankriger/geo3D_wrkshp.git
```
Install the required software environment with Conda. If you do not have Conda, install it by following these instructions (see here). Then create the environment, this can take a few minutes.
```python
conda env create -n geo3D_wrkshp -f geo3D_wrkshp/.binder/environment.yml
```
Launch a Jupyterlab notebook server from the `geo3D_wrkshp` environment.
```python
conda activate geo3D_wrkshp
jupyter lab
```
Open a web browser and connect to the Jupyterlab provided URL (you should see it in the jupyter lab command outputs), something like: `http://localhost:8888/lab?token=42fac6733c6854578b981bca3abf5152`
Navigate to geo3D_wrkshp/workshop/notebooks/ using the file browser on the left side of the Jupyterlab screen.