# Setup: how to run the tutorial

```{contents}
:local:
```

## Execute with a Binder deployment

Binder is a cloud service that allows users to share reproducible interactive computing environments from code repositories. It is generally used to enable other users to easily run your own code through Jupyter notebooks. It is a really cool service offered for free by several organisations (MyBinder through Jupyter, etc.).

Binder is probably the easiest way to execute notebooks in this repository, as you only have to do one click to arrive in a Jupyterlab with all the necessary libraries. However, the hardware resources of the public Binder are limited; it is not meant for more than 50 instances at a time. Processing power is also limited to 

Each Notebook in the interactiveOnly, LoD1 3D City Model and Spatial Data Science sections, for example have a rocket icon 🚀 at the top, from which you can select the Binder button to just launch the particular Notebook on the Binder service.

Alternatively, you can also directly click on the Binder below button:

Binder

You can then navigate and choose a Notebook, using the file browser, on the left side of the Jupyterlab screen.

## Execute on your own computer

Most parts of this tutorial were designed to run with limited computer resources, so it is possible to run on your laptop. It is a bit more complicated as you will have to install the software environment yourself. 

Steps to run this tutorial on your own computer are listed below and demonstrated through Linux commands only:
```python
git clone the .....
git clone https://github.com/.....git
```
Install the required software environment with Conda. If you do not have Conda, install it by following these instructions (see here). Then create the environment, this can take a few minutes.
```python
conda env create -n geo3D_wrkshp -f geo3D_wrkshp/.binder/environment.yml
```
Launch a Jupyterlab notebook server from this environment.
```python
conda activate geo3D_wrkshp
jupyter lab
```
Open a web browser and connect to the Jupyterlab provided URL (you should see it in the jupyter lab command outputs), something like: http://localhost:8888/lab?token=42fac6733c6854578b981bca3abf5152.
Navigate to foss4g-2022/tutorial/pangeo101/ using the file browser on the left side of the Jupyterlab screen.