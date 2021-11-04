# How to install

Using the terminal, we can use `pip to install the package. pip: `pip install git+https://github.com/nunodsousa/emDDA.git` . 


# emDDA

A computational library for computing electromagnetic scatteringby particles with arbitrary shape and (electric/magnetic) polarizability.

The goal of this code is to have a general high performance library that allow us to compute the scattering properties of electromagnetic particles.

Two libraries are going to be builded. 

The *first library* contains a set of functions that allow us to generate different geometries. Also it allow us to have random systems.

The *second library* constains the electromagnetic library.

## What should we have

1 - A Super-Green tensor generator with Electric, Magnetic and Electric-Magnetic Green Tensors.
1.1 - Compute the numerical and analytical Green Tensor derivative.

2 - The execution can be made through MKL and tensorflow.

3 - Tensorial polarizability (3x3 or 6x6, depending of the Super-Green Tensor)

4 - The medium can be different from vaccum.

5 - Different light sources (several can be included simulataneously) are considered.

6 - Quantities to compute:
  - Cross Sections
  - Emission Decay Rates
  - Electromagnetic Fields
  - Electromagnetic Forces
  - Field Maps
  
## Technical documentation

Link to the PDF with the technical description. 

## Sintax

Here we must define the Sintax of the library.

## Tutorials/Examples

1 - Computation of the cross sections for a Si spherical particle. Representation of the Fields map.

2 - Computation of the decay rate of a random system with exclusion volume.

3 - Computation the magneto-optical activity of a oblate particle.

4 - 

## Luis Froufe proposal

Multi-Layered Magneto Optical Discrete Dipole Approach

This is a collaborative project to develop a numerical kernel based on the discrete dipole approach  (DDA). The goal of this software is provide a consistent, robust and versatile set of numerical routines based on DDA to solve a number of electromagnetic problems. The main capabilities are Electric and magnetic dipoles.

Any permittivity (in particular magneto-optical materials are allowed), (what about permeability?).

multilayered media.

2D periodicity is allowed. ``` NdS note: 3D systems with periodicity in two directions? ```

Any piece included in this software must fulfill a few requirements:

Fortran callable: Many prototypes have been already coded in Fortran03. Ideally, the libraries should be callable also from  C/C++, Matlab and Python. ``` NdS note: No way. We should converge into a single scripting language otherwise we will spend part of our time fixing compatibility problems. The core can be written in one or several languages, but the interface must be only one.```

portable code: The numerical kernel must run in Linux, Mac and Windows.

API: precise description of the API (application programmin interface) is mandatory for each routine.

Documentation: the code must be properly commented and auxiliary documents showing in detail the basis of each algorithm.

# To start using emDDA

## Environments, conda and jupyter instalation
The following commands work exactly the same way in linux and MacOS. In windows I never tried, but I suppose that it is the same. The goal of this subsection is to briefly explain how to use conda environments.

1 - It is necessary to install anaconda (https://anaconda.org). Anaconda is a free distribution of Python for scientic computing. The main advantage is the simplicity of the instalation. With the basic instalation you have most of the necessary libraries. For more info: https://en.wikipedia.org/wiki/Anaconda_(Python_distribution)

2 - Create a conda environment. Conda environments are very similar to Python’s virtual environments. Both serve to help manage dependencies and isolate projects, and they function in a similar way, with one key distinction: conda environments are language agnostic. That is, they support languages other than Python. If you're interested, take a look to this medium article: https://towardsdatascience.com/a-guide-to-conda-environments-bc6180fc533.
A detailed list of conda commands is available here.
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

If you are in a rush, here are the main things you might want to do it.

2.1 - The first step is to create a conda environment. Just open a terminal and write: `conda create -n myenvironmentname numpy scipy matplotlib tensorflow termcolor python=3.7`. This will create a conda environment with *python 3.7* with  libraries *matplotlib*, *tensorflow*, *termcolor*, *scipy* and *numpy*. It will install also many other libraries, such as the *MKL*. You should replace the `myenvironmentname` by the name you want. In my case, I'm using `env_emDDA` as a name.

2.2 - Once we have the environment created, just activate the environment. In the terminal, you must write `conda activate myenvironmentname` and it is done,  i.e. you're inside of the environment. If you don't remember or you want to see the name of your environments in your system, just write `conda env list`.

2.3 - If you want to install another library such as  `statsmodels` in the `myenvironmentname` environment, just write in the terminal `conda install -c anaconda statsmodels`.

2.4 - If you want to deactivate an environment, just write `conda deactivate`

2.5 - To open a Jupyter notebook in the webbrowser inside of the environment: `jupyter notebook`. In my case, since I'm using MacOS Catalina, inside of the .zshrc I have included the `alias jup='jupyter notebook'`, so for me to start a jupyter notebook I just have to write `jup`.

And this is it! Now you are ready to use `emDDA`.

# emDDA library

XXXX
