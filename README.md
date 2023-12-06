# If you're using it...
Please cite the paper:
https://www.nature.com/articles/srep30803

de Sousa, N., Froufe-Pérez, L., Sáenz, J. et al. Magneto-Optical Activity in High Index Dielectric Nanoantennas. Sci Rep 6, 30803 (2016). https://doi.org/10.1038/srep30803

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
