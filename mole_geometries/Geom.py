#
# Simple geometry library
#
# Nuno de Sousa
# December 2019

import numpy as np

def sphere_discretization(N_edge, radius = 1, spatial_translation = [0,0,0]):
    """
    sphere_discretization
    =====================
    Fully vectorized discretization of a sphere.
    description: It generates a cube N_edge elements per edge. We then select the elements with center inside of a
    virtual sphere with radius = radius. In the end we can replace the sphere at some other point in space with
    center = spatial_translation.
    :param N_edge: Number of particles in the edge of the cube
    :param radius: Radius of the virtual sphere used to select the particles.
    :param spatial_translation: Position of the center of the sphere.
    :return: numpy.array(positions), int(Number of particles), float(volume of each element)
    """
    x = np.linspace(0, N_edge - 1, N_edge)
    y = np.linspace(0, N_edge - 1, N_edge)
    z = np.linspace(0, N_edge - 1, N_edge)

    X, Y, Z = np.meshgrid(x, y, z)

    X = X.reshape((np.prod(X.shape),))
    Y = Y.reshape((np.prod(Y.shape),))
    Z = Z.reshape((np.prod(Z.shape),))

    coords = list(zip(X, Y, Z))
    coords = (radius / (N_edge) + 2 * radius / N_edge * np.array(coords)) - radius

    sphere = coords[np.linalg.norm(coords, axis=1) <= radius] + spatial_translation
    N_particles = sphere.shape[0]
    element_volume = (2 * radius / N_edge) ** 3
    return sphere, N_particles, element_volume

def random_cube(N_particles, edge_x = 1, edge_y = 1, edge_z = 1, spatial_translation = [0,0,0]):
    """
    Generate N_particles randomly uniform distributed inside of a box with size (edge_x, edge_y, edge_z).
    As option, it is possible to make a spatial translation of the box giving a spatial_translation vector.

    It returns the positions of the particles
    """
    spatial_translation = np.array(spatial_translation)

    edges = np.array([edge_x, edge_y, edge_z]);

    positions = np.random.rand(N_particles,3)*edges - edges/2 + spatial_translation

    return positions
