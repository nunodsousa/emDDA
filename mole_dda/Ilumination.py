#
#
# Illumination
#
# Nuno de Sousa

import numpy as np


class Ilumination(object):
    """
    Ilumination
    ===========
    """

    E_0i = None  # Background Electromagnetic incident wave (usually a plane wave)

    def plane_wave(self, k, E0_const, positions, k_inc='z', return_f=True):
        """
        Plane Wave generator. This function generates the field experienced by a particle at some position in space.
        The wave must have a constant that characterizes the electric field (E0_const), the wave number and wave vector
        (k and k_inc).

        :param k: (float) Wavenumber
        :param positions: (np.array) positions of the particle
        :param E0_const:(float) Electric Field constant
        :param k_inc: direction of the wave vector
        :return: (np.array) Electric field in each particle with two different polarizations (sorted by x, y, z)
        """

        N_particles = positions.shape[0]

        if (k_inc == 'z'):
            kvector = k * np.array([0, 0, 1])
            Ex = np.zeros((N_particles, 3), dtype=complex)
            Ey = np.zeros((N_particles, 3), dtype=complex)
            i, x_i = np.mgrid[0:N_particles, 0:3]
            j, y_j = np.mgrid[0:N_particles, 0:3]
            Ex[np.where((x_i == 0))] = E0_const * np.exp(1j * np.array([np.dot(kvector, i) for i in positions]))
            Ey[np.where((y_j == 1))] = E0_const * np.exp(1j * np.array([np.dot(kvector, j) for j in positions]))
            Ex = Ex.reshape(3 * N_particles)
            Ey = Ey.reshape(3 * N_particles)
            E = np.transpose(np.vstack([Ex, Ey]))

        if (k_inc == 'y'):
            kvector = k * np.array([0, 1, 0])
            Ex = np.zeros((N_particles, 3), dtype=complex)
            Ez = np.zeros((N_particles, 3), dtype=complex)
            i, x_i = np.mgrid[0:N_particles, 0:3]
            j, z_j = np.mgrid[0:N_particles, 0:3]
            Ex[np.where((x_i == 0))] = E0_const * np.exp(1j * np.array([np.dot(kvector, i) for i in positions]))
            Ez[np.where((z_j == 2))] = E0_const * np.exp(1j * np.array([np.dot(kvector, j) for j in positions]))
            Ex = Ex.reshape(3 * N_particles)
            Ez = Ez.reshape(3 * N_particles)
            E = np.transpose(np.vstack([Ex, Ez]))

        if (k_inc == 'x'):
            kvector = k * np.array([1, 0, 0])
            Ey = np.zeros((N_particles, 3), dtype=complex)
            Ez = np.zeros((N_particles, 3), dtype=complex)
            i, y_i = np.mgrid[0:N_particles, 0:3]
            j, z_j = np.mgrid[0:N_particles, 0:3]
            Ex[np.where((y_i == 1))] = E0_const * np.exp(1j * np.array([np.dot(kvector, i) for i in positions]))
            Ey[np.where((z_j == 2))] = E0_const * np.exp(1j * np.array([np.dot(kvector, j) for j in positions]))
            Ey = Ey.reshape(3 * N_particles)
            Ez = Ez.reshape(3 * N_particles)
            E = np.transpose(np.vstack([Ey, Ez]))

        self.E_0i = E
        if (return_f == True):
            return self.E_0i