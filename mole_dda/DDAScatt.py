###########################################################################
# DDAScatt library
#
# Nuno de Sousa
# December 2019
###########################################################################

import numpy as np
from .CrossSections import CrossSections
from .Ilumination import Ilumination
from .Solvers import Solver
import numba as nb


class DDA(CrossSections, Ilumination, Solver):
    """
    DDA main class

    :param CrossSections class
    :param Ilumination class
    :param Solver class
    """

    SGreenTensor = None

    def __init__(self):
        pass

    def SGreenTensorE_compute(self, positions, k, alpha, epsilon_m=1, method='numpy', return_f=True):
        if (method == 'numpy'):
            print("Using numpy to load matrix.")
            np.seterr(divide='ignore', invalid='ignore')

            N_particles = positions.shape[0]

            G = np.zeros((N_particles, 3, N_particles, 3), dtype=complex)

            i, x_i, j, x_j = np.mgrid[0:N_particles, 0:3, 0:N_particles, 0:3]
            G[np.where((i == j) & (x_i == x_j))] = 1
            R = np.linalg.norm(positions[None, :, :] - positions[:, None, :], axis=-1)
            R = R.reshape(N_particles, 1, N_particles, 1)

            r = positions[None, :, :] - positions[:, None, :]

            krq = (k * R) ** 2
            pf = -k ** 2 / (epsilon_m) * alpha * np.exp(1j * k * R) / (4 * np.pi * R)
            a = 1. + (1j * k * R - 1.) / (krq)
            b = (3. - 3. * 1j * k * R - krq) / (krq)

            comb_r = (r[:, :, :, None] * r[:, :, None, :]).transpose([0, 2, 1, 3])
            G = pf * (b * comb_r / (R ** 2))
            G[np.where(x_i == x_j)] = (G + pf * a)[np.where(x_i == x_j)]
            G[np.where(i == j)] = 0
            G[np.where((i == j) & (x_i == x_j))] = 1

            self.SGreenTensor = G.reshape(N_particles * 3, N_particles * 3)

            if (return_f == True):
                return self.SGreenTensor

        elif (method == 'numba'):
            print('Using numba to load matrix.')
            N_particles = positions.shape[0]
            SGreenTensor = np.eye(N_particles*3,dtype=complex)

            #def SGreenTensor_numba(self, SGreenTensor, N_particles, alpha, positions, k):
            self.SGreenTensor = self.SGreenTensor_numba(SGreenTensor, N_particles, alpha, positions, k)

            #self.SGreenTensor_numba(SGreenTensor, N_particles, alpha, positions, k)

            if(return_f == True):
               return self.SGreenTensor

    @staticmethod
    @nb.jit(nopython=True)
    def SGreenTensor_numba(G_tensor, N_particles, alpha, pos, k):
        for i in range(N_particles):
            for j in range(N_particles): #Do not use symmetrization
                if (i != j):
                    # Do lots of things, here is shown an example.
                    # However you should not be scared because
                    # it only fills the G_tensor
                    # R = np.linalg.norm(np.array(pos[i])-np.array(pos[j]))
                    rx = pos[i][0] - pos[j][0]
                    ry = pos[i][1] - pos[j][1]
                    rz = pos[i][2] - pos[j][2]
                    R = np.sqrt(rx ** 2 + ry ** 2 + rz ** 2)
                    krq = (k * R) ** 2
                    pf = -k ** 2 * alpha * np.exp(1j * k * R) / (4 * np.pi * R)
                    a = 1. + (1j * k * R - 1.) / (krq)
                    b = (3. - 3. * 1j * k * R - krq) / (krq)
                    G_tensor[3 * i, 3 * j] = pf * (a + b * (rx * rx) / (R ** 2))  # Gxx
                    G_tensor[3 * i + 1, 3 * j + 1] = pf * (a + b * (ry * ry) / (R ** 2))  # Gyy
                    G_tensor[3 * i + 2, 3 * j + 2] = pf * (a + b * (rz * rz) / (R ** 2))  # Gzz
                    G_tensor[3 * i, 3 * j + 1] = pf * (b * (rx * ry) / (R ** 2))  # Gxy
                    G_tensor[3 * i, 3 * j + 2] = pf * (b * (rx * rz) / (R ** 2))  # Gxz
                    G_tensor[3 * i + 1, 3 * j] = pf * (b * (ry * rx) / (R ** 2))  # Gyx
                    G_tensor[3 * i + 1, 3 * j + 2] = pf * (b * (ry * rz) / (R ** 2))  # Gyz
                    G_tensor[3 * i + 2, 3 * j] = pf * (b * (rz * rx) / (R ** 2))  # Gzx
                    G_tensor[3 * i + 2, 3 * j + 1] = pf * (b * (rz * ry) / (R ** 2))  # Gzy

        return G_tensor