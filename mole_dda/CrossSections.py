#
#
# Cross Sections library
#
#
# Nuno de Sousa

import numpy as np


class CrossSections(object):
    """
    CrossSections calculation

    It contains the functions:
    - sigma_ext_calc - Calculation of the Extintion Cross Section
    - sigma_scatt_calc - Calculation of the Scattering Cross Section
    - sigma_abs_calc - Calculation of the Absorption Cross Section
    - relative_error_calc - It computes the relative error

    It returns the sections in the same units as the input.
    """

    sigma_scatt = None
    sigma_ext = None
    sigma_abs = None

    def sigma_ext_calc(self, epsilon_m, k, alpha, E0_const, epsilon_0=1, return_f=False):
        """
        It Computes the Extintion Cross Section.

        :param: epsilon_m = dielectric constant of the medium
        :param: k = wavenumber in free space
        :param: alpha = polarizability of the particles (it is a scalar)
        :param: E0_const = Magnitud of the electric field
        :param: epsilon_0 = dielectric permitivity of vacuum.
        :param: return_f = (True or False), if True returns a vector with the sigma_ext, otherwise it only writes in the class.
        """
        p = epsilon_0 * epsilon_m * alpha * self.E_inc

        self.sigma_ext = np.diagonal(
            k / (epsilon_0 * epsilon_m * E0_const ** 2) * np.imag(np.dot(np.conjugate(np.transpose(self.E_0i)), p)))
        if (return_f == True):
            return self.sigma_ext

    def sigma_scatt_calc(self, epsilon_m, k, alpha, E0_const, epsilon_0=1, return_f=False):
        """
        It Computes the Scattering Cross Section.

        :param: epsilon_m = dielectric constant of the medium
        :param: k = wavenumber in free space
        :param: alpha = polarizability of the particles (it is a scalar)
        :param: E0_const = Magnitud of the electric field
        :param: n_particles = Number of particles in the system.
        :param: epsilon_0 = dielectric permitivity of vacuum.
        :return: return_f = (True or False), if True returns a vector with the sigma_scatt, otherwise it only writes in the class.
        """

        n_particles = int(self.SGreenTensor.shape[0] / 3)
        p = epsilon_0 * epsilon_m * alpha * self.E_inc

        # diagonal parts
        diag = np.diag((k / (6 * np.pi)) * np.real(np.dot(np.transpose(np.conjugate(p)), p)))

        # nondiagonal parts
        i, j = np.mgrid[0:3 * n_particles, 0:3 * n_particles]
        tempG = self.SGreenTensor.copy()
        tempG[np.where(i == j)] = 0
        out_diag = np.dot(np.transpose(np.conjugate(p)), np.dot(np.imag(-(epsilon_m / (k ** 2 * alpha)) * tempG), p))
        self.sigma_scatt = k ** 3 / (epsilon_0 * epsilon_m * E0_const) ** 2 * np.diagonal(np.real(out_diag + diag))

        if (return_f == True):
            return self.sigma_scatt

    def sigma_abs_calc(self, epsilon_m, k, alpha, E0_const, epsilon_0=1, return_f=False):
        """
        It Computes the Absorption Cross Section.

        :param: epsilon_m = dielectric constant of the medium
        :param: k = wavenumber in free space
        :param: alpha = polarizability of the particles (it is a scalar)
        :param: alpha0inv = Inverse of the static polarizability
        :param: E0_const = Magnitud of the electric field
        :return: return_f = (True or False), if True returns a vector with the sigma_ext, otherwise it only writes in the class.
        """

        alpha0inv = (1 / alpha + 1j * k ** 3 / (6 * np.pi))

        p = epsilon_m * alpha * self.E_inc
        self.sigma_abs = k / ((epsilon_m) ** 2 * E0_const ** 2) * np.diagonal(
            np.imag(np.dot(np.transpose(p), np.conjugate(alpha0inv * p))))

        if (return_f == True):
            return self.sigma_abs

    def relative_error_calc(self, return_f=False):
        """
        It checks the optical theorem
        """
        if ((self.sigma_abs is not None) & (self.sigma_ext is not None) & (self.sigma_scatt is not None)):
            self.relative_error = np.array((self.sigma_ext - self.sigma_scatt - self.sigma_abs) / self.sigma_ext)

        else:
            raise ValueError("One or several sections are not calculated.")

        if (return_f == True):
            return self.relative_error