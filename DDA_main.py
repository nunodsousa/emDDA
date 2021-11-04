################################################################################################
#
# Development Version
#
#
# Nuno de Sousa
# December 2019
################################################################################################

import mole_geometries.Geom as Geom
import mole_dda.DDAScatt as DDAScatt
import mole_mie.Materials as Materials
import numpy as np
import time

if __name__ == "__main__":


    n_particles = 2000
    lam = 1250
    wavenumber = 2*np.pi/lam

    # Geometry generation
    pos = Geom.random_cube(n_particles, edge_x=1000, edge_y=1000, edge_z=1000)

    # particles' polarizability
    alpha, alpha0inv = Materials.polarization(volume = 5, k = wavenumber, epsilon_part = 12.25, epsilon_medium = 1)

    # DDA
    dda = DDAScatt.DDA()
    start_time = time.perf_counter()
    G1 = dda.SGreenTensorE_compute(positions=pos, k=wavenumber, epsilon_m=1, alpha=alpha, method='numba', return_f=True)
    #G2 = dda.SGreenTensorE_compute(positions=pos, k=wavenumber, epsilon_m=1, alpha=alpha, method='numpy', return_f=True)
    end_time = time.perf_counter()
    print(end_time - start_time, "seconds")

    dda.plane_wave(k=wavenumber, E0_const=1, positions=pos, k_inc='z', return_f=False)
    dda.direct_solver(lib = "tensorflow", return_f=False)

    # Cross Sections
    dda.sigma_ext_calc(epsilon_m=1, k=wavenumber, alpha=alpha, E0_const=1, return_f=False)
    dda.sigma_scatt_calc(epsilon_m=1, k=wavenumber, alpha=alpha, E0_const=1, epsilon_0=1, return_f=False)
    dda.sigma_abs_calc(epsilon_m=1, k=wavenumber, alpha=alpha, E0_const=1, epsilon_0=1, return_f=False)

    print('sigma_ext = ', dda.sigma_ext)
    print('sigma_scatt = ', dda.sigma_scatt)
    print('sigma_abs = ', dda.sigma_abs)
    print('relative error = ', dda.relative_error_calc(return_f=True))