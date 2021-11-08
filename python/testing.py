import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

from mastu_exhaust_analysis.read_efit import read_uda
from mastu_exhaust_analysis.fluxsurface_tracer import trace_flux_surface
from pyEquilibrium.equilibrium import equilibrium


def get_field_line_coords_demo():
    """
    Taken from JRH example in mastu_exhasut_analysis/fluxsurface_tracer.py
    """
    # Example script to show how to run fluxsurface_tracer

    efit_data = read_uda('45100', calc_bfield=True)
    i = 100
    start_r = np.linspace(1.40, 1.45, 10)
    start_z = 0.0
    ax = plt.subplot()

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    surfaces_list = []
    for j in np.arange(len(start_r)):
        surfaces_list.append(trace_flux_surface(efit_data, i, start_r[j], start_z, cut_surface=True))
        # ax.plot(surfaces_list[j][0].r, surfaces_list[j][0].z, 'r')

    a = 5
    # ax.plot(surfaces_list[0][0].wall.xy[0], surfaces_list[0][0].wall.xy[1], 'k')
    # ax.set_aspect(1.0)
    plt.show()


def testing_pyequilibrium():

    # Load a MAST shot
    eq = equilibrium(device='MAST', shot=45000, time=0.25)

    a = 5


# def testing_pyfieldlinetracer():
#
#
#     self._psi = getdata("efm_psi(r,z)", shot)



def testing_read_uda():
    out = read_uda('45100', calc_bfield=False)
    r, z, t, psi, psin = out['r'], out['z'], out['t'], out['psi'], out['psiN']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.pcolormesh(r, z, psin[100, ...])
    cbar = fig.colorbar(im, ax=ax)
    ax.set_aspect('equal')
    plt.show()


if __name__ == '__main__':
    # get_field_line_coords_demo()
    testing_pyequilibrium()