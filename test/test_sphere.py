# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 09:36:46 2015

@author: dap124
"""

from __future__ import print_function

import os.path as osp

import matplotlib.pyplot as plt
import numpy as np
from helpers import read_1d_complex, write_1d_complex
from numpy.testing import assert_allclose

import openmodes
from openmodes.basis import DivRwgBasis
from openmodes.constants import c, eta_0
from openmodes.operator import CfieOperator, EfieOperator, MfieOperator
from openmodes.sources import PlaneWaveSource

tests_location = osp.split(__file__)[0]
mesh_dir = osp.join(tests_location, "input", "test_sphere")
reference_dir = osp.join(tests_location, "reference", "test_sphere")


def sphere_extinction_analytical(freqs, r):
    """Analytical expressions for a PEC sphere's extinction for plane wave
    with E = 1V/m

    Parameters
    ----------
    freqs : ndarray
        Frequencies at which to calculate
    r : real
        Radius of sphere
    """
    from scipy.special import spherical_yn

    N = 40

    k0r = freqs * 2 * np.pi / c * r
    # scs = np.zeros(len(k0r))
    # scs_modal = np.zeros((len(k0r), N))
    ecs_kerker = np.zeros(len(k0r))

    for count, x in enumerate(k0r):
        jn, jnp, yn, ynp = sph_jnyn(N, x)
        h2n = jn - 1j * yn
        h2np = jnp - 1j * ynp
        a_n = ((x * jnp + jn) / (x * h2np + h2n))[1:]
        b_n = (jn / h2n)[1:]
        # scs[count] = 2*np.pi*sum((2*np.arange(1, N+1)+1)*(abs(a_n)**2 + abs(b_n)**2))/x**2 #
        # scs_modal[count] = 2*np.pi*(2*np.arange(1, N+1)+1)*(abs(a_n)**2 + abs(b_n)**2)/x**2 #
        ecs_kerker[count] = (
            2
            * np.pi
            * np.real(np.sum((2 * np.arange(1, N + 1) + 1) * (a_n + b_n)))
            / x**2
        )

    return ecs_kerker * r**2 / eta_0


def test_extinction_all(
    plot_extinction=False, skip_asserts=True, write_reference=False
):
    "Extinction of a PEC sphere with EFIE, MFIE, CFIE"

    tests = (
        ("EFIE", EfieOperator, "extinction_efie.npy"),
        ("MFIE", MfieOperator, "extinction_mfie.npy"),
        ("CFIE", CfieOperator, "extinction_cfie.npy"),
    )

    for operator_name, operator_class, reference_filename in tests:
        print(operator_name)

        sim = openmodes.Simulation(
            name="horseshoe_extinction",
            basis_class=DivRwgBasis,
            operator_class=operator_class,
        )

        radius = 5e-3
        sphere = sim.load_mesh(osp.join(mesh_dir, "sphere.msh"))
        sim.place_part(sphere)

        num_freqs = 2
        freqs = np.linspace(5e9, 20e9, num_freqs)

        extinction = np.empty(num_freqs, np.complex128)

        e_inc = np.array([1, 0, 0], dtype=np.complex128)
        k_hat = np.array([0, 0, 1], dtype=np.complex128)
        pw = PlaneWaveSource(e_inc, k_hat)

        for freq_count, s in sim.iter_freqs(freqs):
            Z = sim.impedance(s)
            V = sim.source_vector(pw, s)
            V_E = sim.source_vector(pw, s, extinction_field=True)
            extinction[freq_count] = np.vdot(V_E, Z.solve(V))
            print(extinction[freq_count])

        extinction_filename = osp.join(reference_dir, reference_filename)

        if write_reference:
            # generate the reference extinction solution
            write_1d_complex(extinction_filename, extinction)

        extinction_ref = read_1d_complex(extinction_filename)

        if not skip_asserts:
            assert_allclose(extinction, extinction_ref, rtol=1e-3)

        if plot_extinction:
            # to plot the generated and reference solutions

            # calculate analytically
            extinction_analytical = sphere_extinction_analytical(freqs, radius)
            plt.figure(figsize=(8, 6))
            plt.plot(freqs * 1e-9, extinction.real)
            plt.plot(freqs * 1e-9, extinction_ref.real, "--")
            plt.plot(freqs * 1e-9, extinction_analytical, "x")
            plt.plot(freqs * 1e-9, extinction.imag)
            plt.plot(freqs * 1e-9, extinction_ref.imag, "--")
            plt.xlabel("f (GHz)")
            plt.legend(
                (
                    "Calculated (Re)",
                    "Reference (Re)",
                    "Analytical (Re)",
                    "Calculated (Im)",
                    "Reference (Im)",
                ),
                loc="right",
            )
            plt.title("Extinction with operator %s" % operator_name)
            plt.ylim(ymin=0)
            plt.show()


if __name__ == "__main__":
    test_extinction_all(plot_extinction=False, skip_asserts=True)
