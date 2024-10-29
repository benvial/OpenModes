#!/usr/bin/env python
# coding: utf-8


"""
Calculating Extinction
=======================

A simple example showing how to use OpenModes. A pair of split-ring resonators are 
created, and their extinction cross section is calculated.

"""


####################################################################################
# setup 2D and 3D plotting
from openmodes.ipython import matplotlib_defaults

matplotlib_defaults()
import matplotlib.pyplot as plt

# the numpy library contains useful mathematical functions
import numpy as np

# import useful python libraries
import os.path as osp

# import the openmodes packages
import openmodes
from openmodes.sources import PlaneWaveSource


####################################################################################
# Then we create a `Simulation` object which holds all the data for the simulation
# we are going to run. Since we are running in the IPython/Jupyter notebook, we pass
# `notebook=True` to enable in browser 3D plots and progress sliders.


sim = openmodes.Simulation(notebook=False)


####################################################################################
# Next we load the geometry file. This geometry file is written in the scripting
# language of [gmsh](http://geuz.org/gmsh), a program which converts the geometry
# into a surface mesh. The geometry file used is a split ring resonator provided
# with OpenModes. The installation location of these provided geometries is found
# in `openmodes.geometry_dir`. Also, these included geometries are written so that
# some of their geometric parameters can be modified. In this case the inner radius
# is set to 2.5mm and the outer radius to 4mm.
#
# The mesh density may be specified in the geometry file, but can be over-ridden
# with the parameter `mesh_tol`. Be careful with this parameter, setting it too
# small can result in very long computations.


filename = osp.join(openmodes.geometry_dir, "SRR.geo")
mesh_tol = 1e-3
outer_radius = 4e-3
srr = sim.load_mesh(
    filename,
    mesh_tol,
    parameters={"inner_radius": 2.5e-3, "outer_radius": outer_radius},
)


####################################################################################
# Now we place the parts within the simulation. By default the parts are placed at
# the origin, so after placement we need to move one of them to the desired distance
# from the other. Before moving the second ring, it is rotated $180^\circ$ to create
# a broadside-coupled configuration.


ring1 = sim.place_part(srr)
ring2 = sim.place_part(srr)
ring2.rotate(axis=[0, 0, 1], angle=180)
ring2.translate([0, 0, 2e-3])


####################################################################################
# To check that the parts have been placed in the correct location, we can visualise
# them using the provided `plot_3d` function. This creates a 3d view of the objects
# within the browser. You can control the view point with the mouse
# - Hold the left button to rotate the view
# - Scroll the wheel to zoon
# - Hold the right button to pan


sim.plot_3d(style='wireframe')



####################################################################################
# We will calculate the extinction cross section for this pair of SRRs, assuming
# that they are perfectly conducting. We want to excite it with a $y$ polarised
# plane wave, propagating in the $x$ direction.
#
# $$\mathbf{E}_{inc} = \hat{\mathbf{y}}\exp\left(-\frac{s}{c}\hat{\mathbf{x}}\cdot\mathbf{r}\right)$$
#
# We will calculate at 401 frequencies between 5 and 10 GHz.


####################################################################################
# the frequency range over which to calculate
freqs = np.linspace(5e9, 10e9, 401)

####################################################################################
# construct the source plane wave with given polarisation and propagation direction
e_inc = np.array([0, 1, 0])
k_dir = np.array([1, 0, 0])

####################################################################################
# The optional p_inc parameter ensures that the incident field is normalised to 1W/m^2
plane_wave = PlaneWaveSource(e_inc, k_dir, p_inc=1.0)


####################################################################################
# create an empty array to hold extinction data
extinction_single = np.empty(len(freqs), np.float64)
extinction_pair = np.empty(len(freqs), np.float64)


####################################################################################
# Now loop through all the frequencies. For notational convenience, the time
# dependence is assumed as $\exp(s t)$ with complex frequency $s$.
#
# At each frequency the impedance matrix $Z$ is calculated, as is the source term
# $V$ due to the plane wave. By default these matrices and vectors are composite
# objects accounting for all the parts. We can select one of these matrices to find
# the response of an isolated part, or we can combine them together to find the
# response of the entire system.
#
# In both cases, the extinction cross-section $\sigma_{ext}$ is found from the
# impedance matrix $Z$ and driving term $V$ as $\sigma_{ext}(s) = V^{*}(s)\cdot Z(s)
# \cdot V(s)$


for freq_count, s in sim.iter_freqs(freqs):
    Z = sim.impedance(s)
    V = sim.source_vector(plane_wave, s)

    Z_single = Z[ring1, ring1]
    V_single = V[:, ring1]

    # calculate the extinction only of one ring
    extinction_single[freq_count] = np.vdot(V_single, Z_single.solve(V_single)).real

    # calculate the extinction of the system of two rings
    extinction_pair[freq_count] = np.vdot(V, Z.solve(V)).real

####################################################################################
# Now we plot the extinction cross-section as a function of frequency.
# We use the cross section to a disc with the same radius as the SRR as a
# reference to express it in the normalised form of extinction efficiency $Q_{ext}$.


area = np.pi * outer_radius**2

plt.figure()
plt.plot(freqs * 1e-9, extinction_single / area, label="single")
plt.plot(freqs * 1e-9, extinction_pair / area, label="pair")
plt.legend(loc="upper right")
plt.axis("tight")
plt.xlabel("freq (GHz)")
plt.ylabel("$Q_{ext}$")
plt.title("Extinction efficiency")
plt.show()


####################################################################################
# This figure shows the fundamental resonance of a single ring. By bringing the two
# rings near to each other, this resonance splits into two coupled modes. Subsequent
# examples will show how such modes and be explicitly modelled, and serve as the
# basis for convenient semi-analytical models of these structures.
