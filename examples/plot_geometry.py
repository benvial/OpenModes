#!/usr/bin/env python
# coding: utf-8

"""
Working with the included and custom geometries
================================================

For convenience, a number of common meta-atom geometries are included with OpenModes. 
The code below shows a list of those which are currently available

"""

import openmodes
import os
import os.path as osp

os.listdir(openmodes.geometry_dir)


#################################################################################
# Geometric parameters
# ------------------------------------
# 
# Most of these geometries are parameterised, so that you can specify certain
# dimensions when you load the file. For example, when loading a cross, we can
# specify the width and height (in metres). If we don't specify them, then
# default values will be used, which in most cases result in meta-atoms of
# around 10mm in size.


sim = openmodes.Simulation()
cross_filename = osp.join(openmodes.geometry_dir, "cross.geo")
mesh = sim.load_mesh(cross_filename, parameters = {'width': 4e-3, 'height': 20e-3}, mesh_tol=0.5e-3)
sim.place_part(mesh)

#################################################################################
# We can plot the resulting geometry. Try tweaking the parameters in the call to `load_mesh`, and you should be able to see the changes in the plot below. If the browser-based plot doesn't work for you, try one of the other plot types described in [How to create 3D plots](How%20to%20create%203D%20plots.ipynb).


sim.plot_3d(style='wireframe')

#################################################################################
# Mesh tolerance
# ------------------------------------
# 
# An important additional option to play with is `mesh_tol`. This specifies how
# densely your object will be meshed. To see this in a 3D plot, pass the option
# `wireframe=True`, or click on the check-box.
# 
# The example above with the cross has a very dense mesh, which will
# unnecessarily slow down the computation. Note that the time taken for
# certain computations is $N^3$, where $N$ is approximately equal to the number
# of triangles in the mesh. We can see the number of triangles with the
# following command.


len(mesh.polygons)

#################################################################################
# In general, you should start working with a coarser mesh, to get some intial
# results. Once you have found the results that you want, re-run everything with
# increased mesh density until the results have converged sufficiently.



################################################################################# 
# Which geometry parameters can be modified?
# ------------------------------------
# 
# Unless you wrote the geometry file yourself, you may not know the names of the
# parameters which can be tweaked. However, this is fairly easy to see if you
# inspect the geometry file, which is written in plain text.


with open(cross_filename) as infile:
    print(infile.read())


#################################################################################
# In gmsh, it is possible to check if a parameter is defined with the
# `exists` command. So you can see that if `width` or `height` are not already
# defined, they will be given some default values. If they are specified in the
# `parameters` of `load_mesh`, then the given values will be used.
# 

#################################################################################
# Creating your own geometry
# ------------------------------------
# 
# Geometries for OpenModes can be created with the free CAD program [gmsh](http://geuz.org/gmsh). This is quite a powerful tool for creating geometries, although it takes some getting used to. It is a script-based system, but it does come with a GUI to show you the resulting geometry. You should use this GUI to check that your geometry is correct, and to debug the inevitable errors which will occur. Only once you can succesfully generate the 2D mesh from within gmsh should you attempt to use your geometry file with OpenModes.
# 
# A geometry can either be a zero-thickness sheet (often used for thin metal
# layers at microwave frequencies), or it may be a fully three-dimensional
# object. When creating a 3D object, it is important that the mesh for each part
# is water-tight. Otherwise you may get some very strange numerical results.
# 
# If you have created a geometry file which may be useful to other users, then
# you are encouraged to send it to me for inclusion in OpenModes (see the
# documentation for contact details).

