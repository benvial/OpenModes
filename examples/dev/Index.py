#!/usr/bin/env python
# coding: utf-8

# # OpenModes examples
# 
# Example notebooks for the [OpenModes](http://davidpowell.github.io/OpenModes/) mode
# solver for open electromagnetic resonators such as metamaterials.
# 
These examples are written for the latest version of OpenModes, and may not work with older
versions. Please see the [updating instructions](http://openmodes.readthedocs.org/en/latest/install.html#upgrading-to-a-newer-version)
in the documentation.

## Performing calculations with OpenModes
These examples deal with the simplest case of metal structures for microwave frequencies which can be approximated as perfect electric conductors (PEC).
* [Calculating Extinction Cross-Sections](Calculating%20Extinction.ipynb) - This is the simplest calculation, which shows how a split ring interacts with the incident field, and providing a reference solution for more advanced models.
* [Visualising Modes](Visualising%20Modes.ipynb) - In this example, the modes of a split ring are found, giving their (complex) resonant frequencies, and the currents and charges of the modes.
* [Constructing Scalar Models](Constructing%20Scalar%20Models.ipynb) - The modes are used to construct simple models, which model the split rings quite accurately over a broad bandwidth.
* [Modelling Coupled Elements](Modelling%20Coupled%20Elements.ipynb) - A model is constructed for a pair of split rings which are coupled by their near fields. 

## Working with dielectric structures
* [Dielectric disc](Dielectric%20disc.ipynb) - Solve the modes of a silicon nano-disk, operating in the optical frequency range. Utilises a surface equivalent problem to represent the fields within objects.

## Understanding the features of OpenModes
* [How to create 3D plots](How%20to%20create%203D%20plots.ipynb) - This example shows all the different ways you can view the 3D solution on an element, including interactive in-browser plots
* [Using and creating geometric shapes](Using%20and%20creating%20geometric%20shapes.ipynb) - OpenModes is quite general, working not just with split-rings, but with almost any particle shape. This shows how to use different geometries, and to create your own. 

## Running the examples

If you are viewing this page through the website http://nbviewer.ipython.org, then you can see the results, but cannot run any calculations. To run the examples yourself and play with them, please follow the instructions in the [documentation](http://openmodes.readthedocs.org/en/latest/examples.html#running-examples).


# In[ ]:




