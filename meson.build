project('openmodes','fortran','c',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)


py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()



#########
# Compiler-specific options (for OpenMP and other optimizations)
# Fortran compiler options
openmp_flags = ['-fopenmp']
c_openmp_flags = ['-fopenmp']
fortran_opt_flags = ['-g', '-O3', '-fimplicit-none']
# intel_fortran_flags = ['-openmp', '-O2', '/fpe:0', '/fp:fast=2']
omp = dependency('openmp')
dependencies = [py_dep,omp]

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)


core_source = custom_target('coremodule.c',
  input : ['src/common.f90','src/rwg.f90','src/.f2py_f2cmap'],
  output : ['coremodule.c', 'core-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT0@', '@INPUT1@', '-m', 'core', '--lower',
  '--f2cmap', '@INPUT2@', 
  '--overwrite-signature', 'only: set_threads get_threads face_integrals_hanninen z_efie_faces_self z_efie_faces_mutual arcioni_singular z_mfie_faces_self z_mfie_faces_mutual face_integrals_yla_oijala :']
)

py.extension_module('core',
  [
    'src/common.f90',
    'src/rwg.f90',
    core_source,
  ],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : dependencies,
  fortran_args: fortran_opt_flags + openmp_flags,
  c_args: c_openmp_flags,
  install : true,
  subdir: 'openmodes'
)


dunavant_source = custom_target('dunavantmodule.c',
  input : ['src/dunavant.f90','src/.f2py_f2cmap'],
  output : ['dunavantmodule.c', 'dunavant-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT0@', '-m', 'dunavant', '--lower', 
  '--f2cmap', '@INPUT1@', 
  '--overwrite-signature', 'only: dunavant_order_num dunavant_rule :']
)

py.extension_module('dunavant',
  [
    'src/dunavant.f90',
    dunavant_source,
  ],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : dependencies,
  fortran_args: fortran_opt_flags + openmp_flags,
  c_args: c_openmp_flags,
  install : true,
  subdir: 'openmodes'
)

python_sources = ['openmodes/integration.py', 'openmodes/__init__.py', 'openmodes/sources.py', 'openmodes/ipython.py', 'openmodes/model.py', 'openmodes/impedance.py', 'openmodes/helpers.py', 'openmodes/material.py', 'openmodes/array.py', 'openmodes/modes.py', 'openmodes/simulation.py', 'openmodes/parts.py', 'openmodes/visualise.py', 'openmodes/eig.py', 'openmodes/constants.py', 'openmodes/basis.py', 'openmodes/version.py', 'openmodes/multipole.py']

py.install_sources(
  python_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'openmodes'  # Folder relative to site-packages to install to
)

python_sources = ['openmodes/operator/__init__.py', 'openmodes/operator/operator.py', 'openmodes/operator/rwg.py', 'openmodes/operator/penetrable.py', 'openmodes/operator/singularities.py', 'openmodes/operator/pec.py']

py.install_sources(
  python_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'openmodes/operator'  # Folder relative to site-packages to install to
)

python_sources = ['openmodes/external/__init__.py', 'openmodes/external/ordered_set.py', 'openmodes/external/point_in_polygon.py']

py.install_sources(
  python_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'openmodes/external'  # Folder relative to site-packages to install to
)


python_sources = ['openmodes/templates/three_js_plot.html']
py.install_sources(
  python_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'openmodes/templates'  # Folder relative to site-packages to install to
)

python_sources = ['openmodes/static/three_js_plot.js']
py.install_sources(
  python_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'openmodes/static'  # Folder relative to site-packages to install to
)


python_sources = ['openmodes/mesh/__init__.py', 'openmodes/mesh/mesh.py', 'openmodes/mesh/freecad.py', 'openmodes/mesh/gmsh.py']
py.install_sources(
  python_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'openmodes/mesh'  # Folder relative to site-packages to install to
)

python_sources = ['openmodes/geometry/cross.geo', 'openmodes/geometry/canonical_spiral.geo', 'openmodes/geometry/SRR.geo', 'openmodes/geometry/cylinder_hole.geo', 'openmodes/geometry/sphere.geo', 'openmodes/geometry/single.geo', 'openmodes/geometry/single-bent.geo', 'openmodes/geometry/circle.geo', 'openmodes/geometry/solid_L.geo', 'openmodes/geometry/rectangle.geo', 'openmodes/geometry/horseshoe_rect.geo', 'openmodes/geometry/closed_ring.geo', 'openmodes/geometry/isosceles.geo', 'openmodes/geometry/torus.geo', 'openmodes/geometry/box.geo', 'openmodes/geometry/cylinder_rounded.geo', 'openmodes/geometry/ellipsoid.geo', 'openmodes/geometry/asymmetric_ring.geo', 'openmodes/geometry/cylinder_hollow.geo', 'openmodes/geometry/v_antenna.geo', 'openmodes/geometry/circled_cross.geo', 'openmodes/geometry/elliptical_cylinder.geo']

py.install_sources(
  python_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'openmodes/geometry'  # Folder relative to site-packages to install to
)
