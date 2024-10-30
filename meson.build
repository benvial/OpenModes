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
c_flags = ['-fopenmp', '-Wno-cast-function-type','-Wno-unused-parameter', 
'-Wno-unused-variable','-Wno-missing-field-initializers']
fortran_opt_flags = ['-g', '-O3', '-fimplicit-none','-Wno-unused-variable']
# intel_fortran_flags = ['-openmp', '-O2', '/fpe:0', '/fp:fast=2']
omp = dependency('openmp', required: false)
dependencies = [py_dep, omp]

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)


subdir('openmodes')