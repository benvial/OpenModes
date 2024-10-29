

core_source = custom_target('coremodule.c',
  input : ['common.f90','rwg.f90','.f2py_f2cmap'],
  output : ['coremodule.c', 'core-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT0@', '@INPUT1@', '-m', 'core', '--lower',
  '--f2cmap', '@INPUT2@', '--build-dir', 'openmodes/src', '--quiet',
  '--overwrite-signature', 'only: set_threads get_threads face_integrals_hanninen z_efie_faces_self z_efie_faces_mutual arcioni_singular z_mfie_faces_self z_mfie_faces_mutual face_integrals_yla_oijala :']
)

py.extension_module('core',
  [
    'common.f90',
    'rwg.f90',
    core_source,
  ],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : dependencies,
  fortran_args: fortran_opt_flags + openmp_flags,
  c_args: c_flags,
  install : true,
  subdir: 'openmodes'
)

dunavant_source = custom_target('dunavantmodule.c',
  input : ['dunavant.f90','.f2py_f2cmap'],
  output : ['dunavantmodule.c', 'dunavant-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT0@', '-m', 'dunavant', '--lower', 
  '--f2cmap', '@INPUT1@',  '--build-dir', 'openmodes/src', '--quiet',
  '--overwrite-signature', 'only: dunavant_order_num dunavant_rule :']
)

py.extension_module('dunavant',
  [
    'dunavant.f90',
    dunavant_source,
  ],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : dependencies,
  fortran_args: fortran_opt_flags + openmp_flags,
  c_args: c_flags,
  install : true,
  subdir: 'openmodes'
)
