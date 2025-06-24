py_subdir = 'diskwarp'
python_sources = [
  '__init__.py',
  'helper.py',
  ]
py3.install_sources(python_sources, subdir: py_subdir)

fortran_module_name = '_fortran'
fortran_cwrapper = custom_target(fortran_module_name,
  input: ['fortran.f90'],
  output: [fortran_module_name+'module.c', fortran_module_name + '-f2pywrappers2.f90'],
  command: [py3, '-m', 'numpy.f2py', '--lower', '-m', fortran_module_name, 
  '@INPUT@', '--build-dir', py_subdir]
)
py3.extension_module(fortran_module_name,
  'fortran.f90',
  fortran_cwrapper,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
  subdir: py_subdir
)