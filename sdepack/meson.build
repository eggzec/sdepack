
add_languages('fortran')

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

f2py = find_program('f2py', required: true)

incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true
).stdout().strip()

incdir_f2py = run_command(py,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check: true
).stdout().strip()

fortranobject_c = run_command(py,
  ['-c', 'import numpy.f2py; import os; print(os.path.join(numpy.f2py.get_include(), "fortranobject.c"))'],
  check: true
).stdout().strip()

py_mod_name = 'sdepack'

fortran_sources = run_command(py, ['-c',
    'import glob; print(" ".join(glob.glob("*.f")))'
], check: true).stdout().strip().split()
message('Fortran sources: ', fortran_sources)

f2py_target = custom_target(
  py_mod_name + '_f2py',
  input: fortran_sources,
  output: ['_' + py_mod_name + 'module.c', '_' + py_mod_name + '-f2pywrappers.f'],
  command: [
    f2py,
    '@INPUT@',
    '-m', '_' + py_mod_name,
    '--lower',
    '--build-dir', meson.current_build_dir()
  ]
)

py.extension_module(
  '_' + py_mod_name,
  fortran_sources + [f2py_target, fortranobject_c],
  c_args: ['-I' + incdir_numpy, '-I' + incdir_f2py],
  dependencies: py_dep,
  install: true
)

python_sources = run_command(py, ['-c',
    'import glob; print(" ".join(glob.glob("*.py")))'
], check: true).stdout().strip().split()
message('Python sources: ', python_sources)

py.install_sources(
  python_sources,
  subdir: py_mod_name
)
