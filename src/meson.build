
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

fortran_sources = [
  'r8_normal.f',
  'r8_uniform.f',
  'rk1_ti_solve.f',
  'rk1_tv_solve.f',
  'rk2_ti_solve.f',
  'rk2_tv_solve.f',
  'rk3_ti_solve.f',
  'rk4_ti_solve.f',
  'rk4_tv_solve.f',
  'output.f'
]

pyf_sources = [
  'sdepack.pyf'
]

message('Fortran sources: ', fortran_sources + pyf_sources)

f2py_target = custom_target(
  py_mod_name + '_f2py',
  input: pyf_sources,
  output: [py_mod_name + 'module.c', py_mod_name + '-f2pywrappers.f'],
  command: [f2py, '@INPUT@', '--lower', '--build-dir', meson.current_build_dir()]
)

if host_machine.system() == 'windows'
  fortran_link_args = ['-static-libgfortran', '-static-libgcc', '-static-libquadmath']
else
  fortran_link_args = []
endif

py.extension_module(
  py_mod_name,
  fortran_sources + [f2py_target, fortranobject_c],
  c_args: ['-I' + incdir_numpy, '-I' + incdir_f2py],
  link_args: fortran_link_args,
  dependencies: py_dep,
  install: true
)
