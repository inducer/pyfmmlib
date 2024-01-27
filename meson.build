project('pyfmmlib', 'c', 'fortran',
  license: 'MIT',
  meson_version: '>=1.0.0',
  default_options : ['warning_level=2', 'fortran_std=legacy'],
)

# {{{ gather includes

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

# NOTE: meson does not like absolute paths, so we dance around it
incdir_numpy = run_command(py,
  [
    '-c',
    '''import os
import numpy as np
try:
  incdir = os.path.relpath(np.get_include())
except Exception:
  incdir = np.get_include()
print(incdir)
  '''], check: true).stdout().strip()
inc_numpy = include_directories(incdir_numpy)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)

# }}}

# {{{ generate vec wrappers

generate_wrappers = files('vec_wrappers.py')
vec_wrappers = custom_target('vec_wrappers',
  output: 'vec_wrappers.f90',
  command: [py, generate_wrappers, '@OUTPUT@'])

# }}}

# {{{ generate f2py wrappers

f2py_name = '_internal'
f2py_module_c = f'@f2py_name@module.c'
f2py_wrapper_f = f'@f2py_name@-f2pywrappers.f'

generate_f2py = files('generate-f2py.py')
f2py_dep = custom_target('generate_f2py',
  input: 'wrappers.pyf.mako',
  output: [f2py_module_c, f2py_wrapper_f],
  command: [py, generate_f2py, '@INPUT@', '-o', '@OUTDIR@'])

# }}}


# {{{ fortranobject

f2py_object_c = 'fortranobject.c'
f2py_object_static = static_library('_fortranobject',
  incdir_f2py / f2py_object_c,
  dependencies: py_dep,
  include_directories: [inc_numpy, inc_f2py])

f2py_object_dep = declare_dependency(
  link_with: f2py_object_static,
  include_directories: [inc_numpy, inc_f2py])

# }}}

# {{{ fmmlib

fmmlib2d_sources = [
  'fmmlib2d/src/cdjseval2d.f',
  # 'fmmlib2d/src/d2mtreeplot.f',
  # 'fmmlib2d/src/d2tstrcr_omp.f',
  'fmmlib2d/src/hank103.f',
  'fmmlib2d/src/hfmm2dpart.f',
  'fmmlib2d/src/l2dterms.f',
  'fmmlib2d/src/lfmm2dpart.f',
  # 'fmmlib2d/src/second-r8.f',
  'fmmlib2d/src/cfmm2dpart.f',
  'fmmlib2d/src/d2tstrcr.f',
  'fmmlib2d/src/h2dterms.f',
  'fmmlib2d/src/helmrouts2d.f',
  'fmmlib2d/src/hfmm2drouts.f',
  'fmmlib2d/src/laprouts2d.f',
  'fmmlib2d/src/lfmm2drouts.f',
  'fmmlib2d/src/prini.f',
  'fmmlib2d/src/rfmm2dpart.f',
  'fmmlib2d/src/zfmm2dpart.f',
]

fmmlib2d_static = static_library('fmmlib2d_static', fmmlib2d_sources, pic: true)

fmmlib3d_sources = [
  'fmmlib3d/src/cdjseval3d.f',
  'fmmlib3d/src/dfft.f',
  'fmmlib3d/src/h3dtrirouts.f',
  'fmmlib3d/src/hfmm3drouts.f',
  'fmmlib3d/src/l3dtrans.f',
  'fmmlib3d/src/legeexps.f',
  'fmmlib3d/src/lfmm3dtria.f',
  'fmmlib3d/src/prinm.f',
  'fmmlib3d/src/rotviarecur3.f',
  'fmmlib3d/src/triahquad.f',
  'fmmlib3d/src/triquadflatlib.f',
  # 'fmmlib3d/src/d3mtreeplot.f',
  'fmmlib3d/src/h3dterms.f',
  'fmmlib3d/src/helmrouts3d.f',
  'fmmlib3d/src/hfmm3dtria.f',
  'fmmlib3d/src/l3dtrirouts.f',
  'fmmlib3d/src/lfmm3dpart.f',
  'fmmlib3d/src/projections.f',
  # 'fmmlib3d/src/second-r8.f',
  'fmmlib3d/src/triasymq.f',
  'fmmlib3d/src/yrecursion.f',
  'fmmlib3d/src/d3tstrcr.f',
  'fmmlib3d/src/h3dtrans.f',
  'fmmlib3d/src/hfmm3dpart.f',
  'fmmlib3d/src/l3dterms.f',
  'fmmlib3d/src/laprouts3d.f',
  'fmmlib3d/src/lfmm3drouts.f',
  'fmmlib3d/src/prini.f',
  'fmmlib3d/src/rotproj.f',
  'fmmlib3d/src/triagauc.f',
  'fmmlib3d/src/trilib.f',
]

fmmlib3d_static = static_library('fmmlib3d_static', fmmlib3d_sources, pic: true)

# }}}

# {{{ extension

py.extension_module(f2py_name,
  [f2py_dep, vec_wrappers],
  link_with: [fmmlib2d_static, fmmlib3d_static],
  subdir: 'pyfmmlib',
  link_language: 'fortran',
  dependencies: [f2py_object_dep],
  include_directories: [incdir_numpy, incdir_f2py],
  install: true)

py.install_sources([
  'pyfmmlib/__init__.py',
  'pyfmmlib/version.py',
  ],
  subdir: 'pyfmmlib')

# }}}
