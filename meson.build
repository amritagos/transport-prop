project('transportProp', 'cpp', version:'0.1.1', default_options : ['warning_level=3', 'cpp_std=c++20', 'optimization=3'])

compiler = meson.get_compiler('cpp')

cpp_args = []
cpp_args += compiler.get_supported_arguments([
  '-Wno-unused-local-typedefs',  # Ignore unused local typedefs warnings
  '-Wno-array-bounds',           # Suppress out-of-bounds array access warnings
  '-ffast-math',                 # Enable faster, non-IEEE math calculations
  '-fno-finite-math-only',       # Allow Inf and NaN
  # These are based on recommendations from
  # https://youtu.be/_enXuIxuNV4?si=LvtMqPwJ6jYbDY66
  '-fno-semantic-interposition', # Assume no interposition for module functions
  '-fno-plt',                    # Avoid PLT for function calls within shared libs
  '-Bsymbolic',                  # Resolve symbols to local definitions
])

if cpp_args.length() > 0
  message('Adding compiler flags', cpp_args)
endif

inc = include_directories([['include', 'thirdparty']])

deps = []

fastcpp_dep = declare_dependency(include_directories : inc)
deps += fastcpp_dep

# Python module
py = import('python').find_installation('python', modules: ['numpy'])

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()
inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)

python_deps = deps + [
  py.dependency(),
  dependency('pybind11'),
  np_dep
]

# for the bindings
py.extension_module(
  'fastcpp',
  sources : [
    'python_bindings/bindings.cpp',
    'src/generic.cpp',
    'src/rdf.cpp',
    'src/tcf.cpp'
  ],
  # include_directories: ,
  dependencies: python_deps,
  cpp_args : cpp_args,
  install: true,
  subdir: 'transportProp' # path in site-packages
)

# pyflowy main package
py.install_sources([
    'python/transportProp/__init__.py',
    'python/transportProp/cli.py',
    'python/transportProp/io.py',
    'python/transportProp/misc.py',
    'python/transportProp/msd.py',
    'python/transportProp/structures.py',
    'python/transportProp/util.py',
    'python/transportProp/vacf.py',
  ],
  pure: false, # install next to compiled extension
  subdir: 'transportProp' # path in site-packages
)

