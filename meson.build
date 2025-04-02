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
    'python_bindings/fastcpp.cpp',
    'src/generic.cpp',
    'src/rdf.cpp',
    'src/tcf.cpp'
  ],
  # include_directories: ,
  dependencies: python_deps,
  cpp_args : cpp_args,
  install: true,
  subdir: 'transportProp/'
)

# pyflowy main package
py.install_sources([
    'transportProp/__init__.py',
    'transportProp/cli.py',
    'transportProp/io.py',
    'transportProp/misc.py',
    'transportProp/msd.py',
    'transportProp/structures.py',
    'transportProp/util.py',
    'transportProp/vacf.py',
  ],
  pure: false, # install next to compiled extension
  subdir: 'transportProp/'
)

