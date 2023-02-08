project('fastcpp', 'cpp',
 version : run_command('git', 'rev-parse',  '--short', 'HEAD').stdout().strip(),
  default_options : ['warning_level=3', 'cpp_std=c++20', 'debug=true', 'buildtype=debug'])
  # default_options : ['warning_level=3', 'cpp_std=c++20', 'debug=false', 'buildtype=release'])

py = import('python').find_installation()
python_dep = py.dependency()
cpplibs = shared_library('cppfuncs',
          sources : [
            'tcf.cpp',
            'generic.cpp',
            'rdf.cpp',
                    ],
          dependencies : [dependency('fmt'),
                          dependency('eigen3'),
                         declare_dependency(link_args: '-lstdc++')],
          install : true)
pyb11f = py.extension_module('fastcpp',
          sources : [
            'binding.cpp'
                    ],
          dependencies : [dependency('fmt'),
                          dependency('python3'),
                          dependency('pybind11'),
                          dependency('eigen3'),
                         declare_dependency(link_args: '-lstdc++')],
                         link_with: [cpplibs],
          install : true)
