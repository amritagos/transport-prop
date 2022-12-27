project('fastcpp', 'cpp',
 version : run_command('git', 'rev-parse',  '--short', 'HEAD').stdout().strip(),
  default_options : ['warning_level=3', 'cpp_std=c++20', 'debug=true', 'buildtype=debug'])

py = import('python').find_installation()
python_dep = py.dependency()
pyb11f = py.extension_module('fastcpp',
          sources : [
            'tcf.cpp',
            'binding.cpp'
                    ],
          dependencies : [dependency('fmt'),
                          dependency('python3'),
                          dependency('pybind11'),
                          dependency('eigen3'),
                         declare_dependency(link_args: '-lstdc++')],
          install : true)
