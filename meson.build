project(
  'Farseal',
  'fortran',
  version: '0.1',
  default_options: ['werror=true'],
)

install = not (meson.is_subproject() and get_option('default_library') == 'static')

fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

if fc_id == 'gcc'
  # Silence Nuisance GFortran Warnings
  compiler_args = [
    '-Wno-unused-dummy-argument',
    '-Wno-uninitialized',
    '-Wfatal-errors',
  ]

  foreach arg: compiler_args
    add_project_arguments(arg, language: 'fortran')
  endforeach
endif

librsb_dep = dependency('librsb', required: true)

srcs = []
subdir('src')
subdir('test')

# Build/Install Farseal
farseal_lib = library(
  meson.project_name(),
  sources: srcs,
  install: install,
)

pkg = import('pkgconfig')
pkg.generate(
  farseal_lib,
  description: 'Fortran spARse Simulated annEALing library',
)

# Export Farseal as a Dependency for Other Projects
farseal_inc = farseal_lib.private_dir_include()
farseal_dep = declare_dependency(
  link_with: farseal_lib,
  include_directories: farseal_inc,
)
