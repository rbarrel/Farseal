project(
  'Farseal',
  'fortran',
  version: '0.1',
)

fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

if fc_id == 'gcc'
  # Silence Nuisance GFortran Warnings
  compiler_args = [
    '-Wno-unused-dummy-argument',
    '-Wno-uninitialized',
  ]

  foreach arg: compiler_args
    add_project_arguments(arg, language: 'fortran')
  endforeach
endif

librsb_dep = dependency('librsb', required: true)

srcs = []
subdir('src')
subdir('test')
