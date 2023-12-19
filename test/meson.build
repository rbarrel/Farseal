testdrive_dep = dependency('test-drive', required: false)

if testdrive_dep.found()
  tests = [
    'CoolingTest',
  ]

  test_srcs = files(
    'main.f90'
  )

  foreach t: tests
    test_srcs += files('@0@.f90'.format(t.underscorify()))
  endforeach

  farseal_test = executable(
    'farseal_test',
    sources: [srcs, test_srcs],
    dependencies: [librsb_dep, testdrive_dep],
  )

  foreach t: tests
    test('Test @0@'.format(t.to_upper()), farseal_test, args: t)
  endforeach
endif
