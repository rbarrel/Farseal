# Farseal

A Fortran library allowing for simulated annealing computation on enormous datasets by utilizing sparse matrix storage and computation. 

The library name is a combination of the words "**F**ortran sp**AR**se **S**imulated ann**EAL**ing". The code is based heavily on the work of Nicholas Herring with ([OpenFSAM](https://github.com/nfherrin/OpenFSAM)).

# Installation

Requirements:

Description | Suggested Package | Optional?
--- | --- | ---
A Modern Fortran Compiler (F2008 or later) | [GNU Fortran (13.2.1)](https://gcc.gnu.org/fortran/) | No
A BLAS Compliant Sparse Matrix Library | [librsb (1.3.0)](https://librsb.sourceforge.net/) | No
A Build System | [The Meson Build System (1.3.0)](https://mesonbuild.com/) | No
Backend for Build System | [Ninja Build (1.11.1)](https://ninja-build.org/) | No
A Fortran Testing Library | [test-drive (0.4.0)](https://github.com/fortran-lang/test-drive) | Yes

Clone the repository then run `meson setup build` and `meson install -C build` to install Farseal as a library on your system.

If you want to be able to run the test suite, install `test-drive` using `meson`.
