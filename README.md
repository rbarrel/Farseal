# Farseal

A Fortran library allowing for simulated annealing computation on enormous datasets by utilizing sparse matrix storage and computation. 

The name is a combination of the words "Fortran spARse Simulated annEALing". Based heavily on the work of Nicholas Herring ([OpenFSAM](https://github.com/nfherrin/OpenFSAM)).

# Installation

Requirements:

Description | Suggested Package | Optional?
--- | --- | ---
[A Fortran Compiler](https://gcc.gnu.org/fortran/) | GFortran (13.2.1) | No
[A BLAS Compliant Sparse Matrix Library](https://librsb.sourceforge.net/) | librsb (1.3.0) | No
[A Build System](https://mesonbuild.com/) | The Meson Build System (1.3.0) | Yes
[Backend for Build System](https://ninja-build.org/) | Ninja Build (1.11.1) | Yes
[A Fortran Testing Library](https://github.com/fortran-lang/test-drive) | test-drive (0.4.0) | Yes

Clone the repository then run `meson setup build` and `meson install -C build` to install Farseal as a library on your system.
