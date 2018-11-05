# FUCHSIA

*Fuchsia* reduces differential equations for Feynman master
integrals to an epsilon form.

This is a new version of *Fuchsia*, written in C++ using GiNaC.
This version is still under development, with the main missing
feature being the support for unfactorizable polynomials of powers
higher than one in the denominators. Otherwise this version is
faster and easies to use than the old one.

The previous (Python) version can be found at [1].

When compiled, *Fuchsia* becomes a single executable. See the manual
below for its usage. Also see the manual for the previous version
at [2] for a discussion of the algorithms used. A precompiled
and statically linked version of this executable can be found
in the releases section on Github.

* [1] https://github.com/gituliar/fuchsia/
* [2] https://arxiv.org/abs/1701.04269

# MANUAL

## NAME

`fuchsia` -- transform linear differential equations into epsilon form.

## SYNOPSYS

`fuchsia` [options] **command** *args* ...

## COMMANDS

* **show** [-x *name*] *matrix*

  Show a description of a given matrix.

* **reduce** [-x *name*] [-e *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Find an epsilon form of the given matrix. Internally
  this is a combination of **reduce-diagonal-blocks**,
  **fuchsify-off-diagonal-blocks** and **factorize**.

* **reduce-diagonal-blocks** [-x *name*] [-e *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Transform the matrix into block-triangular form and reduce the diagonal
  blocks into epsilon form.

* **fuchsify-off-diagonal-blocks** [-x *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Transform the off-diagonal blocks of a block-triangular matrix into
  Fuchsian form, assuming the diagonal blocks are already in epsilon
  form, thus making the whole matrix normalized Fuchsian.

* **factorize** [-x *name*] [-e *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Find a transformation that will make a given normalized Fuchsian matrix
  proportional to the infinitesimal parameter.

* **fuchsify** [-x *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Find a transformation that will transform a given matrix into Fuchsian
  form. This is less efficient than block-based commands, because it
  effectively treats the whole matrix as one big block.

* **normalize** [-x *name*] [-e *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Find a transformation that will transform a given Fuchsian matrix into
  normalized form. This is less efficient than block-based commands,
  because it effectively treats the whole matrix as one big block.

* **sort** [-m *path*] [-t *path*] [-i *path*] *matrix*

  Find a block-triangular form of the given matrix by shuffling.

* **transform** [-x *name*] [-m *path*] *matrix* *transform* ...

  Transform a given matrix using a given transformation.

* **changevar** [-x *name*] [-y *name*] [-m *path*] *matrix* *expr*

  Perform a change of variable from x to y, such that x=expr(y).

* **suggest-changevar** [-x *name*] [-y *name*] *matrix*

  Suggest a rational change of variable that will transform residue
  eigenvalues of the form n/2+k×eps into n+k×eps, thus making it possible
  to find an epsilon form of the matrix.

## OPTIONS

* -x *name*

  Use this name for the free variable (default: x).

* -y *name*

  Use this name for the new free variable (default: y).

* -e *name*

  Use this name for the infinitesimal parameter (default: eps).

* -m *path*

  Save the resulting matrix into this file.

* -t *path*

  Save the resulting transformation into this file.

* -i *path*

  Save the inverse transformation into this file.

* -C

  Force colored output even stdout is not a tty.

* -P

  Paranoid mode: spend more time checking internal invariants.

* -v

  Print a more verbose log.

* -h

  Show this help message.

* -V

  Print version information.

## ARGUMENTS

* *matrix*

  Read the input matrix from this file.

* *transform*

  Read the transformation matrix from this file.

* *expr*

  Arbitrary expression.

## AUTHORS

Vitaly Magerya
