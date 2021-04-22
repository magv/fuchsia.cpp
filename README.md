# FUCHSIA

*Fuchsia* reduces differential equations for Feynman master
integrals to an epsilon form.

This is a new version of *Fuchsia*, written in C++ using GiNaC.
It is based on the previous (Python) version by O. Gituliar and
V. Magerya, as described in [2]. That version lives on at [1].

This *Fuchsia* is still under development, with the main missing
feature being the support for unfactorizable polynomials of powers
higher than one in the denominators. Otherwise this version is
faster and easies to use than the old one.

When compiled, *Fuchsia* becomes a single executable. A precompiled
and statically linked version of this executable can be found
in the releases section on Github. See the manual below for its
usage. Also see the manual of the previous version at [2], and the
articles by R. N. Lee and A. Pomeransky [3-4] for a discussion
of the algorithms used.

* [1] https://github.com/gituliar/fuchsia/
* [2] https://arxiv.org/abs/1701.04269
* [3] https://arxiv.org/abs/1411.0911
* [4] https://arxiv.org/abs/1707.07856

# MANUAL

## NAME

`fuchsia` -- transform linear differential equations into an epsilon
form.

## SYNOPSYS

`fuchsia` [options] **command** *args* ...

## DESCRIPTION

`fuchsia` transforms systems of linear differential equations,
>***∂/∂x I(x,ε) = M(x,ε) I(x,ε),***

into an epsilon form,
>***∂/∂x J(x,ε) = ε S(x) J(x,ε),***

where ***I*** and ***J*** are column vectors of functions in the original
and epsilon basis, ***M*** is the original matrix, ***ε×S*** is the
matrix in an epsilon form, and ***I*** is related to ***J*** via the
transformation matrix ***T(x,ε)*** such that
>***I = T J.***

In all cases ***M*** can depend on additional symbolic variables, which
are treated as independent constants.

## EXAMPLES

To reduce a single-variable differential system of equations to an
epsilon form, use this:

`fuchsia` **reduce** -x *x* -e *eps* *matrix.orig* -m *matrix.ep* -t *matrix.ep.t* \\\
    -C 2>&1 | tee *matrix.ep.log*

For differential equations in multiple variables this is the usage:

`fuchsia` **reduce** *matrix.x* *matrix.y* -x *x* -x *y* -e *eps* \\\
    -m *matrix.ep.x* -m *matrix.ep.y* -t *matrix.ep.t* \\\
    -C 2>&1 | tee *matrix.ep.log*

## COMMANDS

* **show** [-x *name*] *matrix*

  Show a description of a given matrix.

* **reduce** [-x *name*] [-e *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Find an epsilon form of the given matrix. Internally
  this is a combination of **reduce-diagonal-blocks**,
  **fuchsify-off-diagonal-blocks** and **factorize**.

* **reduce** [-x *name*] ... [-e *name*] [-m *path*] ... [-t *path*] [-i *path*] *matrix* ...

  Find an epsilon form of a given multivariate differential equation
  system. A matching number of *matrix* arguments, -x, and
  -m flags is required.

  The matrices are reduced one by one, and a single transformation
  is computed that simultaneously transforms all of them into an
  epsilon form. It may be best to list the simplest matrix first.

  NOTE: this command is under development, and may fail when it
  shouldn't.

* **reduce-diagonal-blocks** [-x *name*] [-e *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Transform the matrix into a block-triangular form and reduce the
  diagonal blocks into an epsilon form.

* **fuchsify-off-diagonal-blocks** [-x *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Transform the off-diagonal blocks of a block-triangular matrix
  into a Fuchsian form, assuming the diagonal blocks are already in
  an epsilon form, thus making the whole matrix normalized Fuchsian.

* **factorize** [-x *name*] [-e *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Find a transformation that will make a given normalized Fuchsian
  matrix proportional to the infinitesimal parameter.

* **fuchsify** [-x *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Find a transformation that will transform a given matrix into a
  Fuchsian form. This is less efficient than block-based commands,
  because it effectively treats the whole matrix as one big block.

* **normalize** [-x *name*] [-e *name*] [-m *path*] [-t *path*] [-i *path*] *matrix*

  Find a transformation that will transform a given Fuchsian matrix
  into a normalized form. This is less efficient than block-based
  commands, because it effectively treats the whole matrix as one
  big block.

* **sort** [-m *path*] [-t *path*] [-i *path*] *matrix*

  Find a block-triangular form of the given matrix by shuffling.

* **transform** [-x *name*] [-m *path*] *matrix* *transform* ...

  Transform a given matrix using a given transformation.

* **changevar** [-x *name*] [-y *name*] [-m *path*] *matrix* *expr*

  Perform a change of variable from x to y, such that x=expr(y).

* **suggest-changevar** [-x *name*] [-y *name*] *matrix*

  Suggest a rational change of variable that will transform residue
  eigenvalues of the form n/2+k×eps into n+k×eps, thus making it
  possible to find an epsilon form of the matrix.

  Note that some bad eigenvalues disappear when the matrix is
  fuchsified, so this routine is best used after "fuchsia fuchsify".

* **simplify** [-x *name*] ... [-m *path*] ... [-t *path*] [-i *path*] *matrix* ...

  Try to find a transformation that makes a given matrix (or a set
  of matrices) simpler, for some definition of "simple".

## OPTIONS

* -x *name*

  Use this name for the free variable (default: x).

* -0 *expr*

  Set this value for x during multivariate reduction (default: 0).

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

  Force colored output even if stdout is not a tty.

* -P

  Paranoid mode: spend more time checking internal invariants.

* -q

  Print a more quiet log.

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
