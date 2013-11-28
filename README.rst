pyfmmlib
========

pyfmmlib is a Python wrapper for `fmmlib2d
<https://cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html>`_ and `fmmlib3d
<https://cims.nyu.edu/cmcl/fmm3dlib/fmm3dlib.html>`_ implementations of the
`fast multipole method <https://en.wikipedia.org/wiki/Fast_multipole_method>`_ for
`Laplace <https://en.wikipedia.org/wiki/Laplace%27s_equation>`_ and
`Helmholtz <https://en.wikipedia.org/wiki/Helmholtz_equation>`_ potentials by
Zydrunas Gimbutas and Leslie Greengard (and including code by many more people).

This wrapper is far from comprehensive. It just catches the things I ended up
needing. Nonetheless, the FMMs and a fair bit of other useful stuff is accessible.

- Andreas Kloeckner <inform@tiker.net>

.. image:: https://badge.fury.io/py/pyfmmlib.png
    :target: http://pypi.python.org/pypi/pyfmmlib

Installation
------------

To build this, you need

* `numpy <http://numpy.org>`_
* `mako <http://makotemplates.org>`_ (`pip <https://pypi.python.org/pypi/pip>`_ install mako or `ez_install mako`)

Run::

    python setup.py install

as usual and cross your fingers.

If you'd like (rather effective) go-fast stripes,
on GNU platforms, use this::

    ./setup-optimized.sh install

Documentation
-------------

Not much, unfortunately. Here's what I do to figure out how to use stuff::

    >>> import pyfmmlib
    >>> dir(pyfmmlib)
    ['__builtins__', '__doc__', '__file__', '__name__', '__package__', '_add_plot', ...]

    Fish the desired function from this list (let's use 'legefder' as an
    example) and run:

    >>> print pyfmmlib.legefder.__doc__
    legefder - Function signature:
      val,der = legefder(x,pexp,[n])
    Required arguments:
      x : input float
      pexp : input rank-1 array('d') with bounds (n + 1)
    Optional arguments:
      n := (len(pexp)-1) input int
    Return objects:
      val : float
      der : float

This tells you how to call the function from Python.
You can then use grep to fish out the right Fortran source::

    $ grep -icl 'legefder' fmmlib*/*/*.f
    fmmlib3d/src/legeexps.f

Then look at the docs there, and you're in business. No idea what
function name to look for? Just use the same grep procedure to look
for keywords.

Crude, but effective. :)

Two more things:

* Some functions are wrapped with a ``_vec`` suffix. This means they
  apply to whole vectors of arguments at once. They're also parallel
  via OpenMP.

* ``pyfmmlib.fmm_part`` and ``pyfmmlib.fmm_tria`` are (dimension-independent)
  wrappers that make the calling sequence for the FMMs just a wee bit less
  obnoxious.  See ``examples/fmm.py`` for more.

  Here's a rough idea how these are called::

      from pyfmmlib import fmm_part, HelmholtzKernel

      pot, grad = fmm_part("PG", iprec=2, kernel=HelmholtzKernel(5),
              sources=sources, mop_charge=1, target=targets)

  Unlike the rest of the library (which calls directly into Fortran),
  these routines expect ``(n,3)``-shaped (that is, C-Order) arrays.

License
-------

This wrapper is licensed under the MIT license, as below. Beware, though, that
`fmmlib{2,3}d` are licensed under the GNU General Public License v2, which is
more restrictive.

Copyright (C) 2013 Andreas Kloeckner

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
