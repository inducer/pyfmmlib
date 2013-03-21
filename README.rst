pyfmmlib
========

pyfmmlib is a Python wrapper for `fmmlib2d
<https://cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html>`_ and `fmmlib3d
<https://cims.nyu.edu/cmcl/fmm3dlib/fmm3dlib.html>`_ implementations of the
fast multipole method for Laplace and Helmholtz potentials by Zydrunas Gimbutas
and Leslie Greengard (and including code by many more people).

This wrapper is far from comprehensive. It just catches the things I ended up
needing. Nonetheless, the FMMs and a fair bit of other useful stuff is accessible.

- Andreas Kloeckner <inform@tiker.net>

Installation
------------

To build this, you need

* `numpy <http://numpy.org>`_
* `mako <http://makotemplates.org>`_ (`pip install mako`)

Type::

    ./grab-sources.sh

This will download and unpack the fmmlib sources somewhere below here where the
build expects them to be.

Then run::

    python setup.py install

as usual and cross your fingers.

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
You can then use grep to fish out the right Fortran source, look at the docs
there, and you're in business. Crude, but effective.

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
