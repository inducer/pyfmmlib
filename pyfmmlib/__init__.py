from __future__ import division
import pyfmmlib._internal as _int
from pyfmmlib._internal import *  # noqa
from pyfmmlib.version import VERSION_TEXT as __version__  # noqa

import numpy as np


# Map from matrix indices (i,j) in Hessian into output array of
# FMM 'hess' return array.

hessian_index_lookup = {
        3: {
            # diagonal
            (0, 0): 0,
            (1, 1): 1,
            (2, 2): 2,

            # off-diagonal
            (1, 0): 3,
            (0, 1): 3,
            (2, 0): 4,
            (0, 2): 4,
            (2, 1): 5,
            (1, 2): 5,
        },

        2: {
            (0, 0): 0,
            (1, 0): 1,
            (0, 1): 1,
            (1, 1): 2,
            }
        }


# {{{ kernel classes

class KernelBase(object):
    flag_value = 1

    def __hash__(self):
        return hash((self.__class__,) + (self.__getinitargs__()))

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
                self.__getinitargs__() == other.__getinitargs__()


class LaplaceKernel(KernelBase):
    letter = "l"

    def __getinitargs__(self):
        return ()

    def __str__(self):
        return "lap"

    @property
    def k_arg(self):
        return ()

    def k_kwargs(self):
        return {}

    def evaluate(self, evaluate_func):
        return self


class HelmholtzKernel(KernelBase):
    letter = "h"

    def __init__(self, k):
        self.k = k

    def __getinitargs__(self):
        return (self.k, )

    def __str__(self):
        return "helm(%s)" % self.k

    def evaluate(self, evaluate_func):
        return type(self)(evaluate_func(self.k))

    @property
    def k_arg(self):
        return (self.k, )

    def k_kwargs(self, name="zk"):
        return {name: self.k}


class DifferenceKernel(KernelBase):
    letter = "h"
    flag_value = 2

    def __init__(self, k):
        self.k = k

    def __getinitargs__(self):
        return (self.k, )

    def __str__(self):
        return "diff(%s)" % self.k

    def evaluate(self, evaluate_func):
        return type(self)(evaluate_func(self.k))

    @property
    def k_arg(self):
        return (self.k, )

    def k_kwargs(self, name="zk"):
        return {name: self.k}


def normalize_kernel(kernel):
    if not isinstance(kernel, KernelBase):
        if kernel == 0:
            kernel = LaplaceKernel()
        else:
            kernel = HelmholtzKernel(kernel)

    return kernel

# }}}


def _fmm(dimensions, size, kind, source_args, what, iprec, kernel,
        slp_density=None, dlp_density=None,
        target=None, dipvec=None, force_hess=False,
        mesh=None, debug=False):
    """Internal routine, do not call directly."""

    if dimensions == 3:
        kernel_scale = 1/(4*np.pi)
        grad_scale = -1
    elif dimensions == 2:
        kernel_scale = 1
        grad_scale = 1
    else:
        raise RuntimeError("unsupported dimensionality")

    kernel = normalize_kernel(kernel)

    if not what:
        return ()

    if target is None and any(l.isupper() for l in what):
        raise RuntimeError("must specify target if on-target computation desired")

    if target is None:
        ntarget = 0
        nalloctarget = 1

        target = np.empty((nalloctarget, dimensions))
    else:
        nalloctarget = ntarget = len(target)

    pottarg = np.zeros((nalloctarget), dtype=np.complex128)
    fldtarg = np.zeros((nalloctarget, dimensions), dtype=np.complex128)

    if slp_density is not None:
        ifcharge = 1
    else:
        ifcharge = 0
        slp_density = np.empty(size)

    if dlp_density is not None:
        ifdipole = 1
        if dipvec is None:
            raise TypeError("must specify dipvec if requesting dipole density")
    else:
        ifdipole = 0
        dlp_density = np.empty(size)

    if not hasattr(slp_density, "shape") or slp_density.shape == ():
        slp_density = np.ones(size, dtype=np.float64)*slp_density
    if not hasattr(dlp_density, "shape") or dlp_density.shape == ():
        dlp_density = np.ones(size, dtype=np.float64)*dlp_density

    if dipvec is None:
        dipvec = np.empty((dimensions, size), order="F")
    elif len(dipvec.shape) == 1:
        dipvec_new = np.empty((dimensions, size), order="F")
        if dipvec.dtype == object:
            for i in range(len(dipvec)):
                dipvec_new[i] = dipvec[i]
        else:
            dipvec_new[:] = dipvec[:, np.newaxis]

        dipvec = dipvec_new
    else:
        dipvec = dipvec.T

    use_hessian_codepath = force_hess or "h" in what.lower()

    if dimensions == 2:
        use_hessian_codepath = True

    ifcharge *= kernel.flag_value
    ifdipole *= kernel.flag_value

    if not use_hessian_codepath:
        # {{{ no Hessians, no second derivatives

        # {{{ process kernel argument

        if isinstance(kernel, HelmholtzKernel):
            #kind = kind.replace("tria", "trif")
            pass
        elif isinstance(kernel, DifferenceKernel):
            if kind != "tria":
                raise RuntimeError("difference kernel only supported "
                        "on triangles")
            #kind = "trif"

        if not isinstance(kernel,
                (LaplaceKernel, HelmholtzKernel, DifferenceKernel)):
            raise RuntimeError("unsupported kernel: %s" % kernel)

        # }}}

        args = [iprec] + list(kernel.k_arg) + source_args + [
                ifcharge, slp_density,
                ifdipole, dlp_density, dipvec,
                "p" in what, "g" in what, ntarget, target.T,
                "P" in what, pottarg,
                "G" in what, fldtarg.T]

        routine_name = "%sfmm%dd%starg" % (kernel.letter, dimensions, kind)

        if debug:
            print "ENTER", routine_name
        ier, pot, fld, pottarg, fldtarg = \
            getattr(_int, routine_name)(*args)
        if debug:
            print "LEAVE", routine_name

        result_dict = {
                "p": lambda: kernel_scale*pot,
                "g": lambda: (grad_scale*kernel_scale)*fld.T,
                "P": lambda: kernel_scale*pottarg,
                "G": lambda: (grad_scale*kernel_scale)*fldtarg.T
                }

        # }}}
    else:
        # {{{ with Hessians

        # {{{ process kernel argument

        if not isinstance(kernel, (LaplaceKernel, HelmholtzKernel)):
            raise RuntimeError("unsupported kernel: %s" % kernel)

        # }}}

        hessian_components = {2: 3, 3: 6}
        hesstarg = np.zeros(
                (nalloctarget, hessian_components[dimensions]),
                dtype=np.complex128)

        if dimensions == 3:
            hess_str = "hess"
        else:
            hess_str = ""

        routine_name = "%sfmm%dd%s%starg" % (
                kernel.letter, dimensions, kind, hess_str)
        args = [iprec] + list(kernel.k_arg) + source_args + [
                ifcharge, slp_density,
                ifdipole, dlp_density, dipvec]

        args = args + [
                "p" in what, "g" in what, "h" in what, ntarget, target.T,
                "P" in what, pottarg,
                "G" in what, fldtarg.T,
                "H" in what, hesstarg.T,
                ]

        if debug:
            print "ENTER", routine_name
        ier, pot, fld, hess, pottarg, fldtarg, hesstarg = \
                getattr(_int, routine_name)(*args)
        if debug:
            print "LEAVE", routine_name

        result_dict = {
                "p": lambda: kernel_scale*pot,
                "g": lambda: (grad_scale*kernel_scale)*fld.T,
                "h": lambda: kernel_scale*hess.T,
                "P": lambda: kernel_scale*pottarg,
                "G": lambda: (grad_scale*kernel_scale)*fldtarg.T,
                "H": lambda: kernel_scale*hesstarg.T
                }

        # }}}

    if ier:
        from warnings import warn
        warn("FMM routine '%s' encountered an error (code %d)" % (
            routine_name, ier))

    result = [result_dict[wch]() for wch in what.replace(",", "")]

    if len(what) == 1:
        return result[0]
    else:
        return result


def fmm_part(what, iprec, kernel, sources, mop_charge=None, dip_charge=None,
        target=None, dipvec=None, **kwargs):
    """
    :param what: a string consisting of the following characters, all optional:
        p : compute potential at source
        g : compute gradient at source
        h : compute Hessians at source
        P : compute potential  at specified targets
        G : compute gradient at specified targets
        H : compute Hessians at specified targets

        Results will be returned in the same order as given in *what*.
    :param mop_charge: monopole charge
    :param dip_charge: dipole charge
    """

    _, dimensions = sources.shape
    return _fmm(dimensions, len(sources), "part", [sources.T],
            what, iprec, kernel, mop_charge, dip_charge, target, dipvec,
            **kwargs)


def fmm_tria(what, iprec, kernel, mesh, slp_density=None, dlp_density=None,
        target=None, dipvec=None, **kwargs):
    """
    :param what: a string consisting of the following characters, all optional:
        p : compute potential at source
        g : compute gradient at source
        h : compute Hessians at source
        P : compute potential  at specified targets
        G : compute gradient at specified targets
        H : compute Hessians at specified targets

        Results will be returned in the same order as given in *what*.
    :param mesh: an object expected to have attributes `triangles`, `normals`
        and `centroids`. `triangles` is a `(n,3,3)`-shape array with triangle
        vertices.  `normals` is a `(n,3)`-shape array with normals, and
        `centroids is a `(n,3)`-shape array with centroids.
    """
    force_hess = dipvec is not None

    if dipvec is None:
        dipvec = mesh.normals.T

    _, dimensions = mesh.centroids.shape

    return _fmm(dimensions, len(mesh), "tria",
            [mesh.triangles.T, mesh.normals.T, mesh.centroids.T],
            what, iprec, kernel, slp_density, dlp_density, target, dipvec,
            force_hess=force_hess, mesh=mesh, **kwargs)


# vim: foldmethod=marker
