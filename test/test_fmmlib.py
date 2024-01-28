
import numpy as np
import numpy.linalg as la
from scipy.integrate import dblquad

from pyfmmlib import LaplaceKernel, fmm_part, fmm_tria


def test_fmm():
    for dims in [2, 3]:
        for kernel in [0, 5]:
            sources = np.random.randn(4000, dims)
            dipvec = np.random.randn(4000, dims)
            fmm_part("pg", iprec=1, kernel=kernel, sources=sources,
                    dip_charge=1, dipvec=dipvec,
                    debug=True)

            targ_def = (slice(-3, 3, 20j),)
            targets = np.mgrid[targ_def*dims]
            targets = targets.reshape(dims, -1)

            fmm_part("PG", iprec=1, kernel=kernel,
                    sources=sources, mop_charge=1, target=targets.T,
                    debug=True)


def test_triangle():
    n = 3
    triangles = np.random.rand(n, 3, 3)
    centroids = np.mean(triangles, axis=1)

    normals = np.cross(
        triangles[:, 2]-triangles[:, 0],
        triangles[:, 1]-triangles[:, 0])

    normals /= np.linalg.norm(normals, axis=0)

    charges = np.random.rand(n)

    class Mesh:
        def __init__(self, triangles, normals, centroids):
            self.triangles = triangles
            self.normals = normals
            self.centroids = centroids

        def __len__(self):
            return len(self.triangles)

    # Return the exact value for triangle at index i,
    # using a slow but accurate adaptive integration
    def exact(i):

        target = centroids[i]

        def to_integrate(a, b, j):
            v1, v2, v3 = triangles[j]

            location = v1 + (v3-v1)*a + (v2-v1)*b
            r = np.linalg.norm(location-target)
            area = np.linalg.norm(np.cross(v3-v1, v2-v1))

            return area/(4*np.pi*r)

        potential = 0.

        for j in range(n):
            potential += charges[j]*dblquad(to_integrate, 0, 1, 0,
                lambda y: 1-y, args=(j,))[0]

        return potential

    m = Mesh(triangles, normals, centroids)
    result_fmm = fmm_tria("p", 0, LaplaceKernel(), m, slp_density=charges).real[0]
    result_exact = exact(0)

    assert np.isclose(result_fmm, result_exact)


def test_translations():
    nterms = 15
    zk = 3
    rscale = 1

    n = 40
    # centered at the origin, extent [-.5,.5]
    sources = np.random.uniform(size=(n, 2)) - 0.5
    charges = np.random.uniform(size=n)

    targets_center = np.array([10, 0])
    targets = np.random.uniform(size=(n, 2)) - 0.5 + targets_center

    from pyfmmlib import (
        h2dformmp, h2dlocloc_vec, h2dmpeval_vec, h2dmploc_vec, h2dmpmp_vec,
        h2dtaeval_vec, hpotgrad2dall_vec)

    ref_value, _, _ = hpotgrad2dall_vec(ifgrad=False, ifhess=False,
            sources=sources.T, charge=charges,
            targets=targets.T, zk=zk)

    # {{{ multipole 1

    mp1_center = np.array([0, 0])
    ier, mp1 = h2dformmp(zk, rscale, sources.T, charges, mp1_center, nterms)
    assert ier == 0

    mp1_value, _, _ = h2dmpeval_vec(zk, rscale, mp1_center, mp1, ztarg=targets.T,
            ifgrad=False, ifhess=False)

    assert la.norm(mp1_value - ref_value) / la.norm(ref_value) < 1e-12

    # }}}

    # {{{ multipole 2

    mp2_center = np.array([2, 0])
    mp2 = h2dmpmp_vec(zk, rscale, mp1_center, mp1, rscale, mp2_center, nterms)

    mp2_value, _, _ = h2dmpeval_vec(zk, rscale, mp2_center, mp2, ztarg=targets.T,
            ifgrad=False, ifhess=False)

    assert la.norm(mp2_value - ref_value) / la.norm(ref_value) < 3e-5

    # }}}

    # {{{ local 1

    loc1_center = targets_center - np.array([1, 0])
    loc1 = h2dmploc_vec(zk, rscale, mp2_center, mp2, rscale, loc1_center, nterms)

    loc1_value, _, _ = h2dtaeval_vec(zk, rscale, loc1_center, loc1,
            ztarg=targets.T, ifgrad=False, ifhess=False)

    assert la.norm(loc1_value - ref_value) / la.norm(ref_value) < 3e-5

    # }}}

    # {{{ local 2

    loc2_center = targets_center
    loc2 = h2dlocloc_vec(zk, rscale, loc1_center, loc1, rscale, loc2_center, nterms)

    loc2_value, _, _ = h2dtaeval_vec(zk, rscale, loc2_center, loc2, ztarg=targets.T,
            ifgrad=False, ifhess=False)

    assert la.norm(loc2_value - ref_value) / la.norm(ref_value) < 1e-4

    # }}}


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        exec(sys.argv[1])
    else:
        from pytest import main
        main([__file__])

# vim: fdm=marker
