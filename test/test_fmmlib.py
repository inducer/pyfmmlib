from __future__ import division

from pyfmmlib import fmm_part
import numpy as np
import numpy.linalg as la


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

    from pyfmmlib import (h2dformmp, h2dmpmp_vec, h2dmploc_vec,
            h2dlocloc_vec, h2dtaeval_vec, hpotgrad2dall_vec, h2dmpeval_vec)

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

    assert la.norm(mp2_value - ref_value) / la.norm(ref_value) < 1e-5

    # }}}

    # {{{ local 1

    loc1_center = targets_center - np.array([1, 0])
    loc1 = h2dmploc_vec(zk, rscale, mp2_center, mp2, rscale, loc1_center, nterms)

    loc1_value, _, _ = h2dtaeval_vec(zk, rscale, loc1_center, loc1,
            ztarg=targets.T, ifgrad=False, ifhess=False)

    assert la.norm(loc1_value - ref_value) / la.norm(ref_value) < 1e-5

    # }}}

    # {{{ local 2

    loc2_center = targets_center
    loc2 = h2dlocloc_vec(zk, rscale, loc1_center, loc1, rscale, loc2_center, nterms)

    loc2_value, _, _ = h2dtaeval_vec(zk, rscale, loc2_center, loc2, ztarg=targets.T,
            ifgrad=False, ifhess=False)

    assert la.norm(loc2_value - ref_value) / la.norm(ref_value) < 1e-5

    # }}}


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        exec(sys.argv[1])
    else:
        from py.test.cmdline import main
        main([__file__])




# vim: fdm=marker
