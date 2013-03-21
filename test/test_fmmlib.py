from __future__ import division

from pyfmmlib import fmm_part
import numpy as np




def test_fmm():
    for dims in [2, 3]:
        for kernel in [0, 5]:
            sources = np.random.randn(4000, dims)
            dipvec = np.random.randn(4000, dims)
            fmm_part("pg", iprec=1, kernel=kernel, sources=sources, dip_charge=1, dipvec=dipvec,
                    debug=True)

            targ_def = (slice(-3, 3, 20j),)
            targets = np.mgrid[targ_def*dims]
            targets = targets.reshape(dims, -1)

            fmm_part("PG", iprec=1, kernel=kernel,
                    sources=sources, mop_charge=1, target=targets.T,
                    debug=True)




if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        exec(sys.argv[1])
    else:
        from py.test.cmdline import main
        main([__file__])




# vim: fdm=marker
