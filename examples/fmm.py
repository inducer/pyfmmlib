#!/usr/bin/env python3
"""Computes a Helmholtz potential with the FMM.
"""

import numpy as np

from pyfmmlib import HelmholtzKernel, fmm_part


def main():
    sources = np.random.randn(40, 2)

    targets = np.mgrid[-7:7:400j, -7:7:400j]
    pot_shape = targets.shape[1:]
    targets = targets.reshape(2, -1)

    pot, grad = fmm_part("PG", iprec=2, kernel=HelmholtzKernel(5),
            sources=sources, mop_charge=1, target=targets.T)

    pot = pot.reshape(pot_shape)

    try:
        import matplotlib.pyplot as pt
    except ImportError:
        print("matplotlib not installed, not plotting")
    else:
        pt.imshow(pot.real)
        outfile = "helmholtz-potential.png"
        pt.savefig(outfile)
        print("wrote '%s'" % outfile)


if __name__ == "__main__":
    main()
