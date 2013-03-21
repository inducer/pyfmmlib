from __future__ import division

from pyfmmlib import fmm_part, HelmholtzKernel
import numpy as np
import matplotlib.pyplot as pt

sources = np.random.randn(40, 2)

targets = np.mgrid[-7:7:400j, -7:7:400j]
pot_shape = targets.shape[1:]
targets = targets.reshape(2, -1)

pot, grad = fmm_part("PG", iprec=2, kernel=HelmholtzKernel(5),
        sources=sources, mop_charge=1, target=targets.T)

pot = pot.reshape(pot_shape)

pt.imshow(pot.real)
pt.show()
