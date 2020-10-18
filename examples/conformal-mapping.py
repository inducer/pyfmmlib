#!/usr/bin/env python3
"""Demo of Theodorsen's method for conformal mapping.
"""

__copyright__ = "Copyright (C) 2020 Matt Wala"

__license__ = """
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from functools import partial
import logging

import numpy as np
import pyfmmlib


logger = logging.getLogger(__name__)


# {{{ boundary correspondence

def compute_boundary_correspondence(rho, n, atol=1e-14, maxiter=500):
    """Compute the boundary correspondence using Theodorsen's method.

    The boundary correspondence gives the boundary values f: D -> C of the
    Riemann map of a domain, where D is the unit disk.  The map is such that
    f(0)=0 and f'(0)> 0.

    The domain must be starlike with respect to the origin and obey a
    near-circularity condition to ensure convergence.

    See: https://encyclopediaofmath.org/wiki/Theodorsen_integral_equation

    Arguments:
        rho: A vectorized function for the radial part of the boundary
            parametrization in polar coordinates. Given an angular coordinate
            in [0, 2π] returns the corresponding radial coordinate.
        n: Number of discretization points to use on the unit circle
        atol: Absolute convergence tolerance
        maxiter: Bound on the number of iterations

    Returns:
        Returns the computed boundary correspondence at the discretization
        points defined on the unit circle at angles [0, 2π/n, ..., 2(n-1)π/n].
    """
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    phi = theta

    for niter in range(1, 1 + maxiter):
        phi_prev = phi
        phi = theta + harmonic_conjugate(np.log(rho(phi)))
        if np.max(np.abs(phi - phi_prev)) < atol:
            logger.info("converged after %d iterations", niter)
            break
    else:
        logger.info("stopped after %d iterations", niter)

    return rho(phi) * np.exp(1j * phi)


def harmonic_conjugate(phi):
    """Compute the discrete harmonic conjugate of the boundary values of a
    real-valued function defined on the unit circle at angles
    [0, 2π/n, ..., 2(n-1)π/n], where *n* is the length of the vector *phi*.
    """
    t = -1j * np.fft.rfft(phi)
    t[0] = 0
    return np.fft.irfft(t)

# }}}


# {{{ cauchy fmm

def cauchy_integral(source_density, targets):
    """Compute a Cauchy integral using the FMM.

    Arguments:
        source_density: Boundary values of a function that is complex analytic
            inside the unit disk, at equispaced discretization points on the
            unit circle at angles [0, 2π/n, ..., 2(n-1)π/n], wnere *n* is the
            length of *source_density*.
        targets: A 1D complex array of evaluation targets

    """
    targets = np.array([targets.real, targets.imag]).T

    n = len(source_density)
    z = np.exp(np.linspace(0, 2j * np.pi, n, endpoint=False))
    sources = np.array([z.real, z.imag]).T
    normals = sources
    tangents = np.array([-z.imag, z.real]).T

    weight = 1 / n

    def cauchy_fmm(density):
        # This splits the Cauchy kernel into the sum of two directional source
        # derivatives of the Laplace kernel. See Section 7.5 of Kress, Linear
        # Integral Equations.
        dn_part = pyfmmlib.fmm_part(
            "P", iprec=5, kernel=pyfmmlib.LaplaceKernel(),
            sources=sources,
            target=targets,
            dip_charge=weight * density,
            dipvec=normals)

        dt_part = pyfmmlib.fmm_part(
            "P", iprec=5, kernel=pyfmmlib.LaplaceKernel(),
            sources=sources,
            target=targets,
            dip_charge=weight * density,
            dipvec=tangents)

        return dn_part - 1j * dt_part

    # Barycentric Lagrange interpolation for complex analytic functions, see
    # for instance Section 3.1 of https://arxiv.org/pdf/1410.2187.pdf
    return (
        cauchy_fmm(source_density) / cauchy_fmm(np.ones_like(source_density)))

# }}}


# {{{ spider plot

CMAP = "gist_rainbow"
LINEWIDTH = 1
COLORBAR_LABELS = (r"0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$")
COLORBAR_TICKS = np.array([0, 0.5, 1, 1.5, 2]) * np.pi
# Angles corresponding to spider plot rays from origin
ANGLES = np.linspace(0, 2 * np.pi, 36, endpoint=False)
POINTS_PER_ANGLE = 100
# Radii corresponding to spider plot concentric circles
RADII = np.arange(0.1, 1, 0.1)
POINTS_PER_RADIUS = 500


def spider_plot(correspondence):
    """Produce a spider plot of a conformal mapping."""
    targets = []
    colors = []
    line_segment_sizes = []

    for angle in ANGLES:
        targets.extend(
            np.linspace(0, 1 - 1e-15, POINTS_PER_ANGLE) * np.exp(1j * angle))
        colors.extend([angle] * POINTS_PER_ANGLE)
        line_segment_sizes.append(POINTS_PER_ANGLE)

    circle = np.linspace(0, 2 * np.pi, POINTS_PER_RADIUS)
    for radius in RADII:
        targets.extend(radius * np.exp(1j * circle))
        colors.extend(circle)
        line_segment_sizes.append(POINTS_PER_RADIUS)

    assert len(targets) == len(colors) == sum(line_segment_sizes)

    targets = np.array(targets)
    mapvals = cauchy_integral(correspondence, targets)

    import matplotlib.pyplot as plt
    ax = plt.gca()
    norm = plt.Normalize(0, 2*np.pi)

    def add_line_collection(points, colors):
        # https://matplotlib.org/examples/pylab_examples/multicolored_line.html
        points = np.array([points.real, points.imag]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        from matplotlib.collections import LineCollection
        lc = LineCollection(segments, cmap=CMAP, norm=norm)
        lc.set_array(np.array(colors))
        lc.set_linewidth(LINEWIDTH)

        ax.add_collection(lc)
        return lc

    outer_circle = np.r_[correspondence, correspondence[0]]
    outer_colors = np.linspace(0, 2 * np.pi, len(outer_circle))
    lc = add_line_collection(outer_circle, outer_colors)
    colorbar = plt.colorbar(lc, ticks=COLORBAR_TICKS)
    colorbar.ax.set_yticklabels(COLORBAR_LABELS)

    start = 0
    for size in line_segment_sizes:
        add_line_collection(mapvals[start:start + size],
                            colors[start:start + size])
        start += size

    plt.axis("equal")

# }}}


# {{{ test against reference mapping

@np.vectorize
def ellipkinc(phi, m):
    """Incomplete elliptic integral of the first kind.

    Same interface as *scipy.special.ellipkinc*, but unlike the former, this
    implementation supports complex *phi*.
    """

    # See: http://www.mygeodesy.id.au/documents/Elliptic%20Integrals%20and%20Landen%27s%20Transformation.pdf  # noqa: E501
    m = np.sqrt(m)
    factor = 1 / m
    while not np.isclose(m, 1, atol=1e-15, rtol=1e-15):
        phi = (phi + np.arcsin(m * np.sin(phi))) / 2
        m = 2 * np.sqrt(m) / (1 + m)
        factor *= m

    return np.sqrt(factor) * np.log(np.tan(np.pi / 4 + phi / 2))


def K(x, k):  # noqa: N802
    """Alternative definition for the incomplete elliptic integral of the first
    kind.
    """
    return ellipkinc(np.arcsin(x), k**2)


@np.vectorize
def ellipse_ref(z):
    """Reference conformal map from the disk to the ellipse (x/2)^2 + y^2 = 1.

    See: https://www.emis.de/journals/AASF/Vol31/kanas.pdf
    """

    # xi is the parameter in the equation (u/cosh(xi))^2 + (v/sinh(xi))^2 = 1.
    xi = np.log(3) / 2
    major_radius = np.cosh(xi)
    scale = 2 / major_radius

    # s is the root of mu(t) - 2 * xi in (0, 1),
    # where mu(t) = pi / 2 * K(1, sqrt(1 - t**2)) / K(1, t)
    s = 0.9142838686166854

    z /= np.sqrt(s)
    if np.isreal(z) and np.abs(z) > 1:
        # The elliptic integral has a branch cut for real z with |z| > 1.
        # The output doesn't depend on the choice of side here.
        z += 1e-100j

    return np.sin(np.pi * K(z, s) / (2 * K(1, s))) * scale


def test_versus_ellipse_ref(visualize=True):
    """Tests against an analytical mapping of the ellipse."""
    r = 2

    def rho(x):
        return r / np.sqrt(np.cos(x) ** 2 + (r * np.sin(x)) ** 2)

    n = 2048
    correspondence = compute_boundary_correspondence(rho, n)

    points = np.exp(np.linspace(0, 2j*np.pi, n, endpoint=False))
    ref_correspondence = ellipse_ref(points)

    logger.info("Number of discretization points: %d", n)

    logger.info("Absolute error on boundary: %e",
          np.max(np.abs(correspondence - ref_correspondence)))

    assert np.allclose(correspondence, ref_correspondence)

    test_pts = 0.9 * np.exp(np.linspace(0, 2j * np.pi, 100, endpoint=False))
    map_vals_at_test_pts = cauchy_integral(correspondence, test_pts)
    ref_map_vals_at_test_pts = cauchy_integral(ref_correspondence, test_pts)

    logger.info("Absolute error at test points: %e",
          np.max(np.abs(map_vals_at_test_pts - ref_map_vals_at_test_pts)))

    assert np.allclose(map_vals_at_test_pts, ref_map_vals_at_test_pts)

    if visualize:
        import matplotlib.pyplot as plt
        spider_plot(correspondence)
        outfile = "ellipse-map.png"
        plt.savefig(outfile)
        logger.info("Wrote '%s'", outfile)

# }}}


# {{{ example parametrizations

def rho_ellipse(x, r=2):
    return r / np.sqrt(np.cos(x) ** 2 + (r * np.sin(x)) ** 2)


def rho_n_gon(n, x):
    # https://math.stackexchange.com/q/123548
    c = 2*np.pi
    return 1/np.cos(c/n * (n*x/c - np.floor(n*x/c)) - c/(2*n))


rho_square = partial(rho_n_gon, 4)
rho_pentagon = partial(rho_n_gon, 5)


def rho_peanut(x):
    return ((2*np.cos(x))**2 + np.sin(x)**2)**(1/2)


def rho_epitrochoid(x):
    return 1 + 0.18 * np.cos(5 * x)


def rho_cardioid(x):
    return 1 + 0.6 * np.cos(x)


def rho_cranioid(x):
    return np.sin(x) + 2 * np.sqrt(1 - 0.6 * np.cos(x)**2)

# }}}


# {{{ solve and plot example

def show_example(rho, n=8192):
    import matplotlib.pyplot as plt
    corr = compute_boundary_correspondence(rho, n)
    spider_plot(corr)
    plt.show()

# }}}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    try:
        import matplotlib  # noqa: F401
    except ImportError:
        from warnings import warn
        warn("matplotlib not installed, not visualizing")
        visualize = False
    else:
        visualize = True

    test_versus_ellipse_ref(visualize=visualize)
    # show_example(rho_peanut)
