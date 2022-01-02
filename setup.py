#!/usr/bin/env python

GENERATED_SOURCES = ["wrappers.pyf", "vec_wrappers.f90"]


# {{{ wrapper generation

def generate_wrappers():
    from vec_wrappers import gen_vector_wrappers
    from mako.template import Template
    from codegen_helpers import cpre, cpost

    base_name = "wrappers.pyf"
    mako_name = base_name + ".mako"
    with open(mako_name, "rt") as inf:
        tmpl = Template(inf.read(), uri=mako_name, strict_undefined=True)

    context = dict(cpre=cpre, cpost=cpost, gen_vector_wrappers=gen_vector_wrappers)
    result = tmpl.render(**context)
    open("wrappers.pyf", "wt").write(result)

    open("vec_wrappers.f90", "wt").write(gen_vector_wrappers())

# }}}


def count_down_delay(delay):
    from time import sleep
    import sys

    while delay:
        sys.stdout.write("Continuing in %d seconds...   \r" % delay)
        sys.stdout.flush()
        delay -= 1
        sleep(1)
    print("")


DASH_SEPARATOR = 75 * "-"


def main():
    from aksetup_helper import check_git_submodules

    check_git_submodules()

    import os

    build_mode = os.environ.get("PYFMMLIB_BUILD_MODE", "openmp-opt")

    if build_mode == "openmp-ofast":
        fopt_args = "-Ofast -fopenmp"
        opt_args = "-Ofast"
        extra_link_args = "-fopenmp"

    elif build_mode == "openmp-opt":
        fopt_args = "-O3 -fopenmp"
        opt_args = "-O3 -fopenmp"
        extra_link_args = "-fopenmp"

    elif build_mode == "debug":
        fopt_args = "-g"
        opt_args = "-g"
        extra_link_args = "-g"

    elif build_mode == "setuptools":
        fopt_args = os.environ.get("FOPT", "")
        opt_args = os.environ.get("OPT", "")
        extra_link_args = os.environ.get("extra_link_args", "")

    else:
        raise ValueError("invalid value of $PYFMMLIB_BUILD_MODE")

    # NOTE: gfortran 8.1+ errors on incorrect dummy argument shape
    # https://groups.google.com/forum/#!topic/comp.lang.fortran/x3JnAjRX-KA
    os.environ["FOPT"] = "-std=legacy {}".format(fopt_args)

    os.environ["OPT"] = opt_args
    os.environ["EXTRA_LINK_ARGS"] = extra_link_args

    try:
        import mako  # noqa
    except ImportError:
        print(DASH_SEPARATOR)
        print("Mako is not installed.")
        print(DASH_SEPARATOR)

        all_gen_exist = True
        from os.path import exists

        for s in GENERATED_SOURCES:
            if not exists(s):
                all_gen_exist = False
                break

        if all_gen_exist:
            print("pyfmmlib uses mako [1] to generate its wrappers.")
            print("All the generated files are there, so we'll continue.")
            print(DASH_SEPARATOR)
            print("Hit Ctrl-C now if you'd like to think about the situation.")
            print(DASH_SEPARATOR)

            count_down_delay(delay=5)

        else:
            print(
                "That is a problem because pyfmmlib uses mako [1] "
                "to generate its wrappers."
            )
            print("Try typing")
            print("")
            print("  pip install mako")
            print("")
            print("or")
            print("  ez_install mako")
            print("")
            print(
                "That should install it. If that doesn't work, "
                "go to the mako website,"
            )
            print("download, and install by calling 'python setup.py install'.")
            print("")
            print("[1] http://www.makotemplates.org/")
            print(DASH_SEPARATOR)
            import sys

            sys.exit(1)
    else:
        generate_wrappers()

    # looks pointless, but allows 'python setup.py develop'--do not remove
    import setuptools  # noqa
    from setuptools import find_packages

    from numpy.distutils.core import Extension, setup

    from os.path import basename
    from glob import glob

    source_files = {}

    BLOCKLIST = ["d2tstrcr_omp.f", "second-r8.f"]  # noqa

    for f in glob("fmmlib2d/src/*.f") + glob("fmmlib3d/src/*.f"):
        bn = basename(f)
        if bn in BLOCKLIST:
            continue

        source_files[bn] = f

    source_files = GENERATED_SOURCES + list(source_files.values())

    import os

    extra_link_args = os.environ.get("EXTRA_LINK_ARGS", "").split()
    if extra_link_args == [""]:
        extra_link_args = []

    ver_dic = {}
    with open("pyfmmlib/version.py") as version_file:
        version_file_contents = version_file.read()

    exec(compile(version_file_contents, "pyfmmlib/version.py", "exec"), ver_dic)

    setup(
        name="pyfmmlib",
        version=ver_dic["VERSION_TEXT"],
        description="Python wrappers for particle FMMs",
        long_description=open("README.rst", "rt").read(),
        author="Leslie Greengard, Zydrunas Gimbutas, Andreas Kloeckner",
        author_email="inform@tiker.net",
        license="wrapper: MIT/code: 3-clause BSD",
        url="http://github.com/inducer/pyfmmlib",
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Developers",
            "Intended Audience :: Other Audience",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "License :: OSI Approved :: BSD License",
            "Programming Language :: Fortran",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Mathematics",
            "Topic :: Software Development :: Libraries",
        ],
        packages=find_packages(),
        python_requires="~=3.6",
        setup_requires=[
            # https://github.com/numpy/numpy/issues/20709
            # /!\ also in requirements.txt
            "numpy != 1.22.0",
        ],
        ext_modules=[
            Extension(
                "pyfmmlib._internal",
                source_files,
                extra_link_args=extra_link_args,
            ),
        ],
    )


if __name__ == "__main__":
    main()

# vim: foldmethod=marker
