#!/usr/bin/env python

GENERATED_SOURCES = ["wrappers.pyf", "vec_wrappers.f90"]


# {{{ wrapper generation

def cpre(s):
    if s:
        return ", "+s
    else:
        return s


def cpost(s):
    if s:
        return s + ", "
    else:
        return s


def generate_wrappers():
    from vec_wrappers import gen_vector_wrappers

    from mako.template import Template

    base_name = "wrappers.pyf"
    mako_name = base_name + ".mako"
    tmpl = Template(open(mako_name, "rt").read(),
            uri=mako_name, strict_undefined=True)

    context = dict(cpre=cpre, cpost=cpost,
            gen_vector_wrappers=gen_vector_wrappers)
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
    try:
        import mako  # noqa
    except ImportError:
        print(DASH_SEPARATOR)
        print("Mako is not installed.")
        print(DASH_SEPARATOR)

        from os.path import exists
        if all(exists(s) for s in GENERATED_SOURCES):
            print("pyfmmlib uses mako [1] to generate its wrappers.")
            print("All the generated files are there, so we'll continue.")
            print(DASH_SEPARATOR)
            print("Hit Ctrl-C now if you'd like to think about the situation.")
            print(DASH_SEPARATOR)

            count_down_delay(delay=5)

        else:
            print("That is a problem because pyfmmlib uses mako [1] "
                    "to generate its wrappers.")
            print("Try typing")
            print("")
            print("  pip install mako")
            print("")
            print("or")
            print("  ez_install mako")
            print("")
            print("That should install it. If that doesn't work, "
                    "go to the mako website,")
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

    from numpy.distutils.core import Extension, setup

    from os.path import exists
    if not exists("fmmlib2d") or not exists("fmmlib3d"):
        print(DASH_SEPARATOR)
        print("Missing fmmlib sources")
        print(DASH_SEPARATOR)
        print("These can be downloaded by running ./grab-sources.sh.")
        print("Note that this will fail on Windows--it's a shell script.")
        print(DASH_SEPARATOR)
        print("I will try to do that after a short delay, unless you stop me.")
        print(DASH_SEPARATOR)
        count_down_delay(delay=10)

        from subprocess import check_call
        check_call(["./grab-sources.sh"])

    from os.path import basename
    from glob import glob
    source_files = {}

    BLACKLIST = ["d2tstrcr_omp.f", "second-r8.f"]  # noqa

    for f in glob("fmmlib2d/src/*.f") + glob("fmmlib3d/src/*.f"):
        bn = basename(f)
        if bn in BLACKLIST:
            continue

        source_files[bn] = f

    source_files = GENERATED_SOURCES + list(source_files.values())

    import os
    extra_link_args = os.environ.get("EXTRA_LINK_ARGS", "").split()
    if extra_link_args == [""]:
        extra_link_args = []

    conf = {}
    execfile("pyfmmlib/version.py", conf)
    setup(name="pyfmmlib",
          version=conf["VERSION_TEXT"],
          description="Python wrappers for particle FMMs",
          long_description=open("README.rst", "rt").read(),
          author="Leslie Greengard, Zydrunas Gimbutas, Andreas Kloeckner",
          author_email="inform@tiker.net",
          license="wrapper: MIT/code: GPL2",
          url="http://github.com/inducer/pyfmmlib",
          classifiers=[
              'Development Status :: 4 - Beta',
              'Intended Audience :: Developers',
              'Intended Audience :: Other Audience',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: MIT License',
              'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
              'Programming Language :: Fortran',
              'Programming Language :: Python',
              'Topic :: Scientific/Engineering',
              'Topic :: Scientific/Engineering :: Mathematics',
              'Topic :: Software Development :: Libraries',
              ],

          packages=["pyfmmlib"],
          install_requires=[
              "pytest>=2",
              ],
          ext_modules=[
              Extension(
                  "pyfmmlib._internal",
                  source_files,
                  extra_link_args=extra_link_args,
                  ),
              ]
          )


if __name__ == '__main__':
    main()

# vim: foldmethod=marker
