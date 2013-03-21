#!/usr/bin/env python

import distribute_setup
distribute_setup.use_setuptools()

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
    tmpl = Template(open(mako_name, "rt").read(), uri=mako_name, strict_undefined=True)

    context = dict(cpre=cpre, cpost=cpost,
            #paren=paren,
            #colons=colons, to_colons=to_colons, index=index,
            #var_dim=var_dim, gen_indices=gen_indices,
            #global_state=global_state, use_modules=options.use_modules,
            #debug_verbosity=options.debug_verbosity
            gen_vector_wrappers=gen_vector_wrappers
            )
    result = tmpl.render(**context)
    open("wrappers.pyf", "wt").write(result)

# }}}

def main():
    generate_wrappers()

    import glob
    import setuptools
    from numpy.distutils.core import Extension, setup

    from os.path import exists
    if not exists("fmmlib2d") or not exists("fmmlib3d"):
        print("------------------------------------------------------------")
        print("Missing fmmlib sources")
        print("------------------------------------------------------------")
        print("You probably need to run ./grab-sources.sh.")
        print("------------------------------------------------------------")
        return

    from os.path import basename
    from glob import glob
    source_files = {}

    BLACKLIST = ["d2tstrcr_omp.f"]

    for f in glob("fmmlib2d/src/*.f") + glob("fmmlib3d/src/*.f"):
        bn = basename(f)
        if bn in BLACKLIST:
            continue

        source_files[bn] = f

    source_files = ["wrappers.pyf"] + list(source_files.values())

    conf = {}
    execfile("pyfmmlib/version.py", conf)
    setup(name="pyfmmlib",
          version=conf["VERSION_TEXT"],
          description="Python wrappers for particle FMMs",
          long_description=open("README.rst", "rt").read(),
          author="Leslie Greengard, Zydrunas Gimbutas, Andreas Kloeckner",
          author_email="inform@tiker.net",
          license = "wrapper: MIT/code: GPL2",
          url="http://mathema.tician.de/software/pymetis",
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

          packages = [ "pyfmmlib" ],
          ext_modules = [
            Extension(
              "pyfmmlib._internal",
              source_files,
              ),
            ]
         )




if __name__ == '__main__':
    main()

# vim: foldmethod=marker
