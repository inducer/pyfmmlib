import os
import subprocess
import sys


def generate_pyf(infile, outdir):
    import pathlib

    infile = pathlib.Path(infile)
    if not infile.exists():
        raise FileNotFoundError(f"Input file does not exist: '{infile}'")

    from mako.template import Template

    from codegen_helpers import cpost, cpre
    from vec_wrappers import gen_vector_wrappers

    infile = infile.resolve()
    outfile = pathlib.Path(outdir) / (infile.with_suffix("").name)

    with open(infile, encoding="utf-8") as inf:
        tpl = Template(inf.read(), uri=str(infile), strict_undefined=True)

    result = tpl.render(
        cpre=cpre,
        cpost=cpost,
        gen_vector_wrappers=gen_vector_wrappers)

    with open(outfile, "w", encoding="utf-8") as outf:
        outf.write(result)

    return str(outfile)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str,
                        help="Path to the input file")
    parser.add_argument("-o", "--outdir", type=str,
                        help="Path to the output directory")
    args = parser.parse_args()

    if not args.infile.endswith((".pyf", ".pyf.mako", ".f.mako")):
        raise ValueError(f"Input file has unknown extension: {args.infile}")

    outdir_abs = os.path.join(os.getcwd(), args.outdir)

    # Write out the .pyf/.f file
    if args.infile.endswith(".pyf.mako"):
        fname_pyf = generate_pyf(args.infile, outdir_abs)
    else:
        fname_pyf = args.infile

    # Now invoke f2py to generate the C API module file
    if args.infile.endswith((".pyf.mako", ".pyf")):
        p = subprocess.Popen([sys.executable, "-m", "numpy.f2py", fname_pyf,
                            "--build-dir", outdir_abs],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            cwd=os.getcwd())
        out, err = p.communicate()
        if not (p.returncode == 0):
            raise RuntimeError(f"f2py failed!\n"
                            f"{out}\n"
                            r"{err}")


if __name__ == "__main__":
    main()
