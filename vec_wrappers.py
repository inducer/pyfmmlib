import re

def parse_args(args):
    args_re = re.compile(r"^\s*([a-z]+(?:\s*\*\s*[0-9]+)?)\s+(.*)\s*$")
    array_re = re.compile(r"([a-zA-Z][a-zA-Z0-9]*)\(([-+*a-z:0-9,]+)\)")
    scalar_re = re.compile(r"([a-zA-Z][a-zA-Z0-9]*)")

    for line in args.split("\n"):
        if not line.strip():
            continue

        args_match = args_re.match(line)
        assert args_match, line
        type_str = args_match.group(1)
        names_and_shapes = args_match.group(2)

        while names_and_shapes.strip():
            array_match = array_re.match(names_and_shapes)
            scalar_match = scalar_re.match(names_and_shapes)
            if array_match is not None:
                yield type_str, array_match.group(1), tuple(
                        x.strip() for x in array_match.group(2).split(","))
                names_and_shapes = names_and_shapes[array_match.end():]
            elif scalar_match is not None:
                yield type_str, scalar_match.group(1), ()
                names_and_shapes = names_and_shapes[scalar_match.end():]
            else:
                raise RuntimeError("arg parsing did not understand: %s"
                        % names_and_shapes)

            names_and_shapes = names_and_shapes.strip()
            if names_and_shapes.startswith(","):
                names_and_shapes = names_and_shapes[1:]
                names_and_shapes = names_and_shapes.strip()
            else:
                if names_and_shapes:
                    raise RuntimeError("comma expected")




def get_vector_wrapper(func_name, args, out_args, vec_func_name=None,
        arg_order=None, too_many_ok=False):
    if vec_func_name is None:
        vec_func_name = func_name+"_vec"

    if isinstance(args, str):
        args = list(parse_args(args))

    if arg_order is not None:
        if isinstance(arg_order, str):
            arg_order = [x.strip() for x in arg_order.split(",")]

        arg_dict = dict((name, (type_, shape))
                for type_, name, shape in args)

        args = []
        for arg in arg_order:
            type_, shape = arg_dict.pop(arg)
            args.append((type_, arg, shape))

        if not too_many_ok and arg_dict:
            raise RuntimeError("too many args: %s" % ",".join(arg_dict))

    all_args = set(name for type_, name, shape in args)
    in_args = all_args-set(out_args)

    yield "subroutine %s(%s, nvcount)" % (
            vec_func_name, ", ".join(name for type_, name, shape in args))

    yield "  implicit none"
    yield "  integer, intent(in) :: nvcount"
    yield "  integer ivcount"

    for type_, name, shape in args:
        intent = "in" if name in in_args else "out"
        if shape:
            yield "  %s, intent(%s) :: %s(%s)" % (
                    type_, intent, name, ",".join(str(si) for si in shape))
        else:
            yield "  %s, intent(%s) :: %s" % (type_, intent, name)


    # assemble call_args
    call_args = []

    def gen_first_index(shape_dim):
        if str(shape_dim) == "nvcount":
            return "ivcount"

        colon_idx = str(shape_dim).find(":")
        if colon_idx != -1:
            return shape_dim[:colon_idx]
        else:
            return "1"

    for type_, name, shape in args:
        if not shape or not any("nvcount" in shape_dim for shape_dim in shape):
            call_args.append(name)
        else:
            call_args.append("%s(%s)" % (
                name, ", ".join(gen_first_index(shape_dim)
                    for shape_dim in shape)))

    # generate loop
    yield ""
    yield "  !$omp parallel do"
    yield "  do ivcount = 1, nvcount"

    call_args_txt = ""
    call_args_line = "      "
    while call_args:
        if len(call_args_line) >= 70:
            call_args_txt += call_args_line + " &\n"
            call_args_line = "      "

        call_args_line += call_args.pop(0)
        if call_args:
            call_args_line += ", "

    call_args_txt += call_args_line

    yield "    call %s( &\n%s &\n      )" % (func_name, call_args_txt)
    yield "  enddo"
    yield "  !$omp end parallel do"
    yield ""

    yield "  return"
    yield "end"
    yield ""

def gen_vector_wrappers():
    result = []

    def gen_vector_wrapper(*args, **kwargs):
        for line in get_vector_wrapper(*args, **kwargs):
            result.append(line)

    gen_vector_wrapper("triangle_norm", """
            real*8 triangles(3,3,nvcount)
            real*8 trinorm(3,nvcount)
            """, ["trinorm"])

    gen_vector_wrapper("triangle_area", """
            real*8 triangles(3,3,nvcount)
            real*8 triarea(nvcount)
            """, ["triarea"])

    gen_vector_wrapper("ylgndr", """
            integer nmax
            real *8 x(nvcount), y(0:nmax,0:nmax,nvcount)
            """, ["y"])

    for dp_or_no in ["", "_dp"]:
        for what in ["l", "h"]:
            for dims in [2,3]:
                if dims == 2:
                    fld_or_grad = "grad"
                    hess_dims = 3
                    hess_or_no = ",hess"
                    ifhess_or_no = ",ifhess"
                else:
                    fld_or_grad = "fld"
                    hess_dims = 6
                    hess_or_no = ""
                    ifhess_or_no = ""

                if what == "l":
                    wavek_or_no = ""
                else:
                    wavek_or_no = ",wavek"

                if dp_or_no:
                    charge_or_dip = "dipstr,dipvec"
                else:
                    charge_or_dip = "charge"

                gen_vector_wrapper("%(what)spot%(fld_or_grad)s%(dims)ddall%(dp_or_no)s" 
                        % locals(), 
                        """
                      integer if%(fld_or_grad)s,ifhess,nsources
                      real *8 sources(%(dims)d,nsources),targets(%(dims)d,nvcount)
                      complex *16 charge(nsources)
                      complex *16 dipstr(nsources)
                      real*8 dipvec(%(dims)d,nsources)
                      complex *16 wavek,pot(nvcount),%(fld_or_grad)s(%(dims)d,nvcount)
                      complex *16 hess(%(hess_dims)d,nvcount)
                      """ % locals(), ["pot", fld_or_grad, "hess"], 
                      arg_order=("if%(fld_or_grad)s%(ifhess_or_no)s,sources,%(charge_or_dip)s,nsources,targets%(wavek_or_no)s,"
                          "pot,%(fld_or_grad)s%(hess_or_no)s")
                      % locals(), too_many_ok=True)

    gen_vector_wrapper("h3dtaeval", """
            complex*16 wavek
            real*8 rscale
            real*8 center(3)
            complex*16 locexp(0:nterms,-nterms:nterms)
            integer nterms
            real*8 ztarg(3,nvcount)
            complex*16 pot(nvcount)
            integer iffld
            complex*16 fld(3,nvcount)
            integer ier(nvcount)
            """, ["ier", "pot", "fld"])

    for what, extra_args in [
            ("l", ""),
            ("h", "complex*16 wavek")
            ]:
        gen_vector_wrapper("%s2dtaeval" % what, """
                %s
                real*8 rscale
                real*8 center(2)
                complex*16 mpole(-nterms:nterms)
                integer nterms
                real*8 ztarg(2,nvcount)
                complex*16 pot(nvcount)
                integer ifgrad
                complex*16 grad(2,nvcount)
                integer ifhess
                complex*16 hess(3,nvcount)
                """ % extra_args, ["pot", "grad", "hess"])

        hess_output = ""
        taeval_out_args = ["pot", "grad", "hess", "ier"]
        taeval_func_name ="%s3dtaeval" % what
        if what == "l":
            hess_output = """
                integer ifhess
                complex*16 hess(6,nvcount)
                """

            taeval_func_name += "hess"
            taeval_out_args.append("hess")

        gen_vector_wrapper(taeval_func_name, """
                %s
                real*8 rscale(nvcount)
                real*8 center(3,nvcount)
                complex *16 mpole(0:nterms,-nterms:nterms,nvcount)
                integer nterms
                real*8 ztarg(3,nvcount)
                complex*16 pot(nvcount)
                integer ifgrad
                complex*16 grad(3,nvcount)
                %s
                integer ier
                """ % (extra_args, hess_output), taeval_out_args,
                vec_func_name=taeval_func_name+"_1tgtperexp")


    gen_vector_wrapper("hank103", """
            complex*16 z(nvcount), h0(nvcount), h1(nvcount)
            integer ifexpon
            """, ["h0", "h1"])

    gen_vector_wrapper("legefder", """
            real*8 x(nvcount)
            real*8 val(nvcount)
            real*8 der(nvcount)
            real*8 pexp(n+1)
            integer n
            """, ["val", "der"])

    return "\n".join(result)
