import re
import sys
from mako.template import Template
import functools
import shlex


if sys.version_info < (3,):
    string_types = basestring  # noqa
else:
    string_types = str


# {{{ generation helpers

def parse_args(args):
    args_re = re.compile(r"^\s*([a-z]+(?:\s*\*\s*[0-9]+)?)\s+(.*)\s*$")
    array_re = re.compile(r"([a-zA-Z][_a-zA-Z0-9]*)\(([-+*_a-zA-Z:0-9,() ]+)\)$")
    scalar_re = re.compile(r"([a-zA-Z][_a-zA-Z0-9]*)$")

    for line in args.split("\n"):
        if not line.strip():
            continue

        args_match = args_re.match(line)
        assert args_match, line
        type_str = args_match.group(1)
        names_and_shapes = args_match.group(2)

        array_match = array_re.match(names_and_shapes)
        scalar_match = scalar_re.match(names_and_shapes)
        if array_match is not None:
            shape = tuple(x.strip() for x in array_match.group(2).split(","))

            yield type_str, array_match.group(1), shape
        elif scalar_match is not None:
            yield type_str, scalar_match.group(1), ()
        else:
            raise RuntimeError("arg parsing did not understand: %s"
                    % names_and_shapes)


INDIRECT_MARKER = "*INDIRECT"
MANY_MARKER = "*INDIRECT_MANY"
INDIRECT_MARKERS = [INDIRECT_MARKER, MANY_MARKER]


def shape_has_indirect(shape):
    return any(s_i in INDIRECT_MARKERS for s_i in shape)


def with_sub(name, sub):
    if not sub:
        return name
    else:
        return "%s(%s)" % (name, ", ".join(sub))


def pad_fortran(line, width):
    line += ' ' * (width - 1 - len(line))
    line += '&'
    return line


def wrap_line_base(line, level=0, width=80, indentation='    ',
                   pad_func=lambda string, amount: string,
                   lex_func=functools.partial(shlex.split, posix=False)):
    """
    The input is a line of code at the given indentation level. Return the list
    of lines that results from wrapping the line to the given width. Lines
    subsequent to the first line in the returned list are padded with extra
    indentation. The initial indentation level is not included in the input or
    output lines.

    The `pad_func` argument is a function that adds line continuations. The
    `lex_func` argument returns the list of tokens in the line.
    """
    tokens = lex_func(line)
    resulting_lines = []
    at_line_start = True
    indentation_len = len(level * indentation)
    current_line = ''
    padding_width = width - indentation_len
    for index, word in enumerate(tokens):
        has_next_word = index < len(tokens) - 1
        word_len = len(word)
        if not at_line_start:
            next_len = indentation_len + len(current_line) + 1 + word_len
            if next_len < width or (not has_next_word and next_len == width):
                # The word goes on the same line.
                current_line += ' ' + word
            else:
                # The word goes on the next line.
                resulting_lines.append(pad_func(current_line, padding_width))
                at_line_start = True
                current_line = indentation
        if at_line_start:
            current_line += word
            at_line_start = False
    resulting_lines.append(current_line)
    return resulting_lines


wrap_line = functools.partial(wrap_line_base, pad_func=pad_fortran)


def generate_loop(func_name, args, out_args, has_indirect_many,
        output_reductions, tmp_init):
    ind = 4*" "
    yield ind + "do ivcount = 1, nvcount"

    if has_indirect_many:
        for type_, name, shape in args:
            if shape and MANY_MARKER in shape:
                yield (ind + "  ncsr_count = "
                        "%(name)s_starts(ivcount+1) "
                        "- %(name)s_starts(ivcount)"
                        % {"name": name})

                break
                # FIXME: Check that other starts yield the
                # same count.

        yield ind + "  do icsr = 0, ncsr_count-1"

    # {{{ assemble call_args

    call_args = []

    def gen_first_index(name, shape_dim):
        if str(shape_dim) == "nvcount":
            return "ivcount"

        if str(shape_dim) == INDIRECT_MARKER:
            return "%s_offsets(ivcount)" % name
        if str(shape_dim) == MANY_MARKER:
            return ("%(name)s_offsets(%(name)s_starts(ivcount) + icsr)"
                    % {"name": name})

        colon_idx = str(shape_dim).find(":")
        if colon_idx != -1:
            return shape_dim[:colon_idx]
        else:
            return "1"

    for type_, name, shape in args:
        if has_indirect_many and name in out_args and MANY_MARKER not in shape:
            call_args.append(with_sub(name + "_tmp",
                [gen_first_index(name, shape_dim)
                    for shape_dim in shape
                    if shape_dim != "nvcount"]))
        elif not ("nvcount" in shape or shape_has_indirect(shape)):
            call_args.append(name)
        else:
            call_args.append("%s(%s)" % (
                name, ", ".join(gen_first_index(name, shape_dim)
                    for shape_dim in shape)))

    # }}}

    call_ind = ind + "    "

    if has_indirect_many:
        for type_, name, shape in args:
            if (has_indirect_many and
                    name in out_args and
                    MANY_MARKER not in shape):
                tmp = name + "_tmp"

                if name in tmp_init:
                    yield call_ind + "%s = %s" % (tmp, tmp_init[name])

    for l in wrap_line(
            "%scall %s(%s)" % (call_ind, func_name, ", ".join(call_args)),
            indentation="  "):
        yield call_ind + l

    if has_indirect_many:
        for type_, name, shape in args:
            if (has_indirect_many and
                    name in out_args and
                    MANY_MARKER not in shape):
                tgt_sub = [
                        ":" if shape_dim != "nvcount" else "ivcount"
                        for shape_dim in shape
                        ]

                tgt = with_sub(name, tgt_sub)
                tmp = name + "_tmp"

                out_red = output_reductions[name]
                if out_red == "sum":
                    yield call_ind + "%s = %s + %s" % (tgt, tgt, tmp)
                elif out_red == "max":
                    yield call_ind + "%s = max(%s, %s)" % (tgt, tgt, tmp)

                else:
                    raise ValueError("invalid output reduction: %s" % out_red)

        yield ind + "  enddo"

    yield ind + "enddo"


def get_vector_wrapper(func_name, args, out_args, vec_func_name=None,
        arg_order=None, too_many_ok=False,
        output_reductions=None, tmp_init=None, omp_chunk_size=10):
    if vec_func_name is None:
        vec_func_name = func_name+"_vec"

    # {{{ process args/arg_order

    if isinstance(args, string_types):
        args = list(parse_args(args))

    if arg_order is not None:
        if isinstance(arg_order, string_types):
            arg_order = [x.strip() for x in arg_order.split(",")]

        arg_dict = dict((name, (type_, shape))
                for type_, name, shape in args)

        args = []
        for arg in arg_order:
            type_, shape = arg_dict.pop(arg)
            args.append((type_, arg, shape))

        if not too_many_ok and arg_dict:
            raise RuntimeError("too many args: %s" % ",".join(arg_dict))

    del arg_order
    del too_many_ok

    # }}}

    all_args = set(name for type_, name, shape in args)
    in_args = all_args-set(out_args)

    has_indirect = False
    has_indirect_many = False

    passed_args_names = []
    for type_, name, shape in args:
        passed_args_names.append(name)
        if shape_has_indirect(shape):
            passed_args_names.append(name+"_offsets")
            has_indirect = True
            if MANY_MARKER in shape:
                has_indirect_many = True
                passed_args_names.append(name+"_starts")

    assert not (has_indirect_many and output_reductions is None)
    assert not (not has_indirect_many and output_reductions is not None)
    assert not (not has_indirect_many and tmp_init is not None)

    if tmp_init is None:
        tmp_init = {}

    # {{{ code generation

    for l in wrap_line(
            "subroutine %s(%s)" % (
                vec_func_name,
                ", ".join(passed_args_names + ["nvcount"]))):
        yield l

    yield "  implicit none"
    yield "  integer, intent(in) :: nvcount"
    yield "  integer ivcount"

    if has_indirect_many:
        yield "  integer ncsr_count"
        yield "  integer icsr"

    for type_, name, shape in args:
        intent = "in" if name in in_args else "out"
        if shape:
            processed_shape = ["0:*" if s_i in INDIRECT_MARKERS else s_i
                    for s_i in shape]

            if (has_indirect_many and
                    name in out_args and
                    MANY_MARKER not in shape):
                yield "  %s %s(%s)" % (
                        type_, name, ",".join(str(si) for si in processed_shape))
                yield "  !f2py intent(in,out) %s" % name

                tmp_shape = [s_i for s_i in shape if s_i != "nvcount"]
                yield "  %s :: %s" % (
                        type_, with_sub(name+"_tmp", tmp_shape))
            else:
                yield "  %s, intent(%s) :: %s(%s)" % (
                        type_, intent, name, ",".join(
                            str(si) for si in processed_shape))

            if INDIRECT_MARKER in shape:
                yield "  integer, intent(in) :: %s_offsets(nvcount)" % name
            if MANY_MARKER in shape:
                yield "  integer, intent(in) :: %s_offsets(0:*)" % name
                yield "  integer, intent(in) :: %s_starts(nvcount+1)" % name

        else:
            yield "  %s, intent(%s) :: %s" % (type_, intent, name)

    extra_omp = ""
    if has_indirect:
        extra_omp = " schedule(dynamic, %d)" % omp_chunk_size
    else:
        extra_omp = " schedule(static, %d)" % omp_chunk_size

    if has_indirect_many:
        private_vars = ["icsr", "ncsr_count"]
        for type_, name, shape in args:
            if shape and name in out_args and MANY_MARKER not in shape:
                private_vars.append(name + "_tmp")

        extra_omp += " private(%s)" % ", ".join(private_vars)

    shared_vars = ["nvcount"]
    for type_, name, shape in args:
        shared_vars.append(name)
        if shape and MANY_MARKER in shape:
            shared_vars.append(name + "_offsets")
            shared_vars.append(name + "_starts")
        if shape and INDIRECT_MARKER in shape:
            shared_vars.append(name + "_offsets")

    extra_omp += " shared(%s)" % ", ".join(shared_vars)

    # generate loop
    yield ""
    yield "  if (nvcount .le. %d) then" % omp_chunk_size
    for l in generate_loop(func_name, args, out_args, has_indirect_many,
            output_reductions, tmp_init):
        yield l
    yield "  else"
    for l in wrap_line(
            "!$omp parallel do default(none)" + extra_omp,
            indentation="!$omp "):
        yield "    " + l

    for l in generate_loop(func_name, args, out_args, has_indirect_many,
            output_reductions, tmp_init):
        yield l
    yield "    !$omp end parallel do"
    yield "  endif"

    yield "  return"
    yield "end"
    yield ""

    # }}}

# }}}


def gen_vector_wrappers():
    result = []

    def gen_vector_wrapper(*args, **kwargs):
        for line in get_vector_wrapper(*args, **kwargs):
            result.append(line)

    import codegen_helpers as cgh

    # {{{ helpers

    gen_vector_wrapper("triangle_norm", """
            real*8 triangles(3,3,nvcount)
            real*8 trinorm(3,nvcount)
            """, ["trinorm"])

    gen_vector_wrapper("triangle_area", """
            real*8 triangles(3,3,nvcount)
            real*8 triarea(nvcount)
            """, ["triarea"])

    # }}}

    # {{{ special functions

    gen_vector_wrapper("ylgndr", """
            integer nmax
            real *8 x(nvcount)
            real *8 y(0:nmax,0:nmax,nvcount)
            """, ["y"])

    gen_vector_wrapper("hank103", """
            complex*16 z(nvcount)
            complex*16 h0(nvcount)
            complex*16 h1(nvcount)
            integer ifexpon
            """, ["h0", "h1"])

    gen_vector_wrapper("legefder", """
            real*8 x(nvcount)
            real*8 val(nvcount)
            real*8 der(nvcount)
            real*8 pexp(n+1)
            integer n
            """, ["val", "der"])

    # }}}

    # {{{ direct evaluation

    for dp_or_no in ["", "_dp"]:
        for what in ["l", "h"]:
            for dims in [2, 3]:
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
                    wavek_or_no = ",zk"

                if dp_or_no:
                    charge_or_dip = "dipstr,dipvec"
                else:
                    charge_or_dip = "charge"

                gen_vector_wrapper(
                        "%(what)spot%(fld_or_grad)s%(dims)ddall%(dp_or_no)s"
                        % locals(),
                        """
                        integer if%(fld_or_grad)s
                        integer ifhess
                        integer nsources
                        real *8 sources(%(dims)d,nsources)
                        real *8 targets(%(dims)d,nvcount)
                        complex *16 charge(nsources)
                        complex *16 dipstr(nsources)
                        real*8 dipvec(%(dims)d,nsources)
                        complex *16 zk
                        complex *16 pot(nvcount)
                        complex *16 %(fld_or_grad)s(%(dims)d,nvcount)
                        complex *16 hess(%(hess_dims)d,nvcount)
                        """ % locals(), ["pot", fld_or_grad, "hess"],
                        arg_order=(
                            "if%(fld_or_grad)s%(ifhess_or_no)s,sources,"
                            "%(charge_or_dip)s,"
                            "nsources,targets%(wavek_or_no)s,"
                            "pot,%(fld_or_grad)s%(hess_or_no)s") % locals(),
                        too_many_ok=True)

    # }}}

    # {{{ formta

    for dp_or_no in ["", "_dp"]:
        for dims in [2, 3]:
            for eqn in [cgh.Helmholtz(dims)]:
                gen_vector_wrapper("h%ddformta%s" % (dims, dp_or_no),
                Template("""
                        integer ier(nvcount)
                        complex*16 zk
                        real*8 rscale(nvcount)
                        real *8 sources(${dims},*INDIRECT)

                        % if dp_or_no:
                            complex *16 dipstr(*INDIRECT)
                            real *8 dipvec(${dims}, *INDIRECT)
                        % else:
                            complex *16 charge(*INDIRECT)
                        % endif

                        integer nsources(nvcount)
                        real*8 center(${dims}, nvcount)
                        integer nterms
                        complex*16 locexp(${eqn.expansion_dims("nterms")},nvcount)
                        """, strict_undefined=True).render(
                            dims=dims,
                            eqn=eqn,
                            dp_or_no=dp_or_no,
                            ), ["ier", "locexp"])

    # }}}

    # {{{ {ta,mp}eval

    for expn_type in ["ta", "mp"]:
        gen_vector_wrapper("h3d%seval" % expn_type, """
                complex*16 zk
                real*8 rscale
                real*8 center(3)
                complex*16 expn(0:nterms,-nterms:nterms)
                integer nterms
                real*8 ztarg(3,nvcount)
                complex*16 pot(nvcount)
                integer iffld
                complex*16 fld(3,nvcount)
                integer ier(nvcount)
                """, ["ier", "pot", "fld"])

    for what, extra_args in [
            ("l", ""),
            ("h", "complex*16 zk")
            ]:
        for expn_type in ["ta", "mp"]:
            gen_vector_wrapper("%s2d%seval" % (what, expn_type), """
                    %s
                    real*8 rscale
                    real*8 center(2)
                    complex*16 expn(-nterms:nterms)
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
        taeval_func_name = "%s3dtaeval" % what
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
                complex *16 expn(0:nterms,-nterms:nterms,nvcount)
                integer nterms
                real*8 ztarg(3,nvcount)
                complex*16 pot(nvcount)
                integer ifgrad
                complex*16 grad(3,nvcount)
                %s
                integer ier
                """ % (extra_args, hess_output), taeval_out_args,
                vec_func_name=taeval_func_name+"_1tgtperexp")

    # }}}

    # {{{ translation operators

    for dims in [2, 3]:
        for eqn in [cgh.Laplace(dims), cgh.Helmholtz(dims)]:
            if eqn.lh_letter() == "l" and dims == 3:
                continue

            for xlat in ["mpmp", "mploc", "locloc"]:
                func_name = "%s%dd%s" % (eqn.lh_letter(), dims, xlat)

                if eqn.lh_letter() == "h" and dims == 3:
                    func_name += "quadu"

                args_template = Template("""
                    ${ extra_args }

                    real*8 rscale1(${input_dim})
                    real*8 center1(${dims}, ${input_dim})
                    complex*16 expn1(${expn_dims_1}, ${input_dim})
                    integer nterms1

                    real*8 rscale2(nvcount)
                    real*8 center2(${dims}, nvcount)
                    complex*16 expn2(${expn_dims_2}, nvcount)
                    integer nterms2

                    %if lh_letter == "h" and dims == 3:
                        real*8 radius(nvcount)
                        real*8 xnodes(nquad)
                        real*8 wts(nquad)
                        integer nquad
                        integer ier(nvcount)
                    %endif

                    """, strict_undefined=True)

                for (
                        vec_func_name,
                        input_dim,
                        output_reductions,
                        tmp_init,
                        ) in [
                        (func_name + "_vec", "nvcount",
                            None, None),
                        (func_name + "_imany", "*INDIRECT_MANY",
                            {"expn2": "sum", "ier": "max"},
                            {"ier": "0"}),
                        ]:
                    args = args_template.render(
                        lh_letter=eqn.lh_letter(),
                        dims=dims,
                        expn_dims_1=eqn.expansion_dims("nterms1"),
                        expn_dims_2=eqn.expansion_dims("nterms2"),
                        extra_args=eqn.in_arg_decls(with_intent=False),
                        input_dim=input_dim,
                        )
                    gen_vector_wrapper(func_name, args, ["ier", "expn2"],
                            output_reductions=output_reductions,
                            tmp_init=tmp_init,
                            vec_func_name=vec_func_name)

    # }}}

    result.append("! vim: filetype=fortran")

    return "\n".join(result)


if __name__ == "__main__":
    print(gen_vector_wrappers())

# vim: foldmethod=marker
