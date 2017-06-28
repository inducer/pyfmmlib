from __future__ import division


def hess_size(d):
    return (d*d+d)//2


# {{{ equations

class Equation(object):
    def __init__(self, dimensions):
        self.dimensions = dimensions

    def lh_letter(self):
        raise NotImplementedError

    def k_arg(self):
        raise NotImplementedError

    local_expansion_type = "complex*16"


class Laplace(Equation):
    def in_arg_list(self):
        return ""

    def in_arg_decls(self, with_intent=True):
        return ""

    def expansion_dims(self, nterms="nterms"):
        if self.dimensions == 2:
            return "0:%(nterms)s" % locals()
        elif self.dimensions == 3:
            return "0:%(nterms)s,-(%(nterms)s):%(nterms)s" % locals()
        else:
            raise RuntimeError("invalid dimensionality")

    def local_eval_cache_dims(self, nterms="nterms"):
        if self.dimensions == 2:
            return "0:%(nterms)s" % locals()
        else:
            raise NotImplementedError("unimplemented dimensionality")

    def lh_letter(self):
        return "l"

    def k_arg(self):
        return ""


class Helmholtz(Equation):
    def in_arg_list(self):
        return "zk"

    def in_arg_decls(self, with_intent=True):
        if with_intent:
            return """
                complex*16, intent(in) :: zk
                """
        else:
            return """
                complex*16 zk
                """

    def expansion_dims(self, nterms="nterms"):
        if self.dimensions == 2:
            return "-(%(nterms)s):%(nterms)s" % locals()
        elif self.dimensions == 3:
            return "0:%(nterms)s,-(%(nterms)s):%(nterms)s" % locals()
        else:
            raise RuntimeError("invalid dimensionality")

    def local_eval_cache_dims(self, nterms="nterms"):
        if self.dimensions == 2:
            return "-(%(nterms)s)-2:(%(nterms)s)+2" % locals()
        else:
            raise NotImplementedError("unimplemented dimensionality")

    def lh_letter(self):
        return "h"

    def k_arg(self):
        return "zk, "

# }}}


class SumComputation:
    def __init__(self, eqn, outputs, source_derivative=False):
        """
        :param source_derivative: a :class:`bool`.
        """

        self.eqn = eqn
        self.outputs = outputs
        self.source_derivative = source_derivative

        if source_derivative and isinstance(eqn, Laplace) and eqn.dimensions == 2:
            raise NotImplementedError("2D double-layer Laplace")

    def eqn_and_dim(self):
        return "%s%dd" % (
                self.eqn.lh_letter(),
                self.dimensions)

    def what(self):
        if self.source_derivative:
            src_deriv = "_sd"
        else:
            src_deriv = ""

        return "%s%s%s_%dd" % (
                self.eqn.lh_letter(),
                self.outputs, src_deriv,
                self.dimensions)

    @property
    def dimensions(self):
        return self.eqn.dimensions

    def output_arg_names_and_dims(self, force_pot=False):
        saw_pot = False

        for o in self.outputs:
            if o == "p":
                saw_pot = True
                yield "pot", ()
            elif o == "f":
                yield "fld", (self.dimensions,)
            elif o == "g":
                yield "grad", (self.dimensions,)
            elif o == "h":
                yield "hess", (hess_size(self.dimensions),)
            else:
                raise ValueError("unknown output name")

        if force_pot and not saw_pot:
            yield "pot", ()

    def get_output_kwargs(self, suffix="", force_pot=False):
        return dict(
                (name, name+suffix)
                for name, dims in
                self.output_arg_names_and_dims(force_pot=force_pot))
