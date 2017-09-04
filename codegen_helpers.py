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

    def expansion_dims(self, nterms):
        if self.dimensions == 2:
            return "0:%(nterms)s" % locals()
        elif self.dimensions == 3:
            return "0:%(nterms)s,-(%(nterms)s):%(nterms)s" % locals()
        else:
            raise RuntimeError("invalid dimensionality")

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

    def expansion_dims(self, nterms):
        if self.dimensions == 2:
            return "-(%(nterms)s):%(nterms)s" % locals()
        elif self.dimensions == 3:
            return "0:%(nterms)s,-(%(nterms)s):%(nterms)s" % locals()
        else:
            raise RuntimeError("invalid dimensionality")

    def lh_letter(self):
        return "h"

    def k_arg(self):
        return "zk, "

# }}}


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
