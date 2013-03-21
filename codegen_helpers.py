from __future__ import division

# {{{

def hess_size(d):
    return (d*d+d)//2

# }}}

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

    def form_local_expansion(self, center, locexp, sources, charge, nsources,
            nterms="nterms", rscale="rscale", zk="zk", ier_value=2,
            source_derivative_vec=None):
        lh_letter = self.eqn.lh_letter()
        k_arg = self.eqn.k_arg()

        preamble = ""
        suffix = ""

        if self.source_derivative:
            suffix = "_dp"

            if source_derivative_vec is None:
                raise ValueError("source derivative vector not available")
                direction = source_derivative_vec

            charge = charge + ", " + source_derivative_vec

        if self.dimensions == 2:
            routine_name = "%(lh_letter)s2dformta%(suffix)s" % locals()
            return """
                %(preamble)s

                call %(routine_name)s(ier, %(k_arg)s%(rscale)s, &
                    %(sources)s, &
                    %(charge)s, &
                    %(nsources)s, &
                    %(center)s, %(nterms)s, %(locexp)s)

                if (ier.ne.0) then
                    write(*,*) '%(routine_name)s failed'
                    ier=%(ier_value)s
                    return
                endif
                """ % locals()
        else:
            raise RuntimeError("invalid dimensionality")

    def eval_local_expansion(self,
            center, locexp, target,
            pot, grad="zdummy1", hess="zdummy1",
            nterms="nterms", zk="zk", rscale="rscale", 
            cache_valid=None, cache=None, residuals=None):
        ifgrad = int(grad != "zdummy1")
        ifhess = int(hess != "zdummy1")

        lh_letter = self.eqn.lh_letter()
        k_arg = self.eqn.k_arg()

        extra_args = ", ier"
        suffix = ""
        extra_stmts = ""
        if cache_valid is not None:
            suffix = "_cache"
            extra_args += ", %s, %s" % (cache_valid, cache)
            extra_stmts = "if (ier.ne.0) return"

        if residuals:
            suffix += "_estimate"
            extra_args += ", %s, %s" % residuals

        var_dict = locals()

        if self.dimensions == 2:
            result = """
                ifgrad = %(ifgrad)s
                ifhess = %(ifhess)s
                call nsq_%(lh_letter)s2dtaeval%(suffix)s(%(k_arg)s%(rscale)s, &
                    %(center)s, %(locexp)s, %(nterms)s, &
                    %(target)s,  &
                    %(pot)s, ifgrad, %(grad)s, ifhess, &
                    %(hess)s%(extra_args)s)
                %(extra_stmts)s
                """ % var_dict

            if isinstance(self, Laplace):
                result += """
                    %(pot)s = real(%(pot)s)
                    """ % var_dict
                if ifgrad:
                    result += """
                        %(grad)s = real(%(grad)s)
                        """ % var_dict
                if ifhess:
                    result += """
                        %(hess)s = real(%(hess)s)
                        """ % var_dict

            return result
        else:
            raise RuntimeError("invalid dimensionality")

    def compute_direct_sum(self,
            sources, charge, nsources, target,
            pot, grad="0", hess="0", zk="zk",
            source_derivative_vec=None):
        ifgrad = int(grad != "0")
        ifhess = int(hess != "0")

        lh_letter = self.eqn.lh_letter()
        k_arg = self.eqn.k_arg()

        suffix = ""
        if self.source_derivative:
            suffix = "_dp"

            if source_derivative_vec is None:
                raise ValueError("source_derivative_vec not available")

            charge = charge + ", " + source_derivative_vec

        if self.dimensions == 2:
            return """call %(lh_letter)spotgrad2dall%(suffix)s(%(ifgrad)s, %(ifhess)s, &
                %(sources)s, &
                %(charge)s, &
                %(nsources)s, &
                %(target)s, %(k_arg)s&
                %(pot)s, %(grad)s, %(hess)s)""" % locals()
        else:
            raise RuntimeError("invalid dimensionality")
