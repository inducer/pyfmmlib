! -*- f90 -*-
! If you're editing this, make sure the file extension is .pyf.mako, not .pyf.
! If it is .pyf, you are editing a generated file, and your changes will be
! overwritten.

python module _internal
  interface
      subroutine jfuns2d(ier,nterms,z,scale,fjs,ifder,fjder, &
          lwfjs,iscale,ntop)
        ! implicit real *8 (a-h,o-z)
        integer, intent(out) :: ier
        integer, intent(in) :: nterms
        complex*16, intent(in) :: z
        real*8, intent(in) :: scale
        integer, intent(in),check(lwfjs>=nterms+2) :: lwfjs
        complex*16, intent(out) :: fjs(0:lwfjs)
        integer, intent(in) :: ifder
        complex*16, intent(out) :: fjder(0:lwfjs)
        complex*16, intent(cache,hide) :: iscale(0:lwfjs)
        complex*16, intent(out) :: ntop
    end subroutine

    ! {{{ formta entrypoints

    % for dp_or_no in ["", "_dp"]:
      % for what in ["l", "h"]:
        % for dims in [2, 3]:
          subroutine ${what}${dims}dformta${dp_or_no}(ier, &
            % if what == "h":
              zk, &
            %endif
            rscale,source, &
            % if dp_or_no and not (what=="l" and dims == 2):
              dipstr,  dipvec, &
            % else:
              charge, &
            % endif
            ns, center, &
            nterms,locexp)
              intent(in) rscale,sources,charge,ns,center,nterms
              intent(out) ier,locexp

              implicit real *8 (a-h,o-z)
              complex *16 zk,charge(ns)
              dimension center(${dims}),source(${dims},ns),zdiff(${dims})
              real *8 dipvec(${dims},ns)
              complex *16 dipstr(ns)
              % if dims == 2:
                % if what == "l":
                  complex*16 locexp(0:nterms)
                % else:
                  complex*16 locexp(-nterms:nterms)
                % endif
              % else:
                complex*16 locexp(0:nterms,-nterms:nterms)
              % endif
          end subroutine
        %endfor
      % endfor
    % endfor

    ! }}}

    subroutine legeexps(itype,n,x,u,v,whts)
      integer, intent(in) :: itype, n
      real*8, intent(out) :: x(n),whts(n),u(n,n),v(n,n)
    end subroutine

    subroutine legefder(x,val,der,pexp,n)
      implicit real *8 (a-h,o-z)
      real *8 pexp(n+1)
      intent(out) val, der
      intent(in) x, pexp, n
    end subroutine

    ! {{{ fmm entrypoints

    <%! import codegen_helpers as cgh %>

    % for dim in [2, 3]:
    % for eqn in [cgh.Laplace(dim), cgh.Helmholtz(dim)]:

    <%

    if dim == 3:
        kinds = [
          "part",
          #"tria"
          ]
    else:
        kinds = ["part"]

    %>

    % for kind in kinds:

    <%
    if kind == "tria":
        iter_variants = ["", "iter"]
    else:
        iter_variants = [""]
    %>

    % for iter in iter_variants:

        <%

        have_dipvec = not (eqn.lh_letter() == "l" and dim == 2)

        suffix = "targ"

        has_hess = not (dim == 3 and eqn.lh_letter() == "h") and not iter

        has_hess_suffix = has_hess and (dim == 3)

        if has_hess_suffix:
            suffix = "hess"+suffix

        name_kind = kind

        %>

        subroutine ${eqn.lh_letter()}fmm${dim}d${name_kind}${iter}${suffix}( &
                ier, iprec, ${eqn.in_arg_list()|cpost} &
                nsource, &
                % if kind == "part":
                    source, &
                % elif kind == "tria":
                    triaflat, trianorm, source, &
                % endif
                ifcharge, charge, ifdipole, dipstr, &
                % if have_dipvec:
                    dipvec, &
                % endif
                ifpot, pot, iffld, fld, &
                % if has_hess:
                    ifhess, hess, &
                % endif
                ntarget, target, &
                ifpottarg, pottarg, iffldtarg, fldtarg &
                % if has_hess:
                    , ifhesstarg, hesstarg &
                % endif
                % if iter:
                    , icomp_type,wsave,lwsave,lused &
                % endif
                )
            implicit none

            integer, intent(out) :: ier
            integer, intent(in) :: iprec

            ${eqn.in_arg_decls()}

            integer, intent(in) :: nsource
            % if kind == "tria":
                real*8, intent(in) :: triaflat(3, 3, nsource)
                real*8, intent(in) :: trianorm(3, nsource)
            % endif
            real*8, intent(in) :: source(${dim},nsource)

            integer, intent(in) :: ifcharge
            complex*16, intent(in) :: charge(nsource)

            integer, intent(in) :: ifdipole
            complex*16, intent(in) :: dipstr(nsource)
            % if have_dipvec:
              real*8, intent(in) :: dipvec(${dim},nsource)
            % endif

            integer, intent(in) :: ifpot
            complex*16, intent(out) :: pot(nsource)
            integer, intent(in) :: iffld
            complex*16, intent(out) :: fld(${dim},nsource)
            % if has_hess:
                integer, intent(in) :: ifhess
                complex*16, intent(out) :: hess(${cgh.hess_size(dim)},nsource)
            % endif

            integer, intent(in) :: ntarget

            integer, intent(in) :: ifpottarg, iffldtarg
            real*8, intent(in) :: target(${dim}, ntarget)
            complex*16, intent(in,out) :: pottarg(ntarget)
            complex*16, intent(in,out) :: fldtarg(${dim},ntarget)

            % if has_hess:
                integer, intent(in) :: ifhesstarg
                complex*16, intent(in,out) :: hesstarg(${cgh.hess_size(dim)},ntarget)
            % endif

            % if iter:
                integer, intent(in) :: icomp_type,lwsave
                real*8, intent(in,out) :: wsave(*)
                integer, intent(out) :: lused
            % endif

            required ntarget

            check((!ifpottarg && !iffldtarg) || (shape(target,0)==${dim} && shape(target,1) == ntarget))  target
            check((!ifpottarg) || (shape(pottarg,0)==ntarget))  pottarg
            check((!iffldtarg) || (shape(fldtarg,0)==${dim} && shape(fldtarg,1) == ntarget))  fldtarg
            check((!ifhesstarg) || (shape(hesstarg,0)==${cgh.hess_size(dim)} && shape(hesstarg,1) == ntarget))  hesstarg

            check(!ifcharge || (shape(charge,0) == nsource))  charge
            depend(nsource)  charge
            check(!ifdipole || (shape(dipstr,0) == nsource))  dipstr
            depend(nsource)  dipstr

            ! F2PY workaround: pottarg, fldtarg must be input because f2py
            ! refuses to allocate zero-size output arrays.
            !
            ! This also means that these arrays might end up being 1 long
            ! even if ntarget == 0--but that is only allowed if the
            ! corresponding if*targ flag is off.

        end subroutine

    % endfor
    % endfor
    % endfor
    % endfor

    ! }}}

    subroutine l3dtaevalhess(rscale,center,mpole,nterms,ztarg, &
                        pot,iffld,fld,ifhess,hess,ier)
        implicit none
        integer, intent(in) :: ier,nterms,iffld,ifhess
        real *8, intent(in) :: rscale, center(3),ztarg(3)
        complex *16, intent(in) :: mpole(0:nterms,-nterms:nterms)
        complex *16, intent(out) :: pot,fld(3),hess(6)
        integer, intent(out) :: ier
    end subroutine

    ! }}}

    ! {{{ generated vectorized wrappers

    ${gen_vector_wrappers()}

    ! }}}

  end interface
end python module

! vim: filetype=fortran:foldmethod=marker
