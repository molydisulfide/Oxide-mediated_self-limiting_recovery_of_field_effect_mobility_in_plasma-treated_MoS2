from cvxopt import matrix
from cvxopt.lapack import *
from cvxopt.blas import *
from math import sqrt

# 2009_05_11 - removed unnecessary declaration of '_y'

"""
function [ x, istop, itn, anorm, acond, rnorm, xnorm ] = ...
           symmlq( n, b, aprodname, msolvename, iw, rw,  ...
                   precon, shift, show, check, maxit, rtol )
%        [ x, istop, itn, anorm, acond, rnorm, xnorm ] = ...
%          symmlq( n, b, aprodname, msolvename, iw, rw,  ...
%                  precon, shift, show, check, maxit, rtol )
%
%  SYMMLQ is designed to solve the system of linear equations Ax = b
%  where A is an n by n symmetric matrix and b is a given vector.
%  A is accessed by means of a function call of the form
%               y = aprod ( n, x, iw, rw )
%  which must return the product y = Ax for any given vector x.
%  A positive-definite preconditioner M = C C' may optionally
%  be specified.  If precon is true, a function call of the form
%               y = msolve( n, x, iw, rw )
%  must solve the system My = x.
%  WARNING:   The files containing the functions 'aprod' and 'msolve'
%             must not be called aprodname.m or msolvename.m !!!!
%
%  For further information, type    help symdoc.

%  07 Jun 1989: Date of Fortran 77 version written by
%               Michael Saunders, Stanford University.
%  15 May 1990: MATLAB m-file symmlq.m derived from Fortran version
%               by Gerald N. Miranda Jr, UCSD.
%  02 Oct 1997: Move to CG point only if it is better than LQ point.
%               For interest, print qrnorm (= rnorm for minres).
%               Note that cgnorm is always one step ahead of qrnorm.
%  20 Oct 1999: Bug.  alfa1 = 0 caused Anorm = 0, divide by zero.
%               Need to estimate Anorm from column of Tk.
%  ------------------------------------------------------------------
%
%  Initialize                               
"""
def symmlq(b,A,maxit=None,n=None,rtol=1e-9,M=None, show=False, shift=0.0, Gmat=None,eps=2.2e-16):
    """
    Set parameters
    """
    if n is None:     n     = len(b)
    if maxit is None: maxit = 2*n
    if M is None:        precon=False
    else:                precon=True


    first = 'Enter SYMMLQ.   ';
    last  = 'Exit  SYMMLQ.   ';
    space = ' ';
    msg={
       -1:' beta2 = 0.  If M = I, b and x are eigenvectors    ',
        0:' beta1 = 0.  The exact solution is  x = 0          ',
        1:' Requested accuracy achieved, as determined by rtol ',
        2:' Reasonable accuracy achieved, given eps           ',
        3:' x has converged to an eigenvector                 ',
        4:' acond has exceeded 0.1/eps                        ',
        5:' The iteration limit was reached                   ',
        6:' aprod  does not define a symmetric matrix         ',
        7:' msolve does not define a symmetric matrix         ',
        8:' msolve does not define a pos-def preconditioner   '}
    if show:
        print ""
        print first+ 'Solution of symmetric Ax = b'
        print 'n      =  %3g     precon =  %4g           shift  =  %23.14e' % (n,precon,shift )
        print 'maxit =  %3g     eps    =  %11.2e    rtol   =  %11.2e\n' % (maxit,eps,rtol )


    
    istop  = 0;   ynorm  = 0;    w = matrix(0., (n,1));  acond = 0;
    itn    = 0;   xnorm  = 0;    x = matrix(0., (n,1)); 
    anorm  = 0;   rnorm  = 0;    v = matrix(0., (n,1));
    done=False


    y      =+ b;     r1     = +b;
    if precon: M(y)
    b1     = y[0];  beta1  = dotu(r1, y);

    # no checking of M

    if beta1 <  0: istop = 8;  show = True;  done = True;
    if beta1 == 0:             show = True;  done = True;

    if beta1 > 0:
        beta1  = sqrt( beta1 );
        s      = 1 / beta1;
        #v      = s * y;
        copy(y,v)
        scal(s,v)
        
        A(v,y)
        #y    = (- shift) * v + y;
        axpy(v,y,-shift)
        alfa = dotu(v, y);
        #y    = (- alfa / beta1) * r1 + y;
        axpy(r1,y,-alfa/beta1)

        z  = dotu(v, y);
        s  = dotu(v, v);
        #y  = (- z / s) * v + y;
        axpy(v,y,-z/s)
        r2 = +y;
        if precon: M(y)
        oldb   = beta1;
        beta   = dotu(r2, y);
        if beta < 0: istop = 8; show = True;  done = True;
        #  Cause termination (later) if beta is essentially zero.

        beta  = sqrt( beta );
        if beta <= eps: istop = -1;

        #  See if the local reorthogonalization achieved anything.

        denom = sqrt( s ) * nrm2( r2 )  +  eps;
        s     = z / denom;
        t     = dotu(v, r2);
        t     = t / denom;

        if show:
            print ""
            print 'beta1 =  %10.2e   alpha1 =  %9.2e'% (beta1, alfa );
            print '(v1, v2) before and after  %14.2e'% s ;
            print 'local reorthogonalization  %14.2e'% t;
 

        #  Initialize other quantities.
        cgnorm = beta1;     rhs2   = 0;       tnorm  = alfa**2 + beta**2;
        gbar   = alfa;      bstep  = 0;       ynorm2 = 0;
        dbar   = beta;      snprod = 1;       gmax   = abs( alfa ) + eps;
        rhs1   = beta1;     x1cg   = 0;       gmin   = gmax;
        qrnorm = beta1;
    # end  % of beta1 > 0

    if show:
        head1 = '   Itn     x(1)(cg)  normr(cg)  r(minres)';
        head2 = '    bstep    anorm    acond';
        print head1+head2

        str1 =  '%6g %12.5e %10.3e'% (itn, x1cg, cgnorm);
        str2 =  ' %10.3e  %8.1e' %    (qrnorm, bstep/beta1 );
        print str1+str2

    
    """
      ------------------------------------------------------------------
      Main iteration loop.
      ------------------------------------------------------------------
      Estimate various norms and test for convergence.
    """
    
    if not done:
        while itn < maxit:
            itn    = itn  +  1;
            anorm  = sqrt( tnorm  );
            ynorm  = sqrt( ynorm2 );
            epsa   = anorm * eps;
            epsx   = anorm * ynorm * eps;
            epsr   = anorm * ynorm * rtol;
            diag   = gbar;

            if diag == 0: diag = epsa; 

            lqnorm = sqrt( rhs1**2 + rhs2**2 );
            qrnorm = snprod * beta1;
            cgnorm = qrnorm * beta / abs( diag );


            """
            %     Estimate  Cond(A).
            %     In this version we look at the diagonals of  L  in the
            %     factorization of the tridiagonal matrix,  T = L*Q.
            %     Sometimes, T(k) can be misleadingly ill-conditioned when
            %     T(k+1) is not, so we must be careful not to overestimate acond.
            """
            if lqnorm < cgnorm:
                acond  = gmax / gmin;
            else:
                denom  = min( gmin, abs( diag ) );
                acond  = gmax / denom;

            zbar   = rhs1 / diag;
            z      = (snprod * zbar + bstep) / beta1;
            x1lq   = x[0] + b1 * bstep / beta1;
            x1cg   = x[0] + w[0] * zbar  +  b1 * z;

            """
            %     See if any of the stopping criteria are satisfied.
            %     In rare cases, istop is already -1 from above (Abar = const * I).
            """
            if istop == 0:
                if itn    >= maxit : istop = 5; 
                if acond  >= 0.1/eps: istop = 4; 
                if epsx   >= beta1  : istop = 3; 
                if cgnorm <= epsx   : istop = 2; 
                if cgnorm <= epsr   : istop = 1; 
                                    
                                    
            prnt = 0;               
            if n      <= 40         :   prnt = 1;
            if itn    <= 10         :   prnt = 1;
            if itn    >= maxit - 10 :   prnt = 1;
            if itn%10 == 0          :   prnt = 1;
            if cgnorm <= 10.0*epsx  :   prnt = 1;
            if cgnorm <= 10.0*epsr  :   prnt = 1;
            if acond  >= 0.01/eps   :   prnt = 1;
            if istop  != 0          :   prnt = 1;


            if show and prnt == 1:
                str1 =  '%6g %12.5e %10.3e'% (itn, x1cg, cgnorm);
                str2 =  ' %10.3e  %8.1e'%    (qrnorm, bstep/beta1);
                str3 =  ' %8.1e %8.1e'%      (anorm, acond );
                print str1+str2+str3

            if istop !=0: break

            s      = 1/beta;
        
            copy(y,v)
            scal(s,v)
            #v      = s * y;
            #y      = feval( aprodname, n, v, iw, rw );
            A(v,y)
            #y      = (- shift) * v + y;
            axpy(v,y,-shift)
            #y      = (- beta / oldb) * r1 + y;
            axpy(r1,y,-beta/oldb)
            alfa   = dotu(v, y);
            #y      = (- alfa / beta) * r2 + y;
            axpy(r2,y,-alfa/beta)
            ###r1     = r2;
            ###r2     = y;
            copy(y,r1)
            _y=r1
            r1=r2
            r2=y
            y=_y
            if precon: M(y)
            oldb   = beta;
            beta   = dotu(r2, y);
            
            if beta < 0: istop = 6;  break
            beta   = sqrt( beta );
            tnorm  = tnorm  +  alfa**2  +  oldb**2  +  beta**2;

            """
                 Compute the next plane rotation for Q.
            """
            gamma  = sqrt( gbar**2 + oldb**2 );
            cs     = gbar / gamma;
            sn     = oldb / gamma;
            delta  = cs * dbar  +  sn * alfa;
            gbar   = sn * dbar  -  cs * alfa;
            epsln  = sn * beta;
            dbar   =            -  cs * beta;
            """
            %     Update  X.
            """
            z      = rhs1 / gamma;
            s      = z*cs;
            t      = z*sn;

            axpy(v,x,t)
            axpy(w,x,s)
            #w      = sn*w - cs*v;
            scal(sn,w)
            axpy(v,w,-cs)
            """
            %     Accumulate the step along the direction  b, and go round again.
            """
            bstep  = snprod * cs * z  +  bstep;
            snprod = snprod * sn;
            gmax   = max( gmax, gamma );
            gmin   = min( gmin, gamma );
            ynorm2 = z**2  +  ynorm2;
            rhs1   = rhs2  -  delta * z;
            rhs2   =       -  epsln * z;
        # end % while
    """
    %  ------------------------------------------------------------------
    %  End of main iteration loop.
    %  ------------------------------------------------------------------
    
    %  Move to the CG point if it seems better.
    %  In this version of SYMMLQ, the convergence tests involve
    %  only cgnorm, so we're unlikely to stop at an LQ point,
    %  EXCEPT if the iteration limit interferes.
    """
    if cgnorm < lqnorm:
        zbar   = rhs1 / diag;
        bstep  = snprod * zbar  +  bstep;
        ynorm  = sqrt( ynorm2  +  zbar**2 );
        # x      = zbar * w + x;
        axpy(w,x,zbar)

    """
    %  Add the step along  b.
    """
    bstep  = bstep / beta1;
    #y      = +b;
    copy(b,y)
    if precon:  M(y) #y =  feval( msolvename, n, b, iw, rw );
    #x      = bstep * y + x;
    axpy(y,x,bstep)
    """
    %  Compute the final residual,  r1 = b - (A - shift*I)*x.
    """
    A(x,y)
    #y      = (- shift) * x + y;
    axpy(x,y,-shift)
    r1     = b - y;
    rnorm  = nrm2 ( r1 );
    xnorm  = nrm2 (  x );
    """
    %  ==================================================================
    %  Display final status.
    %  ==================================================================
    """
    if show:
        print ""
        print last + ' istop   =  %3g               itn   =  %5g' % (istop, itn)
        print last + ' anorm   =  %12.4e      acond =  %12.4e' % (anorm, acond)
        print last + ' rnorm   =  %12.4e      xnorm =  %12.4e' % (rnorm, xnorm)
        print last +  msg[istop]

    return x, istop, itn, anorm, acond, rnorm, xnorm
