/* % -*- mode: C; mode: fold -*-

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 Compute the Jacobi elliptic functions sn(u|m), cn(u|m) and dn(u|m) for 
 argument u (real or complex) and parameter m. 

 usage: [sn,cn,dn] = ellipj(u,m[,tol])
 
 u and can be complex.
 m is restricted to 0 <= m <= 1.
 They can be scalars, matrix and scalar, scalar and matrix,
 column and row, conformant matrices.

 modified so u can be complex.   Leopoldo Cerbaro redbliss@libero.it
 
 Ref: Abramowitz, Milton and Stegun, Irene A
      Handbook of Mathematical Functions, Dover, 1965
      Chapter 16 (Sections 16.4, 16.13 and 16.15)

 Based upon ellipj.m  made by David Billinghurst <David.Billinghurst@riotinto.com>
 and besselj.cc

 Author: Leopoldo Cerbaro <redbliss@libero.it>
 Created: 15 December 2001
*/

#include "oct.h"
#include "lo-ieee.h"  /* for octave_NaN */


static void
gripe_ellipj_arg ( const char *arg)
{
  error ("ellipj: expecting scalar or matrix as %s argument", arg);
}

const double  eps = 2.220446049e-16;
const int  Nmax = 16;

static void
sncndn ( double u, double m, double& sn, double& cn, double& dn, double& err) {
/* real */
double sqrt_eps, m1, t=0., si_u, co_u, se_u, ta_u, b, c[Nmax], a[Nmax], phi;
int n, Nn, ii;

  if (m < 0. || m > 1.) {
     warning ("ellipj: expecting 0. <= m <= 1."); /* -lc- */
     sn = cn = dn = octave_NaN;
     return;
	}
  sqrt_eps = sqrt(eps);
  if (m < sqrt_eps) {
    /*  # For small m, ( Abramowitz and Stegun, Section 16.13 ) */
    /*{{{*/
        si_u = sin(u);
        co_u = cos(u);
        t = 0.25*m*(u-si_u*co_u);
        sn = si_u - t * co_u;
        cn = co_u + t * si_u;
        dn = 1.0 - 0.5*m*si_u*si_u;
/*}}}*/
  } else if ( (1.0 - m) < sqrt_eps ) {
    /*  For m1 = (1-m) small ( Abramowitz and Stegun, Section 16.15 ) */
    /*{{{*/
        m1 = 1.0-m;
        si_u = sinh(u);
        co_u = cosh(u);
        ta_u = tanh(u);
        se_u = 1.0/co_u;
        sn = ta_u + 0.25*m1*(si_u*co_u-u)*se_u*se_u;
        cn = se_u - 0.25*m1*(si_u*co_u-u)*ta_u*se_u;
        dn = se_u + 0.25*m1*(si_u*co_u+u)*ta_u*se_u;
/*}}}*/
  } else {
    /*{{{*/
        /*
        //  Arithmetic-Geometric Mean (AGM) algorithm
        //    ( Abramowitz and Stegun, Section 16.4 )
        */
       
        a[0] = 1.0;
        b    = sqrt(1.0-m);
        c[0] = sqrt(m);
        for (n = 1; n<Nmax; ++n) {
          a[n] = (a[n-1]+b)/2;
          c[n] = (a[n-1]-b)/2;
          b = sqrt(a[n-1]*b);
          if ( c[n]/a[n] < eps) break; 
				}
        if ( n >= Nmax-1) {
           // fprintf(stderr, "Not enough workspace\n");
           err = 1.;
           return;
        }
        Nn = n;
        for ( ii = 1;  n>0;	ii = ii*2, --n) ; // pow(2, Nn)
        phi = ii*a[Nn]*u;
        for ( n = Nn; n > 0; --n) {
          t = phi;
          phi = (asin((c[n]/a[n])* sin(phi))+phi)/2.;
        }
        sn = sin(phi);
        cn = cos(phi);
        dn = cn/cos(t-phi);
/*}}}*/
  }
 return;
}

static void
sncndn ( Complex& u, double m, 
         Complex& sn, Complex& cn, Complex& dn, double& err) {
double m1 = 1.-m, ss1, cc1, dd1;

  sncndn( imag(u), m1, ss1, cc1, dd1, err);
  if ( real(u) == 0.) { 
    /* u is pure imag: Jacoby imag. transf. */
    /*{{{*/
    sn = Complex (0. , ss1/cc1);
    cn = 1/cc1;         //    cn.imag = 0.;
    dn = dd1/cc1;       //    dn.imag = 0.;
    /*}}}*/
  } else {
    /* u is generic complex */
    /*{{{*/
    double ss, cc, dd, ddd;

    sncndn( real(u), m, ss, cc, dd, err);
      ddd = cc1*cc1 + m*ss*ss*ss1*ss1;
      sn = Complex (ss*dd1/ddd, cc*dd*ss1*cc1/ddd); 
      cn = Complex (cc*cc1/ddd, -ss*dd*ss1*dd1/ddd);
      dn = Complex (dd*cc1*dd1/ddd, -m*ss*cc*ss1/ddd);
    /*}}}*/
  }
 return;
}

DEFUN_DLD (ellipj, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {[@var{sn}, @var{cn}, @var{dn}] =} \
ellipj (@var{u}, @var{m}, @var{err})\n\
Compute the Jacobi elliptic functions sn, cn, dn of complex argument and real parameter.\n\
\n\
If @var{m} is a scalar, the results are the same size as @var{u}.\n\
If @var{u} is a scalar, the results are the same size as @var{m}.\n\
If @var{u} is a column vector and @var{m} is a row vector, the\n\
results are matrices with @code{length (@var{u})} rows and\n\
@code{length (@var{m})} columns.  Otherwise, @var{u} and\n\
@var{m} must conform and the results will be the same size.\n\
\n\
The value of @var{u} may be complex.\n\
The value of @var{m} must be 0 <= m <= 1. .\n\
\n\
If requested, @var{err} contains the following status information\n\
and is the same size as the result.\n\
\n\
@enumerate 0\n\
@item\n\
Normal return.\n\
@item\n\
Error---no computation, algorithm termination condition not met,\n\
return @code{NaN}.\n\
@end enumerate\n\
@end deftypefn")
{
  octave_value_list retval;

  int nargin = args.length ();

  if (nargin == 2 ) {
      octave_value u_arg = args(0);
      octave_value m_arg = args(1);

      if (m_arg.is_scalar_type ()) {  // m is scalar
        double  m = args(1).double_value ();

        if (! error_state) {

          if (u_arg.is_scalar_type ()) {   /*  u scalar */
            /*{{{*/
            if (u_arg.is_real_type ()) {  // u real
              double  u = args(0).double_value ();

              if (! error_state) {
                double sn, cn, dn; 
                double err=0;
                octave_value result;

                sncndn(u, m, sn, cn, dn, err);
                retval (0) = sn;
                retval (1) = cn;
                retval (2) = dn;
                if (nargout > 3)
                  retval(3) =  err;
            } else 
                gripe_ellipj_arg ( "first");

            } else {  // u complex
              Complex u = u_arg.complex_value ();

              if (! error_state) {
                Complex sn, cn, dn;
                double err;
                octave_value result;

            		sncndn( u, m, sn, cn, dn, err);

                retval (0) = sn;
                retval (1) = cn;
                retval (2) = dn;
                if (nargout > 3)  retval(3) = err;
              } else
                gripe_ellipj_arg ( "second");
            }
            /*}}}*/
          } else {  /* u is matrix ( m is scalar ) */
            /*{{{*/
            ComplexMatrix u = u_arg.complex_matrix_value ();

            if (! error_state) {
              octave_value result;
              int nr = u.rows ();
              int nc = u.cols ();

              ComplexMatrix sn (nr, nc), cn (nr, nc), dn (nr, nc);
              Matrix err (nr, nc);

              for (int j = 0; j < nc; j++)
                for (int i = 0; i < nr; i++)
                  sncndn (u(i,j), m, sn(i,j), cn(i,j), dn(i,j), err(i,j));

                retval (0) = sn;
                retval (1) = cn;
                retval (2) = dn;
                if (nargout > 3)  retval(3) = err;
            } else
                gripe_ellipj_arg ( "first");
            /*}}}*/
          }
	      } else
            gripe_ellipj_arg ( "second");
     } else { // m is matrix
       Matrix m = args(1).matrix_value ();

       if (! error_state) {
         int mr = m.rows ();
         int mc = m.cols ();

         if (u_arg.is_scalar_type ()) {    /* u is scalar */
           /*{{{*/
           octave_value result;
           int nr = m.rows ();
           int nc = m.cols ();
           Matrix err (nr, nc);

           if (u_arg.is_real_type ()) {
             double  u = u_arg.double_value ();
             Matrix sn (nr, nc), cn (nr, nc), dn (nr, nc);
             if (! error_state) {
               for (int j = 0; j < nc; j++)
                 for (int i = 0; i < nr; i++)
                   sncndn (u, m(i,j), sn(i,j), cn(i,j), dn(i,j), err(i,j));

               retval (0) = sn;
               retval (1) = cn;
               retval (2) = dn;
               if (nargout > 3)  retval(3) = err;
             } else
               gripe_ellipj_arg ( "first");
           } else {
             Complex u = u_arg.complex_value ();
             ComplexMatrix sn (nr, nc), cn (nr, nc), dn (nr, nc);
             if (! error_state) {
               for (int j = 0; j < nc; j++)
                 for (int i = 0; i < nr; i++)
                   sncndn (u, m(i,j), sn(i,j), cn(i,j), dn(i,j), err(i,j));
               retval (0) = sn;
               retval (1) = cn;
               retval (2) = dn;
               if (nargout > 3)  retval(3) = err;
             } else
               gripe_ellipj_arg ( "first");
           }
           /*}}}*/
         } else {    // u is matrix  (m is matrix)
           /*{{{*/
           if (u_arg.is_real_type ()) {  // u real matrix

              Matrix u = u_arg.matrix_value ();
              if (! error_state) {
                int ur = u.rows ();
                int uc = u.cols ();

              if (mr == 1 && uc == 1)	{  // u column, m row
                RowVector rm = m.row (0);
                ColumnVector cu = u.column (0);

                Matrix sn (ur, mc), cn (ur, mc), dn (ur, mc);
                Matrix err(ur,mc);
//               octave_value result;

                for (int j = 0; j < mc; j++)
                  for (int i = 0; i < ur; i++)
                    sncndn (cu(i), rm(j), sn(i,j), cn(i,j), dn(i,j), err(i,j));

                retval (0) = sn;
                retval (1) = cn;
                retval (2) = dn;
                if (nargout > 3)  retval(3) = err;
              } else if (ur == mr && uc == mc)	{
                Matrix sn (ur, mc), cn (ur, mc), dn (ur, mc);
                Matrix err(ur,mc);
//               octave_value result;

                for (int j = 0; j < uc; j++)
                 for (int i = 0; i < ur; i++)
                  sncndn (u(i,j), m(i,j), sn(i,j), cn(i,j), dn(i,j), err(i,j));

                retval (0) = sn;
                retval (1) = cn;
                retval (2) = dn;
                if (nargout > 3)  retval(3) = err;
              } else
                 error("u m invalid");   
              } else
                gripe_ellipj_arg ( "first ");
            } else {  // u complex matrix
              ComplexMatrix u = u_arg.complex_matrix_value ();
              if (! error_state) {
                int ur = u.rows ();
                int uc = u.cols ();

              if (mr == 1 && uc == 1)	{
                RowVector rm = m.row (0);
                ComplexColumnVector cu = u.column (0);

                ComplexMatrix sn (ur, mc), cn (ur, mc), dn (ur, mc);
                Matrix err(ur,mc);
//               octave_value result;

                for (int j = 0; j < mc; j++)
                  for (int i = 0; i < ur; i++)
                    sncndn (cu(i), rm(j), sn(i,j), cn(i,j), dn(i,j), err(i,j));

                retval (0) = sn;
                retval (1) = cn;
                retval (2) = dn;
                if (nargout > 3)  retval(3) = err;
              } else if (ur == mr && uc == mc)	{

                ComplexMatrix sn (ur, mc), cn (ur, mc), dn (ur, mc);
                Matrix err(ur,mc);
//               octave_value result;

                for (int j = 0; j < uc; j++)
                 for (int i = 0; i < ur; i++)
                  sncndn (u(i,j), m(i,j), sn(i,j), cn(i,j), dn(i,j), err(i,j));

                retval (0) = sn;
                retval (1) = cn;
                retval (2) = dn;
                if (nargout > 3)  retval(3) = err;
              } else
                 error("u m invalid");   
              } else
                gripe_ellipj_arg ( "second");
            }
           /*}}}*/
         }
       } else
          gripe_ellipj_arg ( "second");
     }  // m matrix
   } else  // wrong n. of argin
       print_usage ("ellipj");
   return retval;
}



/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
