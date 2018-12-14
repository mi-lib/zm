/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nle_se - nonlinear equation: single nonlinear equation solver.
 */

#include <zm/zm_nle.h>

#define _ZNLE_TEST(__xval,__fval,__f,__priv,__tol,__ret) \
  if( zIsTol( ( __fval = __f( __xval, __priv ) ), __tol ) ) return (__ret)

#define ZNLE_ITERWARN( iter ) do{\
  ZITERWARN( iter );\
  ZRUNWARN( ZM_WARN_NLE_NOROOT );\
} while(0)

static int _zNLE_Bracket(double (* f)(double,void*), double *x1, double *x2, void *priv, double *f1, double *f2, double tol, int iter);

/* _zNLE_Bracket
 * - find initial bracketing values of function.
 *   it returns 0 for success case, or -1 for failure case.
 *   (static)
 */
int _zNLE_Bracket(double (* f)(double,void*), double *x1, double *x2, void *priv, double *f1, double *f2, double tol, int iter)
{
  long count = 0;
  double _x1, _x2;

  zRandInit();
  for( _x1=*x1, _x2=*x2; count < iter;
       _x1=zRandF(*x1,*x2), _x2=zRandF(*x1,*x2), count++ ){
    _ZNLE_TEST( _x1, *f1, f, priv, tol, 1 );
    _ZNLE_TEST( _x2, *f2, f, priv, tol, 2 );
    if( *f1 * *f2 < 0 ){ /* success */
      *x1 = _x1;
      *x2 = _x2;
      return 0;
    }
  }
  ZNLE_ITERWARN( iter );
  return -1; /* failure */
}

/* zNLE_Bisec
 * - bisection method.
 */
double zNLE_Bisec(double (* f)(double,void*), double x1, double x2, void *priv, double tol, int iter)
{
  double xm;
  double f1, f2, fm;
  register int i = 0;

  ZITERINIT( iter );
  switch( _zNLE_Bracket( f, &x1, &x2, priv, &f1, &f2, tol, iter ) ){
  case -1: return 0;
  case 1:  return x1;
  case 2:  return x2;
  default: ;
  }
  ZITERINIT( iter );
  do{
    xm = 0.5 * ( x1 + x2 );
    _ZNLE_TEST( xm, fm, f, priv, tol, xm );
    if( f1*fm < 0 ){
      x2 = xm; f2 = fm;
    } else{
      x1 = xm; f1 = fm;
    }
  } while( i++ < iter );
  ZNLE_ITERWARN( iter );
  return xm;
}

/* zNLE_Secant
 * - Secant method.
 */
#define ZNLE_D 0.1
double zNLE_Secant(double (* f)(double,void*), double x1, double x2, void *priv, double tol, int iter)
{
  double f1, f2, fc, xc = 0 /* dummy */;
  long count;

  zRandInit();
  ZITERINIT( iter );
  for( count=0; ; x1+=zRandF(-ZNLE_D,ZNLE_D) ){
    _ZNLE_TEST( x1, f1, f, priv, tol, x1 );
    _ZNLE_TEST( x2, f2, f, priv, tol, x2 );
    if( !zIsTol( f1 - f2, tol ) ) break;
    if( ++count > iter ){
      ZNLE_ITERWARN( iter );
      return 0;
    }
  }
  for( ; count<iter; count++ ){
    xc = ( f2*x1 - f1*x2 ) / ( f2 - f1 );
    _ZNLE_TEST( xc, fc, f, priv, tol, xc );
    x1 = x2; f1 = f2;
    x2 = xc; f2 = fc;
  }
  ZNLE_ITERWARN( iter );
  return xc;
}

/* zNLE_RF
 * - Regula-Falsi method.
 */
double zNLE_RF(double (* f)(double,void*), double x1, double x2, void *priv, double tol, int iter)
{
  double f1, f2, fc, xc = 0 /* dummy */;
  long count;

  ZITERINIT( iter );
  switch( _zNLE_Bracket( f, &x1, &x2, priv, &f1, &f2, tol, iter ) ){
  case -1: return 0;
  case 1:  return x1;
  case 2:  return x2;
  default: ;
  }
  for( count=0; count<iter; count++ ){
    xc = ( f2*x1 - f1*x2 ) / ( f2 - f1 );
    _ZNLE_TEST( xc, fc, f, priv, tol, xc );
    if( f1*fc < 0 ){
      x2 = xc; f2 = fc;
    } else{
      x1 = xc; f1 = fc;
    }
  }
  ZNLE_ITERWARN( iter );
  return xc;
}

/* zNLE_VDB
 * - Van Wijngaarden-Dekker-Brent method
 */
double zNLE_VDB(double (* f)(double,void*), double x1, double x2, void *priv, double tol, int iter)
{
#define ZNLE_VDB_TOL zTOL
  long count;
  double x3, dx1, dx2, min1, min2;
  double f1, f2, f3, p, q, r, s, tol2, xm;

  ZITERINIT( iter );
  switch( _zNLE_Bracket( f, &x1, &x2, priv, &f1, &f2, tol, iter ) ){
  case -1: return 0;
  case 1:  return x1;
  case 2:  return x2;
  default: ;
  }
  x3 = x2;
  f3 = f2;
  for( count=0; count<iter; count++ ){
    if( f2*f3 > 0 ){
      x3 = x1;
      f3 = f1;
      dx1 = dx2 = x2 - x1;
    }
    if( fabs(f3) < fabs(f2) ){
      x1 = x2; x2 = x3; x3 = x1;
      f1 = f2; f2 = f3; f3 = f1;
    }
    tol2 = 2*zTOL*fabs(x2) + 0.5*ZNLE_VDB_TOL;
    dx1 = dx2 = xm = 0.5*( x3 - x2 );
    if( fabs(xm) <= tol2 || zIsTol(f2,tol) ) return x2;
    if( fabs(dx2) >= tol2 && fabs(f1) > fabs(f2) ){
      s = f2 / f1;
      if( x1 == x3 ){
        p = 2 * xm * s;
        q = 1 - s;
      } else{
        q = f1 / f3;
        r = f2 / f3;
        p = s * ( s*xm*q*(q-r) - (x2-x1)*(r-1) );
        q = (q-1)*(r-1)*(s-1);
      }
      if( p > 0 )
        q = -q;
      else
        p = -p;
      min1 = 3*xm*q - fabs(tol2*q);
      min2 = fabs(dx2*q);
      if( 2*p < zMin(min1,min2) ){
        dx2 = dx1;
        dx1 = p / q;
      } else{
        dx1 = dx2 = xm;
      }
    } else{
      dx1 = dx2 = xm;
    }
    x1 = x2;
    f1 = f2;
    if( fabs(dx1) > tol2 )
      x2 += dx1;
    else
      x2 += zSgn(tol2)*xm;
    f2 = f( x2, priv );
  }
  ZNLE_ITERWARN( iter );
  return x3;
}
