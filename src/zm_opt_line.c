/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_line - optimization tools: line search.
 */

#include <zm/zm_opt.h>

/* zOptLineGSEC
 * - golden section method.
 */
double zOptLineGSEC(double (*eval)(double,void*), double a, double b, void *util, int iter)
{
  register int i;
  double c, d, e1, e2, x;

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    x = 0.5 * ( a + b );
    if( zIsTiny( a - b ) ) return x;
    c = a - Z_GSEC * ( a - b );
    d = b + Z_GSEC * ( a - b );
    e1 = eval( c, util );
    e2 = eval( d, util );
    if( e1 >= e2 ) a = c;
    if( e1 <= e2 ) b = d;
  }
  ZITERWARN( iter );
  return x;
}

/* zOptLineBisec
 * - bisection method.
 */
double zOptLineBisec(double (*eval)(double,void*), double a, double b, void *util, int iter)
{
  register int i;
  double x, e1, e2, e3;

  e1 = eval( a, util );
  e2 = eval( b, util );
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    x = 0.5 * ( a + b );
    if( zIsTiny( a - b ) ) return x;
    e3 = eval( x, util );
    if( e1 > e2 ){
      e1 = e3;
      a = x;
    } else{
      e2 = e3;
      b = x;
    }
  }
  ZITERWARN( iter );
  return x;
}

/* zOptLineBrent
 * - Brent's method.
 */
double zOptLineBrent(double (*eval)(double,void*), double a, double b, void *util, int iter)
{
  register int i;
  double u, v, w, x, xm, eu, ev, ew, ex, dm, da, db;
  double d = 0, s = 0; /* second last step */
  double p, q, r, tol1, tol2;

  x = w = v = 0.5 * ( a + b );
  ex = ew = ev = eval( x, util );
  ZITERINIT( iter );
  for( i=0; i<=iter; i++ ){
    tol2 = 2 * ( tol1 = Z_OPT_EPS * fabs(x) + zTOL );
    dm = x - ( xm = 0.5 * ( a + b ) );
    if( zIsTol( dm, tol2-0.5*(b-a) ) ) return x; /* succeed. */
    da = a - x;
    db = b - x;
    if( zIsTol( s, tol1 ) ) /* golden section method */
      d = Z_GSEC * ( s = ( dm >= 0 ) ? da : db );
    else{ /* try quadratic approximation */
      r = (x-w) * (ex-ev);
      q = (x-v) * (ex-ew);
      p = (x-v) * q - (x-w) * r;
      if( ( q = 2 * (q-r) ) > 0 ) p = -p;
      q = fabs(q);
      if( !zIsTol( p, fabs(0.5*q*s) ) || p <= q * da || p >= q * db )
        /* golden section method */
        d = Z_GSEC * ( s = ( dm >= 0 ) ? da : db );
      else{ /* parabolic interpolation */
        s = d;
        u = x + ( d = p / q );
        if( u - a < tol2 || b - u < tol2 )
          d = ( dm <= 0 ? tol1 : -tol1 );
      }
    }
    /* update */
    u = x + d;
    if( ( eu = eval( u, util ) ) <= ex ){
      ( d >= 0 ) ? ( a = x ) : ( b = x );
      v = w; w = x; x = u;
      ev=ew; ew=ex; ex=eu;
    } else{
      ( d < 0 ) ? ( a = u ) : ( b = u );
      if( eu <= ew || w == x ){
        v = w; w = u;
        ev=ew; ew=eu;
      } else if( eu <= ev || v == x || v == w ){
        v = u;
        ev=eu;
      }
    }
  }
  ZITERWARN( iter );
  return x;
}
