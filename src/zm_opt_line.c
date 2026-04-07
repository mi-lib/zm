/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_line - optimization tools: line search.
 */

#include <zm/zm_opt.h>

static double _zOptLineGoldenSectionUpdate1(double *x1, double xmin, double xmax, double (*eval)(double,void*), void *util){
  return eval( ( *x1 = xmin + zGOLDENRATIO2 * ( xmax - xmin ) ), util );
}
static double _zOptLineGoldenSectionUpdate2(double *x2, double xmin, double xmax, double (*eval)(double,void*), void *util){
  return eval( ( *x2 = xmax - zGOLDENRATIO2 * ( xmax - xmin ) ), util );
}

/* golden section method. */
double zOptLineGoldenSection(double (*eval)(double,void*), double xmin, double xmax, void *util, int iter)
{
  int i;
  double x1, x2, e1, e2;

  e1 = _zOptLineGoldenSectionUpdate1( &x1, xmin, xmax, eval, util );
  e2 = _zOptLineGoldenSectionUpdate2( &x2, xmin, xmax, eval, util );
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    if( zEqual( xmin, xmax, zTOL ) ) goto TERMINATE;
    if( e1 > e2 ){
      xmin = x1;
      x1 = x2; e1 = e2;
      e2 = _zOptLineGoldenSectionUpdate2( &x2, xmin, xmax, eval, util );
    } else
    if( e1 < e2 ){
      xmax = x2;
      x2 = x1; e2 = e1;
      e1 = _zOptLineGoldenSectionUpdate1( &x1, xmin, xmax, eval, util );
    } else{
      xmin = x1;
      xmax = x2;
      e1 = _zOptLineGoldenSectionUpdate1( &x1, xmin, xmax, eval, util );
      e2 = _zOptLineGoldenSectionUpdate2( &x2, xmin, xmax, eval, util );
    }
  }
  ZITERWARN( iter );
 TERMINATE:
  return 0.5 * ( xmin + xmax );
}

/* trisection method. */
double zOptLineTrisection(double (*eval)(double,void*), double xmin, double xmax, void *util, int iter)
{
  int i;
  double x1, x2, e1, e2;

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    if( zEqual( xmin, xmax, zTOL ) ) goto TERMINATE;
    e1 = eval( ( x1 = ( 2 * xmin + xmax ) / 3.0 ), util );
    e2 = eval( ( x2 = ( xmin + 2 * xmax ) / 3.0 ), util );
    if( e1 > e2 ){
      xmin = x1;
    } else
    if( e1 < e2 ){
      xmax = x2;
    } else{
      xmin = x1;
      xmax = x2;
    }
  }
  ZITERWARN( iter );
 TERMINATE:
  return 0.5 * ( xmin + xmax );
}

/* Brent's method. */
double zOptLineBrent(double (*eval)(double,void*), double a, double b, void *util, int iter)
{
  int i;
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
    if( zIsTol( s, tol1 ) ){ /* golden section method */
      d = zGOLDENRATIO2 * ( s = ( dm >= 0 ) ? da : db );
    } else{ /* try quadratic approximation */
      r = (x-w) * (ex-ev);
      q = (x-v) * (ex-ew);
      p = (x-v) * q - (x-w) * r;
      if( ( q = 2 * (q-r) ) > 0 ) p = -p;
      q = fabs(q);
      if( !zIsTol( p, fabs(0.5*q*s) ) || p <= q * da || p >= q * db ){
        /* golden section method */
        d = zGOLDENRATIO2 * ( s = ( dm >= 0 ) ? da : db );
      } else{ /* parabolic interpolation */
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
