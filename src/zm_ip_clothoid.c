/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_clothoid - interpolation: clothoid curve interpolation.
 */

#include <zm/zm_ip.h>
#include <zm/zm_sf.h>

/* argument angle of terminal displacement vector of a clothoid curve segment. */
static double _zClothoidArg(double x0, double y0, double f0, double f1, double fs)
{
  double s, c;

  zFresnelIntgGen( 1, f0, f1 - f0 - fs, fs, &s, &c );
  return atan2( s - y0, c - x0 );
}

/* create a segment of clothoid curve. */
zClothoid *zClothoidCreateSegment(zClothoid *cl, double x0, double y0, double f0, double x1, double y1, double f1)
{
  double a, a0, a1, am, s, c;
  double fs, fs0, fs1;
  int iter = 0;
  register int i;

  a = atan2( y1 - y0, x1 - x0 );
  a0 = _zClothoidArg( x0, y0, f0, f1, ( fs0 =-2*zPI ) );
  a1 = _zClothoidArg( x0, y0, f0, f1, ( fs1 = 2*zPI ) );
  if( a0 > a1 ){
    zSwap( double, a0, a1 );
    zSwap( double, fs0, fs1 );
  }
  ZITERINIT( iter );
  for( fs=0, i=0; i<iter; i++ ){
    am = _zClothoidArg( x0, y0, f0, f1, fs );
    if( zIsTiny( a - am ) ) break;
    if( a > am ){
      fs0 = fs;
      a0 = am;
    } else{
      fs1 = fs;
      a1 = am;
    }
    fs = 0.5 * ( fs0 + fs1 );
  }
  if( i >= iter ) ZITERWARN( iter );
  cl->x0 = x0;
  cl->y0 = y0;
  zFresnelIntgGen( 1, ( cl->f0 = f0 ), ( cl->fc = f1 - f0 - fs ), ( cl->fs = fs ), &s, &c );
  cl->f1 = f1;
  cl->_h = sqrt( ( zSqr( x1 - x0 ) + zSqr( y1 - y0 ) ) / ( s*s + c*c ) );
  return cl;
}

/* x-y values of a clothoid curve. */
bool zClothoidXY(zClothoid *cl, double s, double *x, double *y)
{
  double dx, dy;

  if( !zFresnelIntgGen( s, cl->f0, cl->fc, cl->fs, &dy, &dx ) ) return false;
  *x = cl->x0 + cl->_h * dx;
  *y = cl->x0 + cl->_h * dy;
  return true;
}
