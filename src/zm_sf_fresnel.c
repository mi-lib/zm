/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_sf_fresnel - special functions : Fresnel integral.
 */

#include <zm/zm_sf.h>
#include <zm/zm_complex.h>

#define Z_FRESNEL_FP_MIN ( 1.0e-30 )

/* Fresnel integral by Gilbert formula */
static bool _zFresnelIntgGilbert(double x, double *s, double *c)
{
  int i, n;
  double sum, t, t0;

  *s = sum = 0;
  *c = t = x;
  t0 = zPI_2 * x * x;
  for( n=3, i=1; i<=Z_MAX_ITER_NUM; i++, n+=2 ){
    t *= t0/i;
    sum += ( i % 4 <= 1 ? t : -t ) / n;
    if( t < fabs(sum) * zTOL ) break;
    if( i % 2 ){
      *s = sum;
      sum = *c;
    } else{
      *c = sum;
      sum = *s;
    }
  }
  if( i > Z_MAX_ITER_NUM ){
    ZITERWARN( Z_MAX_ITER_NUM );
    return false;
  }
  return true;
}

/* Fresnel integral by modified Lentz method */
static bool _zFresnelIntgLentz(double x, double *s, double *c)
{
  int i, n;
  double a, pix2, sp, cp;
  zComplex z1, b, ca, cd, cc, d, h, del, cs;
  bool ret = true;

  pix2 = zPI_2 * x * x;
  _zComplexCreate( &z1, 1.0, 0.0 );
  _zComplexCreate( &b, 1.0, -pix2*2 );
  _zComplexCreate( &cc, 1.0/Z_FRESNEL_FP_MIN, 0.0 );
  zComplexCDiv( &z1, &b, &d );
  _zComplexCopy( &d, &h );
  for( n=1, i=2; i<=Z_MAX_ITER_NUM; i++, n+=2 ){
    a = -n * (n+1);
    b.re += 4.0;
    _zComplexMul( &d, a, &cd );
    _zComplexAddDRC( &cd, &b );
    zComplexCDiv( &z1, &cd, &d );
    zComplexInv( &cc, &cd );
    _zComplexMulDRC( &cd, a );
    _zComplexAdd( &cd, &b, &cc );
    _zComplexCMul( &cc, &d, &del );
    zComplexCMulDRC( &h, &del );
    del.re -= 1.0;
    if( zComplexIsTiny( &del ) ) break;
  }
  if( i > Z_MAX_ITER_NUM ){
    ZITERWARN( Z_MAX_ITER_NUM );
    ret = false;
  }
  zSinCos( pix2, &sp, &cp );
  _zComplexCreate( &ca, (h.re + h.im)*x, (-h.re + h.im)*x );
  _zComplexCreate( &cd, cp*ca.re - sp*ca.im, cp*ca.im + sp*ca.re );
  _zComplexCreate( &cs, 0.5*(1-cd.re +cd.im), 0.5*(1-cd.re-cd.im) );
  *c = cs.re;
  *s = cs.im;
  return ret;
}

/* Fresnel integral */
bool zFresnelIntgPI_2(double x, double *s, double *c)
{
  double ax;
  const double xth = 1.5; /* a heuristic threshould */
  bool ret = true;

  if( ( ax = fabs( x ) ) < sqrt( Z_FRESNEL_FP_MIN ) ){ /* avoid failure of convergence */
    *s = 0;
    *c = ax;
  } else if( ax <= xth ){ /* Gilbert formula for a small argument */
    ret = _zFresnelIntgGilbert( ax, s, c );
  } else{ /* modified Lentz method for a large argument */
    ret = _zFresnelIntgLentz( ax, s, c );
  }
  if( x < 0 ){ /* antisymmetry */
    *c = -*c;
    *s = -*s;
  }
  return ret;
}

/* a scaled version of Fresnel integral */
bool zFresnelIntgScale(double x, double f, double *s, double *c)
{
  bool ret;
  double r;

  if( zIsTiny( f ) ){ /* sin 0x = 0, cos 0x = 1 */
    *s = 0;
    *c = x;
    return true;
  }
  r = sqrt( zPI_2/fabs(f) );
  ret = zFresnelIntgPI_2( x/r, s, c );
  *s *= f > 0 ? r : -r;
  *c *= r;
  return ret;
}

/* normalized Fresnel integral */
bool zFresnelIntg(double x, double *s, double *c)
{
  return zFresnelIntgScale( x, 1, s, c );
}

/* generalized Fresnel integral */
bool zFresnelIntgGen(double x, double f0, double f1, double f2, double *s, double *c)
{
  double s1, c1, s2, c2;
  bool ret1 = true, ret2 = true;

  if( zIsTiny( f2 ) ){
    zSinCos( f0, &s1, &c1 );
    if( zIsTiny( f1 ) ){
      *s = x * s1;
      *c = x * c1;
    } else{
      zSinCos( f0 + f1 * x, &s2, &c2 );
      *s =-( c2 - c1 ) / f1;
      *c = ( s2 - s1 ) / f1;
    }
  } else{
    ret1 = zFresnelIntgScale( x+0.5*f1/f2, f2, &s1, &c1 );
    ret2 = zFresnelIntgScale(   0.5*f1/f2, f2, &s2, &c2 );
    s1 -= s2;
    c1 -= c2;
    zSinCos( f0 - 0.25 * f1*f1 / f2, &s2, &c2 );
    *c = c1*c2 - s1*s2;
    *s = c1*s2 + s1*c2;
  }
  return ret1 && ret2;
}
