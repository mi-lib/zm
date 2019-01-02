/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nurbs - NURBS curve.
 */

#include <zm/zm_nurbs.h>

static int _zNURBS1Seg(zNURBS1 *nurbs, double t);

static double _zNURBS1Basis(zNURBS1 *nurbs, double t, int i, int r, int k);
static double _zNURBS1BasisDiff(zNURBS1 *nurbs, double t, int i, int r, int k, int diff);
static double _zNURBS1DenDiff(zNURBS1 *nurbs, double t, int s, int e, int diff);

/* zNURBS1Create
 * - create a 1-dimensional NURBS curve.
 */
bool zNURBS1Create(zNURBS1 *nurbs, zSeq *seq, int dim)
{
  register int i, j;
  zSeqListCell *cp;
  bool ret = true;

  if( zListNum(seq) <= dim ){
    ZRUNERROR( ZM_ERR_NURBS_INVDIM );
    return false;
  }
  nurbs->dim = dim;
  nurbs->knot = zVecAlloc( zListNum(seq)+dim+1 );

  zArrayAlloc( &nurbs->cparray, zNURBSCPCell, zListNum(seq) );
  if( !nurbs->knot || zNURBS1CPNum(nurbs) == 0 ){
    ZALLOCERROR();
    zNURBS1Destroy( nurbs );
    return false;
  }
  /* set knots & assign control points & initialize weight uniformly */
  for( j=0; j<dim/2+1; j++ )
    zNURBS1Knot(nurbs,j) = 0;
  i = 0;
  zListForEachRew( seq, cp ){
    zNURBS1Weight(nurbs,i) = 1.0;
    if( !( zNURBS1CP(nurbs,i) = zVecClone( cp->data.v ) ) )
      ret = false;
    zNURBS1Knot(nurbs,j) = zNURBS1Knot(nurbs,j-1) + cp->data.dt;
    j++;
    i++;
  }
  for( ; j<zNURBS1KnotNum(nurbs); j++ )
    zNURBS1Knot(nurbs,j) = zNURBS1Knot(nurbs,j-1);
  if( !ret )
    zNURBS1Destroy( nurbs );
  return ret;
}

/* zNURBS1Destroy
 * - destroy a 1-dimensional NURBS curve.
 */
void zNURBS1Destroy(zNURBS1 *nurbs)
{
  register int i;

  nurbs->dim = 0;
  zVecFree( nurbs->knot );
  nurbs->knot = NULL;
  for( i=0; i<zNURBS1CPNum(nurbs); i++ )
    zVecFree( zNURBS1CP(nurbs,i) );
  zArrayFree( &nurbs->cparray );
}

/* zNURBS1KnotNormalize
 * - normalize the knot vector of a 1-dimensional NURBS curve.
 */
void zNURBS1KnotNormalize(zNURBS1 *nurbs)
{
  zVecShift( nurbs->knot, -zNURBS1Knot0(nurbs) );
  zVecDivDRC( nurbs->knot, zNURBS1KnotE(nurbs) );
}

/* (static)
 * _zNURBS1Seg
 * - find a knot segment that includes the given parameter.
 */
int _zNURBS1Seg(zNURBS1 *nurbs, double t)
{
  register int i, j, k;

  if( t < zNURBS1Knot0(nurbs) ) return -1;
  if( t >= zNURBS1KnotE(nurbs) ) return -2;
  for( i=0, j=zNURBS1KnotNum(nurbs)-1; ; ){
    while( zNURBS1Knot(nurbs,i+1) == zNURBS1Knot(nurbs,i) && i < j ) i++;
    while( zNURBS1Knot(nurbs,j-1) == zNURBS1Knot(nurbs,j) && j > i ) j--;
    if( ( k = ( i + j ) / 2 ) == i ) break;
    if( zNURBS1Knot(nurbs,k) <= t )
      i = k;
    else
      j = k;
  }
  return i;
}

/* (static)
 * _zNURBS1Basis
 * - basis function of 1-dimensional NURBS.
 */
double _zNURBS1Basis(zNURBS1 *nurbs, double t, int i, int r, int k)
{
  double t1, tr1, b=0;

  if( r == 1 )
    return i == k ? 1 : 0;
  if( i > k - r + 1 ){
    t1  = zNURBS1Knot(nurbs,i);
    tr1 = zNURBS1Knot(nurbs,i+r-1);
    if( tr1 != t1 )
      b += ( t - t1 ) / ( tr1 - t1 ) * _zNURBS1Basis(nurbs,t,i,r-1,k);
  }
  if( i <= k ){
    t1  = zNURBS1Knot(nurbs,i+1);
    tr1 = zNURBS1Knot(nurbs,i+r);
    if( tr1 != t1 )
      b += ( tr1 - t ) / ( tr1 - t1 ) * _zNURBS1Basis(nurbs,t,i+1,r-1,k);
  }
  return b;
}

/* zNURBS1Vec
 * - compute a vector on a 1-dimensional NURBS curve.
 */
zVec zNURBS1Vec(zNURBS1 *nurbs, double t, zVec v)
{
  register int s, e, i;
  double b, den;

  s = _zNURBS1Seg( nurbs, t );
  if( s == -1 )
    return zVecCopy( zNURBS1CP(nurbs,0), v );
  if( s == -2 )
    return zVecCopy( zNURBS1CP(nurbs,zNURBS1CPNum(nurbs)-1), v );
  e = zMin( s+1, zNURBS1CPNum(nurbs) );
  zVecClear( v );
  for( den=0, i=zMax(s-nurbs->dim,0); i<e; i++ ){
    b = zNURBS1Weight(nurbs,i) * _zNURBS1Basis(nurbs,t,i,nurbs->dim+1,s);
    den += b;
    zVecCatNCDRC( v, b, zNURBS1CP(nurbs,i) );
  }
  return zIsTiny(den) ?
    zVecCopy( zNURBS1CP(nurbs,0), v ) : zVecDivDRC( v, den );
}

/* (static)
 * _zNURBS1BasisDiff
 * - derivative of the basis function of 1-dimensional NURBS.
 */
double _zNURBS1BasisDiff(zNURBS1 *nurbs, double t, int i, int r, int k, int diff)
{
  double dt, b=0;

  if( diff == 0 )
    return _zNURBS1Basis( nurbs, t, i, r, k );
  if( diff > nurbs->dim + 1 || diff < 0 ){
    ZRUNERROR( ZM_ERR_NURBS_INVODR );
    return NAN;
  }
  if( i > k - r + 1 ){
    if( !zIsTiny( ( dt = zNURBS1Knot(nurbs,i+r-1) - zNURBS1Knot(nurbs,i) ) ) )
      b += _zNURBS1BasisDiff(nurbs,t,i,r-1,k,diff-1) / dt;
  }
  if( i <= k ){
    if( !zIsTiny( ( dt = zNURBS1Knot(nurbs,i+r) - zNURBS1Knot(nurbs,i+1) ) ) )
      b -= _zNURBS1BasisDiff(nurbs,t,i+1,r-1,k,diff-1) / dt;
  }
  return b * ( r - 1 );
}

/* (static)
 * _zNURBS1DenDiff
 * - derivative of the denominator of 1-dimensional NURBS.
 */
double _zNURBS1DenDiff(zNURBS1 *nurbs, double t, int s, int e, int diff)
{
  register int i;
  double den;

  for( den=0, i=zMax(s-nurbs->dim,0); i<e; i++ )
    den += zNURBS1Weight(nurbs,i) * _zNURBS1BasisDiff(nurbs,t,i,nurbs->dim+1,s,diff);
  return den;
}

/* zNURBS1VecDiff
 * - compute the derivative a 1-dimensional NURBS curve.
 */
zVec zNURBS1VecDiff(zNURBS1 *nurbs, double t, zVec v, int diff)
{
  register int s, e, i;
  double den, b;
  zVec tmp;

  if( diff == 0 )
    return zNURBS1Vec( nurbs, t, v );
  if( diff > nurbs->dim + 1 || diff < 0 ){
    ZRUNERROR( ZM_ERR_NURBS_INVODR );
    return NULL;
  }
  if( ( tmp = zVecAlloc( zVecSize(zNURBS1CP(nurbs,0)) ) ) == NULL ){
    ZALLOCERROR();
    return NULL;
  }
  zVecClear( v );
  s = _zNURBS1Seg( nurbs, t );
  if( s == -1 )
    t = zNURBS1Knot0(nurbs);
  else if( s == -2 )
    t = zNURBS1KnotE(nurbs);
  e = zMin( s+1, zNURBS1CPNum(nurbs) );
  for( den=0, i=zMax(s-nurbs->dim,0); i<e; i++ ){
    b = zNURBS1Weight(nurbs,i) * _zNURBS1BasisDiff(nurbs,t,i,nurbs->dim+1,s,diff);
    den += zNURBS1Weight(nurbs,i) * _zNURBS1Basis(nurbs,t,i,nurbs->dim+1,s);
    zVecCatNCDRC( v, b, zNURBS1CP(nurbs,i) );
  }
  for( i=1; i<diff+1; i++ ){
    if( !zNURBS1VecDiff( nurbs, t, tmp, diff-i ) ) break;
    zVecCatNCDRC( v, -zCombi(diff,i)*_zNURBS1DenDiff(nurbs,t,s,e,i), tmp );
  }
  zVecFree( tmp );
  return zIsTiny(den) ?
    zVecCopy( zNURBS1CP(nurbs,0), v ) : zVecDivDRC( v, den );
}

/* zNURBS1VecNN
 * - nearest neighbor on a 1-dimensional NURBS curve.
 */
#define ZNURBS_NN_DIV 30
double zNURBS1VecNN(zNURBS1 *nurbs, zVec v, zVec nn)
{
  double s1, s2, s1old, s2old, sj;
  double d, dmin1, dmin2;
  zVec vs;
  register int i, j;
  int iter = 0;

  if( !( vs = zVecAlloc( zVecSizeNC(v) ) ) )
    return zNURBS1Knot0(nurbs); /* dummy */
  s1 = zNURBS1Knot0(nurbs);
  s2 = zNURBS1KnotE(nurbs);
  dmin1 = dmin2 = HUGE_VAL;
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    s1old = s1;
    s2old = s2;
    for( j=0; j<=ZNURBS_NN_DIV; j++ ){
      sj = (s2old-s1old)*j/ZNURBS_NN_DIV + s1old;
      zNURBS1Vec( nurbs, sj, vs );
      d = zVecDist( v, vs );
      if( d < dmin1 ){
        dmin2 = dmin1; s2 = s1;
        dmin1 = d;     s1 = sj;
      } else
      if( d < dmin2 ){
        dmin2 = d;
        s2 = sj;
      }
    }
    if( zIsTiny( s1 - s2 ) || zIsTiny( dmin1 - dmin2 ) ) break;
  }
  zNURBS1Vec( nurbs, ( sj = 0.5*(s1+s2) ), nn );
  zVecFree( vs );
  return sj;
}

/* for debug */

/* zNURBS1KnotFWrite
 * - output knots of a 1-dimensional NURBS curve.
 */
void zNURBS1KnotFWrite(FILE *fp, zNURBS1 *nurbs)
{
  register int i;

  fprintf( fp, "[" );
  for( i=0; i<zNURBS1KnotNum(nurbs); i++ ){
    fprintf( fp, " %g", zNURBS1Knot(nurbs,i) );
  }
  fprintf( fp, " ]\n" );
}

/* zNURBS1CPFWrite
 * - outpu control points of a 1-dimensional NURBS curve.
 */
void zNURBS1CPFWrite(FILE *fp, zNURBS1 *nurbs)
{
  register int i;

  for( i=0; i<zNURBS1CPNum(nurbs); i++ ){
    fprintf( fp, "[%03d] (%g) ", i, zNURBS1Weight(nurbs,i) );
    zVecFWrite( fp, zNURBS1CP(nurbs,i) );
  }
}
