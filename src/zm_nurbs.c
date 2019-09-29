/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nurbs - NURBS curve.
 */

#include <zm/zm_nurbs.h>

static void _zNURBSKnotInit(zNURBS *nurbs);
static int _zNURBSSeg(zNURBS *nurbs, double t);

static double _zNURBSBasis(zNURBS *nurbs, double t, int i, int r, int seg);
static double _zNURBSBasisDiff(zNURBS *nurbs, double t, int i, int r, int seg, int diff);
static double _zNURBSDenDiff(zNURBS *nurbs, double t, int s, int diff);

/* set knots & assign control points & initialize weight uniformly. */
void _zNURBSKnotInit(zNURBS *nurbs)
{
  register int j;

  for( j=0; j<=nurbs->dim; j++ )
    zNURBSSetKnot( nurbs, j, 0 );
  for( ; j<=zNURBSCPNum(nurbs); j++ )
    zNURBSSetKnot( nurbs, j, zNURBSKnot(nurbs,j-1) + 1 );
  for( ; j<zNURBSKnotNum(nurbs); j++ )
    zNURBSSetKnot( nurbs, j, zNURBSKnot(nurbs,j-1) );
}

/* create a NURBS curve. */
bool zNURBSCreate(zNURBS *nurbs, zSeq *seq, int dim)
{
  register int i;
  zSeqListCell *cp;
  bool ret = true;

  if( zListSize(seq) <= dim ){
    ZRUNERROR( ZM_ERR_NURBS_INVDIM );
    return false;
  }
  nurbs->dim = dim;
  nurbs->knot = zVecAlloc( zListSize(seq)+dim+1 );

  zArrayAlloc( &nurbs->cparray, zNURBSCPCell, zListSize(seq) );
  if( !nurbs->knot || zNURBSCPNum(nurbs) == 0 ){
    ZALLOCERROR();
    zNURBSDestroy( nurbs );
    return false;
  }
  _zNURBSKnotInit( nurbs );
  i = 0;
  zListForEachRew( seq, cp ){
    zNURBSSetWeight( nurbs, i, 1.0 );
    if( !( zNURBSCP(nurbs,i) = zVecClone( cp->data.v ) ) )
      ret = false;
    i++;
  }
  if( !ret )
    zNURBSDestroy( nurbs );
  return ret;
}

/* destroy a NURBS curve. */
void zNURBSDestroy(zNURBS *nurbs)
{
  register int i;

  nurbs->dim = 0;
  zVecFree( nurbs->knot );
  nurbs->knot = NULL;
  for( i=0; i<zNURBSCPNum(nurbs); i++ )
    zVecFree( zNURBSCP(nurbs,i) );
  zArrayFree( &nurbs->cparray );
}

/* normalize knot vector of a NURBS curve. */
void zNURBSKnotNormalize(zNURBS *nurbs)
{
  zVecShift( nurbs->knot, -zVecElemNC(nurbs->knot,0) );
  zVecDivDRC( nurbs->knot, zVecElemNC(nurbs->knot,zVecSizeNC(nurbs->knot)-1) );
}

/* find a knot segment that includes the given parameter. */
int _zNURBSSeg(zNURBS *nurbs, double t)
{
  register int i, j, k;

  for( i=nurbs->dim, j=zNURBSCPNum(nurbs); ; ){
    while( zNURBSKnot(nurbs,i+1) == zNURBSKnot(nurbs,i) ) i++;
    while( zNURBSKnot(nurbs,j-1) == zNURBSKnot(nurbs,j) ) j--;
    if( j <= i + 1 ) break;
    k = ( i + j ) / 2;
    if( zNURBSKnot(nurbs,k) > t )
      j = k;
    else
      i = k;
  }
  return i;
}

/* basis function of NURBS. */
double _zNURBSBasis(zNURBS *nurbs, double t, int i, int r, int seg)
{
  double t1, tr1, b = 0;

  if( r == 0 )
    return i == seg ? 1 : 0;
  if( i >= seg - r ){
    t1  = zNURBSKnot(nurbs,i);
    tr1 = zNURBSKnot(nurbs,i+r);
    if( tr1 != t1 )
      b += ( t - t1 ) / ( tr1 - t1 ) * _zNURBSBasis(nurbs,t,i,r-1,seg);
  }
  if( i <= seg ){
    t1  = zNURBSKnot(nurbs,i+1);
    tr1 = zNURBSKnot(nurbs,i+r+1);
    if( tr1 != t1 )
      b += ( tr1 - t ) / ( tr1 - t1 ) * _zNURBSBasis(nurbs,t,i+1,r-1,seg);
  }
  return b;
}

/* compute a vector on a NURBS curve. */
zVec zNURBSVec(zNURBS *nurbs, double t, zVec v)
{
  register int s, i;
  double b, den;

  s = _zNURBSSeg( nurbs, t );
  zVecZero( v );
  for( den=0, i=s-nurbs->dim; i<=s; i++ ){
    b = zNURBSWeight(nurbs,i) * _zNURBSBasis(nurbs,t,i,nurbs->dim,s);
    den += b;
    zVecCatNCDRC( v, b, zNURBSCP(nurbs,i) );
  }
  return zIsTiny(den) ?
    zVecCopy( zNURBSCP(nurbs,0), v ) : zVecDivDRC( v, den );
}

/* derivative of the basis function of NURBS. */
double _zNURBSBasisDiff(zNURBS *nurbs, double t, int i, int r, int seg, int diff)
{
  double dt, b = 0;

  if( diff == 0 )
    return _zNURBSBasis( nurbs, t, i, r, seg );
  if( diff > nurbs->dim + 1 || diff < 0 ){
    ZRUNERROR( ZM_ERR_NURBS_INVODR );
    return NAN;
  }
  if( i >= seg - r && ( dt = zNURBSKnot(nurbs,i+r) - zNURBSKnot(nurbs,i) ) != 0 )
    b += _zNURBSBasisDiff(nurbs,t,i,r-1,seg,diff-1) / dt;
  if( i <= seg && ( dt = zNURBSKnot(nurbs,i+r+1) - zNURBSKnot(nurbs,i+1) ) != 0 )
    b -= _zNURBSBasisDiff(nurbs,t,i+1,r-1,seg,diff-1) / dt;
  return b * r;
}

/* derivative of the denominator of NURBS. */
double _zNURBSDenDiff(zNURBS *nurbs, double t, int s, int diff)
{
  register int i;
  double den;

  for( den=0, i=s-nurbs->dim; i<=s; i++ )
    den += zNURBSWeight(nurbs,i) * _zNURBSBasisDiff(nurbs,t,i,nurbs->dim,s,diff);
  return den;
}

/* compute the derivative a NURBS curve. */
zVec zNURBSVecDiff(zNURBS *nurbs, double t, int diff, zVec v)
{
  register int s, i;
  double den, b;
  zVec tmp;

  if( diff == 0 )
    return zNURBSVec( nurbs, t, v );
  if( diff > nurbs->dim + 1 || diff < 0 ){
    ZRUNERROR( ZM_ERR_NURBS_INVODR );
    return NULL;
  }
  if( ( tmp = zVecAlloc( zVecSize(zNURBSCP(nurbs,0)) ) ) == NULL ){
    ZALLOCERROR();
    return NULL;
  }
  zVecZero( v );
  s = _zNURBSSeg( nurbs, t );
  for( den=0, i=s-nurbs->dim; i<=s; i++ ){
    b = zNURBSWeight(nurbs,i) * _zNURBSBasisDiff(nurbs,t,i,nurbs->dim,s,diff);
    den += zNURBSWeight(nurbs,i) * _zNURBSBasis(nurbs,t,i,nurbs->dim,s);
    zVecCatNCDRC( v, b, zNURBSCP(nurbs,i) );
  }
  for( i=1; i<diff+1; i++ ){
    if( !zNURBSVecDiff( nurbs, t, diff-i, tmp ) ) break;
    zVecCatNCDRC( v, -zCombi(diff,i)*_zNURBSDenDiff(nurbs,t,s,i), tmp );
  }
  zVecFree( tmp );
  return zIsTiny(den) ?
    zVecCopy( zNURBSCP(nurbs,0), v ) : zVecDivDRC( v, den );
}

/* nearest neighbor on a NURBS curve. */
#define ZNURBS_NN_DIV 30
double zNURBSVecNN(zNURBS *nurbs, zVec v, zVec nn)
{
  double s1, s2, s1old, s2old, sj;
  double d, dmin1, dmin2;
  zVec vs;
  register int i, j;
  int iter = 0;

  if( !( vs = zVecAlloc( zVecSizeNC(v) ) ) )
    return zNURBSKnotS(nurbs); /* dummy */
  s1 = zNURBSKnotS(nurbs);
  s2 = zNURBSKnotE(nurbs);
  dmin1 = dmin2 = HUGE_VAL;
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    s1old = s1;
    s2old = s2;
    for( j=0; j<=ZNURBS_NN_DIV; j++ ){
      sj = (s2old-s1old)*j/ZNURBS_NN_DIV + s1old;
      zNURBSVec( nurbs, sj, vs );
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
    if( zIsEqual( s1, s2 ) || zIsEqual( dmin1, dmin2 ) ) break;
  }
  zNURBSVec( nurbs, ( sj = 0.5*(s1+s2) ), nn );
  zVecFree( vs );
  return sj;
}

/* for debug */

/* print control points of a NURBS curve out to a file. */
void zNURBSCPFPrint(FILE *fp, zNURBS *nurbs)
{
  register int i;

  for( i=0; i<zNURBSCPNum(nurbs); i++ ){
    fprintf( fp, "[%03d] (%g) ", i, zNURBSWeight(nurbs,i) );
    zVecFPrint( fp, zNURBSCP(nurbs,i) );
  }
}
