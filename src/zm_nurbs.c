/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nurbs - NURBS curve.
 */

#include <zm/zm_nurbs.h>

/* B-spline parameter */

/* allocate B-spline parameter. */
zBSplineParam *zBSplineParamAlloc(zBSplineParam *param, int order, int nc, int slice)
{
  zBSplineParamInit( param );
  if( !( param->knot = zVecAlloc( nc + order + 1 ) ) ) return NULL;
  param->order = order;
  zBSplineParamSetSlice( param, slice != 0 ? slice : ZM_BSPLINE_DEFAULT_SLICE_NUM );
  return param;
}

/* free B-spline parameters. */
void zBSplineParamFree(zBSplineParam *param)
{
  zVecFree( param->knot );
  zBSplineParamInit( param );
}

/* copy B-spline parameters. */
bool zBSplineParamCopy(zBSplineParam *src, zBSplineParam *dest)
{
  if( src->order != dest->order ){
    ZRUNERROR( ZM_ERR_NURBS_ORDERMISMATCH, src->order, dest->order );
    return false;
  }
  if( !zVecSizeIsEqual( src->knot, dest->knot ) ){
    ZRUNERROR( ZM_ERR_NURBS_KNOTNUMMISMATCH, zVecSize(src->knot), zVecSize(dest->knot) );
    return false;
  }
  zBSplineParamSetSlice( dest, src->slice );
  zVecCopyNC( src->knot, dest->knot );
  return true;
}

/* initialize knots of a B-spline parameter. */
void zBSplineParamKnotInit(zBSplineParam *param)
{
  int j, nc;

  nc = zBSplineParamCPNum( param );
  for( j=0; j<=param->order; j++ )
    zBSplineParamSetKnot( param, j, 0 );
  for( ; j<=nc; j++ )
    zBSplineParamSetKnot( param, j, zBSplineParamKnot(param,j-1) + 1 );
  for( ; j<zBSplineParamKnotNum(param); j++ )
    zBSplineParamSetKnot( param, j, zBSplineParamKnot(param,j-1) );
}

/* normalize knot vector of a NURBS curve. */
void zBSplineParamKnotNormalize(zBSplineParam *param)
{
  zVecShift( param->knot, -zVecElemNC(param->knot,0) );
  zVecDivDRC( param->knot, zVecElemNC(param->knot,zVecSizeNC(param->knot)-1) );
}

/* find a knot segment that includes the given parameter. */
int zBSplineParamSeg(zBSplineParam *param, double t)
{
  int i, j, k, nc;

  nc = zBSplineParamCPNum( param );
  for( i=param->order, j=nc; ; ){
    while( zBSplineParamKnot(param,i+1) == zBSplineParamKnot(param,i) ) i++;
    while( zBSplineParamKnot(param,j-1) == zBSplineParamKnot(param,j) ) j--;
    if( j <= i + 1 ) break;
    k = ( i + j ) / 2;
    if( zBSplineParamKnot(param,k) > t )
      j = k;
    else
      i = k;
  }
  return i;
}

/* basis function of B-spline family. */
double zBSplineParamBasis(zBSplineParam *param, double t, int i, int r, int seg)
{
  double t1, tr1, b = 0;

  if( r == 0 )
    return i == seg ? 1 : 0;
  if( i + r >= seg ){
    t1  = zBSplineParamKnot(param,i);
    tr1 = zBSplineParamKnot(param,i+r);
    if( tr1 != t1 )
      b += ( t - t1 ) / ( tr1 - t1 ) * zBSplineParamBasis(param,t,i,r-1,seg);
  }
  if( i <= seg ){
    t1  = zBSplineParamKnot(param,i+1);
    tr1 = zBSplineParamKnot(param,i+r+1);
    if( tr1 != t1 )
      b += ( tr1 - t ) / ( tr1 - t1 ) * zBSplineParamBasis(param,t,i+1,r-1,seg);
  }
  return b;
}

/* derivative of the basis function of B-spline family. */
double zBSplineParamBasisDiff(zBSplineParam *param, double t, int i, int r, int seg, int diff)
{
  double dt, b = 0;

  if( diff == 0 )
    return zBSplineParamBasis( param, t, i, r, seg );
  if( diff > param->order + 1 || diff < 0 ){
    ZRUNERROR( ZM_ERR_NURBS_INVALID_DIFFORDER );
    return NAN;
  }
  if( i >= seg - r && ( dt = zBSplineParamKnot(param,i+r) - zBSplineParamKnot(param,i) ) != 0 )
    b += zBSplineParamBasisDiff(param,t,i,r-1,seg,diff-1) / dt;
  if( i <= seg && ( dt = zBSplineParamKnot(param,i+r+1) - zBSplineParamKnot(param,i+1) ) != 0 )
    b -= zBSplineParamBasisDiff(param,t,i+1,r-1,seg,diff-1) / dt;
  return b * r;
}

/* NURBS */

/* create a NURBS curve. */
bool zNURBSCreate(zNURBS *nurbs, zSeq *seq, int order)
{
  int i;
  zSeqCell *cp;
  bool ret = true;

  if( zListSize(seq) <= order ){
    ZRUNERROR( ZM_ERR_NURBS_INVALID_ORDER );
    return false;
  }
  if( !zBSplineParamAlloc( &nurbs->param, order, zListSize(seq), 0 ) )
    return false;
  zArrayAlloc( &nurbs->cparray, zNURBSCPCell, zListSize(seq) );
  if( zNURBSCPNum(nurbs) == 0 ){
    ZALLOCERROR();
    zNURBSDestroy( nurbs );
    return false;
  }
  zBSplineParamKnotInit( &nurbs->param );
  i = 0;
  zListForEachRew( seq, cp ){
    zNURBSSetWeight( nurbs, i, ZM_NURBS_DEFAULT_CP_WEIGHT );
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
  int i;

  zBSplineParamFree( &nurbs->param );
  for( i=0; i<zNURBSCPNum(nurbs); i++ ){
    zVecFree( zNURBSCP(nurbs,i) );
  }
  zArrayFree( &nurbs->cparray );
}

/* compute a vector on a NURBS curve. */
zVec zNURBSVec(zNURBS *nurbs, double t, zVec v)
{
  int s, i;
  double b, den;

  s = zBSplineParamSeg( &nurbs->param, t );
  zVecZero( v );
  for( den=0, i=s-nurbs->param.order; i<=s; i++ ){
    b = zNURBSWeight(nurbs,i) * zBSplineParamBasis(&nurbs->param,t,i,nurbs->param.order,s);
    den += b;
    zVecCatNCDRC( v, b, zNURBSCP(nurbs,i) );
  }
  return zIsTiny(den) ?
    zVecCopy( zNURBSCP(nurbs,0), v ) : zVecDivDRC( v, den );
}

/* derivative of the denominator of NURBS. */
static double _zNURBSDenDiff(zNURBS *nurbs, double t, int s, int diff)
{
  int i;
  double den;

  for( den=0, i=s-nurbs->param.order; i<=s; i++ )
    den += zNURBSWeight(nurbs,i) * zBSplineParamBasisDiff(&nurbs->param,t,i,nurbs->param.order,s,diff);
  return den;
}

/* compute the derivative a NURBS curve. */
zVec zNURBSVecDiff(zNURBS *nurbs, double t, int diff, zVec v)
{
  int s, i;
  double den, b;
  zVec tmp;

  if( diff == 0 )
    return zNURBSVec( nurbs, t, v );
  if( diff > nurbs->param.order + 1 || diff < 0 ){
    ZRUNERROR( ZM_ERR_NURBS_INVALID_DIFFORDER );
    return NULL;
  }
  if( ( tmp = zVecAlloc( zVecSize(zNURBSCP(nurbs,0)) ) ) == NULL ){
    ZALLOCERROR();
    return NULL;
  }
  zVecZero( v );
  s = zBSplineParamSeg( &nurbs->param, t );

  for( den=0, i=s-nurbs->param.order; i<=s; i++ ){
    b = zNURBSWeight(nurbs,i) * zBSplineParamBasisDiff(&nurbs->param,t,i,nurbs->param.order,s,diff);
    den += zNURBSWeight(nurbs,i) * zBSplineParamBasis(&nurbs->param,t,i,nurbs->param.order,s);
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
  int i, j, iter = 0;

  if( !( vs = zVecAlloc( zVecSizeNC(v) ) ) )
    return zBSplineParamKnotS(&nurbs->param); /* dummy */
  s1 = zBSplineParamKnotS(&nurbs->param);
  s2 = zBSplineParamKnotE(&nurbs->param);

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
    if( zIsEqual( s1, s2, zTOL ) || zIsEqual( dmin1, dmin2, zTOL ) ) break;
  }
  zNURBSVec( nurbs, ( sj = 0.5*(s1+s2) ), nn );
  zVecFree( vs );
  return sj;
}

/* for debug */

/* print control points of a NURBS curve out to a file. */
void zNURBSCPFPrint(FILE *fp, zNURBS *nurbs)
{
  int i;

  for( i=0; i<zNURBSCPNum(nurbs); i++ ){
    fprintf( fp, "[%03d] (%g) ", i, zNURBSWeight(nurbs,i) );
    zVecFPrint( fp, zNURBSCP(nurbs,i) );
  }
}

/* B-spline */

/* compute a vector on a B-spline curve. */
zVec zBSplineVec(zBSpline *bspline, double t, zVec v)
{
  int s, i;
  double b;

  s = zBSplineParamSeg( &bspline->param, t );
  zVecZero( v );
  for( i=s-bspline->param.order; i<=s; i++ ){
    b = zBSplineParamBasis(&bspline->param,t,i,bspline->param.order,s);
    zVecCatNCDRC( v, b, zBSplineCP(bspline,i) );
  }
  return v;
}

/* compute the derivative a B-spline curve. */
zVec zBSplineVecDiff(zBSpline *bspline, double t, int diff, zVec v)
{
  int s, i;
  double b;

  if( diff == 0 )
    return zBSplineVec( bspline, t, v );
  if( diff > bspline->param.order + 1 || diff < 0 ){
    ZRUNERROR( ZM_ERR_NURBS_INVALID_DIFFORDER );
    return NULL;
  }
  zVecZero( v );
  s = zBSplineParamSeg( &bspline->param, t );

  for( i=s-bspline->param.order; i<=s; i++ ){
    b = zBSplineParamBasisDiff( &bspline->param, t, i, bspline->param.order, s, diff );
    zVecCatNCDRC( v, b, zBSplineCP(bspline,i) );
  }
  return v;
}
