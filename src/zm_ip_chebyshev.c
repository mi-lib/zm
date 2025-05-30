/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_chebyshev - interpolation: Chebyshev's interpolation.
 */

#include <zm/zm_ip.h>

/* scaling of abscissa. */
static double _zIPScaleChebyshev(double t, double tmin, double tmax)
{
  return ( 2 * t - ( tmax + tmin ) ) / ( tmax - tmin );
}

/* vector on Chebyshev interpolation. */
static zVec _zIPVecChebyshev(const zIPData *dat, double t, zVec v)
{
  int i;
  double x, f1, f2, f;

  zVecZero( v );
  f1 = 1;
  f2 = x = _zIPScaleChebyshev( t, zIPDataTime(dat,0), zIPDataTime(dat,zIPDataSize(dat)-1) );
  zVecCopyNC( *zArrayElem(&dat->va,0), v );
  zVecCatNCDRC( v, x, *zArrayElem(&dat->va,1) );
  for( i=2; i<zIPDataSize(dat); i++ ){
    f = 2 * x * f2 - f1;
    zVecCatDRC( v, f, *zArrayElem(&dat->va,i) );
    f1 = f2;
    f2 = f;
  }
  return v;
}

/* velocity on Chebyshev interpolation. */
static zVec _zIPVelChebyshev(const zIPData *dat, double t, zVec v)
{
  int i;
  double x, f1, f2, f;
  double dx, df1, df2, df;

  zVecZero( v );
  f1 = 1;
  f2 = x = _zIPScaleChebyshev( t, zIPDataTime(dat,0), zIPDataTime(dat,zIPDataSize(dat)-1) );
  df1 = 0;
  df2 = dx = 2 / ( zIPDataTime(dat,zIPDataSize(dat)-1) - zIPDataTime(dat,0) );
  zVecMul( *zArrayElem(&dat->va,1), df2, v );
  for( i=2; i<zIPDataSize(dat); i++ ){
    f = 2 * x * f2 - f1;
    df = 2 * ( dx * f2 + x * df2 ) - df1;
    zVecCatDRC( v, df, *zArrayElem(&dat->va,i) );
    f1 = f2;
    f2 = f;
    df1 = df2;
    df2 = df;
  }
  return v;
}

/* acceleration on Chebyshev interpolation. */
static zVec _zIPAccChebyshev(const zIPData *dat, double t, zVec v)
{
  int i;
  double x, f1, f2, f;
  double dx, df1, df2, df;
  double ddf1, ddf2, ddf;

  zVecZero( v );
  f1 = 1;
  f2 = x = _zIPScaleChebyshev( t, zIPDataTime(dat,0), zIPDataTime(dat,zIPDataSize(dat)-1) );
  df1 = 0;
  df2 = dx = 2 / ( zIPDataTime(dat,zIPDataSize(dat)-1) - zIPDataTime(dat,0) );
  ddf1 = ddf2 = 0;
  for( i=2; i<zIPDataSize(dat); i++ ){
    f = 2 * x * f2 - f1;
    df = 2 * ( dx * f2 + x * df2 ) - df1;
    ddf = 2 * ( 2 * dx * df2 + x * ddf2 ) - ddf1;
    zVecCatDRC( v, ddf, *zArrayElem(&dat->va,i) );
    f1 = f2;
    f2 = f;
    df1 = df2;
    df2 = df;
    ddf1 = ddf2;
    ddf2 = ddf;
  }
  return v;
}

/* velocity at section on Chebyshev interpolation. */
static zVec _zIPSecVelChebyshev(const zIPData *dat, int i, zVec v)
{
  return _zIPVelChebyshev( dat, zIPDataTime(dat,i), v );
}

/* acceleration at section on Chebyshev interpolation. */
static zVec _zIPSecAccChebyshev(const zIPData *dat, int i, zVec v)
{
  return _zIPAccChebyshev( dat, zIPDataTime(dat,i), v );
}

/* methods */
static zIPCom _zm_ip_com_chebyshev = {
  _zIPVecChebyshev,
  _zIPVelChebyshev,
  _zIPAccChebyshev,
  _zIPSecVelChebyshev,
  _zIPSecAccChebyshev,
};

/* create Chebyshev interpolator. */
bool zIPCreateChebyshev(zIP *ip, const zSeq *seq)
{
  int i, j;
  zSeqCell *cp;
  zMat r, p, v;
  double x, f1, f2, f;
  bool ret = true;

  if( !zIPDataAlloc( &ip->dat, seq ) ) return false;
  cp = zListHead(seq);
  r = zMatAllocSqr( zIPSize(ip) );
  p = zMatAlloc( zIPSize(ip), zVecSizeNC(cp->data.v) );
  v = zMatAlloc( zIPSize(ip), zVecSizeNC(cp->data.v) );
  if( !r || !p || !v ){
    ZALLOCERROR();
    ret = false;
    goto TERMINATE;
  }
  for( i=0; i<zIPSize(ip); i++, cp=zListCellPrev(cp) ){
    x = _zIPScaleChebyshev( zIPTime(ip,i), zIPTime(ip,0), zIPTime(ip,zIPSize(ip)-1) );
    zMatElemNC(r,i,0) = f1 = 1;
    zMatElemNC(r,i,1) = f2 = x;
    for( j=2; j<zIPSize(ip); j++ ){
      zMatElemNC(r,i,j) = f = 2 * x * f2 - f1;
      f1 = f2;
      f2 = f;
    }
    zMatPutRowNC( p, i, cp->data.v );
  }
  zMulInvMatMat( r, p, v );
  for( i=0; i<zIPSize(ip); i++ )
    zMatGetRowNC( v, i, *zArrayElem(&ip->dat.va,i) );

  ip->com = &_zm_ip_com_chebyshev;
 TERMINATE:
  zMatFreeAtOnce( 3, r, p, v );
  return ret;
}
