/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_akima - interpolation: Akima's interpolation (1970).
 */

#include <zm/zm_ip.h>

/* vector on Akima interpolation */
static zVec _zIPVecAkima(const zIPData *dat, double t, zVec v)
{
  int i;
  double d, dt, d2, dt2;

  i = zIPDataSeg( dat, t );
  d2 = zSqr( ( d = zIPDataDeltaTime(dat,i+1) ) );
  dt2 = zSqr( ( dt = t - zIPDataTime(dat,i) ) );
  zVecSubNC( zIPDataSecVec(dat,i+1), zIPDataSecVec(dat,i), v );
  zVecMulDRC( v, dt2*(3-2*dt/d)/d2 );
  zVecCatNCDRC( v, dt*(1-2*dt/d+dt2/d2), *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v, dt2*(dt/d-1)/d, *zArrayElem(&dat->va,i+1) );
  return zVecAddNCDRC( v, zIPDataSecVec(dat,i) );
}

/* velocity on Akima interpolation */
static zVec _zIPVelAkima(const zIPData *dat, double t, zVec v)
{
  int i;
  double d, dt;

  i = zIPDataSeg( dat, t );
  d = zIPDataDeltaTime(dat,i+1);
  dt = ( t - zIPDataTime(dat,i) ) / d;
  zVecSubNC( zIPDataSecVec(dat,i+1), zIPDataSecVec(dat,i), v );
  zVecMulDRC( v, 6*dt*(1-dt)/d );
  zVecCatNCDRC( v, 1-4*dt+3*_zSqr(dt), *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v, dt*(3*dt-2), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* acceleration on Akima interpolation */
static zVec _zIPAccAkima(const zIPData *dat, double t, zVec v)
{
  int i;
  double d, dt, d2;

  i = zIPDataSeg( dat, t );
  d2 = zSqr( ( d = zIPDataDeltaTime(dat,i+1) ) );
  dt = t - zIPDataTime(dat,i);
  zVecSubNC( zIPDataSecVec(dat,i+1), zIPDataSecVec(dat,i), v );
  zVecMulDRC( v, 6*(1-2/d*dt)/d2 );
  zVecCatNCDRC( v, 6/d2*dt-4/d, *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v, 6/d2*dt-2/d, *zArrayElem(&dat->va,i+1) );
  return v;
}

/* velocity at a section on Akima interpolation */
static zVec _zIPSecVelAkima(const zIPData *dat, int i, zVec v)
{
  if( i>=zIPDataSize(dat) || i==0 ) return zVecZero( v );
  return zVecCopyNC( *zArrayElem(&dat->va,i), v );
}

/* acceleration at a section on Akima interpolation */
static zVec _zIPSecAccAkima(const zIPData *dat, int i, zVec v)
{
  double d;
  zVec a1, a2;

  a1 = zVecAlloc( zVecSizeNC(v) );
  a2 = zVecAlloc( zVecSizeNC(v) );
  if( !a1 || !a2 ){
    v = NULL;
    goto TERMINATE;
  }
  if( i >= zIPDataSize(dat) ) return zVecZero( v );
  if( i < zIPDataSize(dat)-1 ){
    d = zIPDataDeltaTime(dat,i+1);
    zVecSubNC( zIPDataSecVec(dat,i+1), zIPDataSecVec(dat,i), a1 );
    zVecMulNCDRC( a1, 6/_zSqr(d) );
    zVecCatNCDRC( a1,-4/d, *zArrayElem(&dat->va,i) );
    zVecCatNCDRC( a1,-2/d, *zArrayElem(&dat->va,i+1) );
    if( i == 0 ){
      zVecCopyNC( a1, v );
      goto TERMINATE;
    }
    i--;
  }
  d = zIPDataDeltaTime(dat,i+1);
  zVecSubNC( zIPDataSecVec(dat,i+1), zIPDataSecVec(dat,i), a2 );
  zVecMulNCDRC( a2, -6/_zSqr(d) );
  zVecCatNCDRC( a2, 4/d, *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( a2, 2/d, *zArrayElem(&dat->va,i+1) );
  zVecAddNC( a1, a2, v );
  zVecMulDRC( v, 0.5 );

 TERMINATE:
  zVecFree( a1 );
  zVecFree( a2 );
  return v;
}

/* methods */
static zIPCom _zm_ip_com_akima = {
  _zIPVecAkima,
  _zIPVelAkima,
  _zIPAccAkima,
  _zIPSecVelAkima,
  _zIPSecAccAkima,
};

/* create an Akima/modified-Akima interpolator */
static bool _zIPCreateAkima(zIP *ip, const zSeq *seq, double (* akima_weight)(zVec,int))
{
  int i, j, n;
  zVec m;
  double dm1, dm2;

  if( !zIPDataAlloc( &ip->dat, seq ) ) return false;
  n = zIPSize(ip);
  if( !( m = zVecAlloc( n+3 ) ) ){
    ZALLOCERROR();
    return false;
  }
  for( i=0; i<zVecSizeNC(zListHead(seq)->data.v); i++ ){
    for( j=2; j<=n; j++ )
      zVecSetElemNC( m, j,
        ( zIPSecVal(ip,j-1,i)-zIPSecVal(ip,j-2,i) ) / zIPDeltaTime(ip,j-1) );
    zVecSetElemNC( m,   0, 3*zVecElem(m,2)-2*zVecElem(m,3) );
    zVecSetElemNC( m,   1, 2*zVecElem(m,2)-  zVecElem(m,3) );
    zVecSetElemNC( m, n+1, 2*zVecElem(m,n)-2*zVecElem(m,n-1) );
    zVecSetElemNC( m, n+2, 3*zVecElem(m,n)-  zVecElem(m,n-1) );
    for( j=0; j<n; j++ ){
      dm1 = akima_weight( m, j );
      dm2 = akima_weight( m, j+2 );
      zVecArrayElem(&ip->dat.va,j,i) = zIsTiny(dm1) && zIsTiny(dm2) ?
        0.5*(dm1+dm2) : ( dm2*zVecElem(m,j+1) + dm1*zVecElem(m,j+2) ) / ( dm1 + dm2 );
    }
  }
  zVecFree( m );

  ip->com = &_zm_ip_com_akima;
  return true;
}

/* weight function for Akima interpolator. */
static double _zIPAkimaWeight(zVec m, int j)
{
  return fabs( zVecElem(m,j+1) - zVecElem(m,j) );
}

/* weight function for modified Akima interpolator. */
static double _zIPModifiedAkimaWeight(zVec m, int j)
{
  return fabs( zVecElem(m,j+1) - zVecElem(m,j) ) + 0.5 * fabs( zVecElem(m,j+1) + zVecElem(m,j) );
}

/* create an Akima interpolator */
bool zIPCreateAkima(zIP *ip, const zSeq *seq)
{
  return _zIPCreateAkima( ip, seq, _zIPAkimaWeight );
}

/* create a modified Akima interpolator */
bool zIPCreateModifiedAkima(zIP *ip, const zSeq *seq)
{
  return _zIPCreateAkima( ip, seq, _zIPModifiedAkimaWeight );
}
