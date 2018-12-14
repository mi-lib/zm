/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_minv - linear equation: determinant and inverse matrix.
 */

#include <zm/zm_le.h>

static zMat _zMulInvMat(zMat m1, zMat m2, zMat m, zIndex idx, zVec s);
static void _zBalancingMatDST(zMat m1, zMat m2, zVec s);

/* zMatDetDST
 * - determinant of matrix (destructive).
 */
double zMatDetDST(zMat m, zIndex idx)
{
  register int i, j, k, p, q;
  double det = 1.0;

  zIndexOrder( idx, 0 );
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    p = zIndexElem( idx, i );
    if( p != zPivoting( m, idx, i, i ) ){
      det = -det;
      p = zIndexElem( idx, i );
    }
    if( zIsTiny( ( det *= zMatElem( m, p, i ) ) ) ) return 0;
    for( j=i+1; j<zMatRowSizeNC(m); j++ ){
      q = zIndexElem( idx, j );
      for( k=i+1; k<zMatRowSizeNC(m); k++ )
        zMatElem(m,q,k) -=
          zMatElem(m,p,k) / zMatElem(m,p,i) * zMatElem(m,q,i);
    }
  }
  return det;
}

/* zMatDet
 * - determinant of matrix.
 */
double zMatDet(zMat m)
{
  int n;
  double det = 0;
  zMat mcp;
  zIndex idx;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return 0;
  }
  n = zMatRowSize( m );
  mcp = zMatClone( m );
  idx = zIndexCreate( n );
  if( mcp && idx )
    det = zMatDetDST( mcp, idx );

  zMatFree( mcp );
  zIndexFree( idx );
  return det;
}

/* zMatAdj
 * - adjoint matrix.
 */
zMat zMatAdj(zMat m, zMat adj)
{
  register int i, j, k, l, u, v;
  zMat smat;
  zIndex idx;

  if( !zMatIsSqr(m) || !zMatIsSqr(adj) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( !zMatSizeIsEqual(m,adj) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  smat = zMatAllocSqr( zMatRowSizeNC(m)-1 );
  idx = zIndexCreate( zMatRowSizeNC(smat) );
  if( !smat || !idx ){
    ZALLOCERROR();
    adj = NULL;
    goto TERMINATE;
  }
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    for( j=0; j<zMatColSizeNC(m); j++ ){
      for( k=u=0; k<zMatRowSizeNC(m); k++ ){
        if( k == i ) continue;
        for( l=v=0; l<zMatColSizeNC(m); l++ ){
          if( l == j ) continue;
          zMatSetElem( smat, u, v, zMatElem(m,k,l) );
          v++;
	}
        u++;
      }
      zMatSetElem( adj, j, i,
        (i+j)%2 ? -zMatDetDST( smat, idx ) : zMatDetDST( smat, idx ) );
    }
  }
 TERMINATE:
  zMatFree( smat );
  zIndexFree( idx );
  return adj;
}

/* (static)
 * _zBalancingMatDST
 * - directly make a matrix row-balanced and column-balanced.
 */
void _zBalancingMatDST(zMat m1, zMat m2, zVec s)
{
  register int i, j;
  double *mp1, *mp2, tmp;

  /* column balancing */
  if( s )
    for( i=0; i<zMatColSizeNC(m1); i++ ){
      zVecSetElem( s, i, fabs( zMatElem(m1,0,i) ) );
      for( j=1; j<zMatRowSizeNC(m1); j++ ){
        tmp = fabs( zMatElem(m1,j,i) );
        if( tmp > zVecElem(s,i) ) zVecSetElem( s, i, tmp );
      }
      if( zVecElem(s,i) == 0 ) continue;
      /* inverse column-balancing factor */
      zVecSetElem( s, i, 1.0 / zVecElem(s,i) );
      for( j=0; j<zMatRowSizeNC(m1); j++ )
        zMatElem(m1,j,i) *= zVecElem(s,i);
    }
  /* row balancing */
  for( mp1=zMatBuf(m1), mp2=zMatBuf(m2), i=0; i<zMatRowSizeNC(m1); mp1+=zMatColSizeNC(m1), mp2+=zMatColSizeNC(m2), i++ ){
    if( ( tmp = zDataAbsMax( mp1, zMatColSizeNC(m1), NULL ) ) == 0 )
      continue;
    zRawVecDivDRC( mp1, tmp, zMatColSizeNC(m1) );
    zRawVecDivDRC( mp2, tmp, zMatColSizeNC(m2) );
  }
}

/* (static)
 * _zMulInvMat
 * - inner operation of zMulInvMatMat and zMulMatInvMat.
 */
zMat _zMulInvMat(zMat m1, zMat m2, zMat m, zIndex idx, zVec s)
{
  register int i, j, k;
  int n, p, q;
  double head;
  double x;

  n = zMatRowSizeNC( m );
  _zBalancingMatDST( m1, m2, s );
  /* forward elimination */
  for( i=0; i<n; i++ ){
    p = zPivoting( m1, idx, i, i );
    if( ( head = zMatElem(m1,p,i) ) == 0 ){
      ZRUNERROR( ZM_ERR_LE_SINGULAR );
      return NULL;
    }
    head = 1.0 / head;
    zMatSetElem( m1, p, i, 1 );
    for( j=i+1; j<n; j++ )
      zMatElem(m1,p,j) *= head;
    for( j=0; j<zMatColSizeNC(m2); j++ )
      zMatElem(m2,p,j) *= head;
    for( j=i+1; j<n; j++ ){
      q = zIndexElem( idx, j );
      if( !zIsTiny( head = zMatElem(m1,q,i) ) ){
        for( k=i+1; k<n; k++ )
          zMatElem(m1,q,k) -= zMatElem(m1,p,k) * head;
        for( k=0; k<zMatColSizeNC(m2); k++ )
          zMatElem(m2,q,k) -= zMatElem(m2,p,k) * head;
      }
      zMatSetElem( m1, q, i, 0 );
    }
  }
  /* backward elimination */
  for( i=n-1; i>=0; i-- ){
    p = zIndexElem( idx, i );
    for( j=0; j<zMatColSizeNC(m2); j++ ){
      x = zMatElem( m2, p, j );
      for( k=n-1; k>i; k-- )
        x -= zMatElem(m1,p,k)*zMatElem(m,k,j);
      zMatSetElem( m, i, j, x );
    }
  }
  if( s )
    for( i=0; i<n; i++ )
      zRawVecMulDRC( &zMatElem(m,i,0), zVecElem(s,i), zMatColSizeNC(m) );
  return m;
}

/* zMulInvMatMatNC
 * - multiplication of inverse matrix and matrix without checking sizes.
 */
zMat _zMulInvMatMatNC(zMat m1, zMat m2, zMat m)
/* m = m1^-1 m2 */
{
  zMat mcp1, mcp2;
  zVec s;
  zIndex idx;

  mcp1 = zMatClone( m1 );
  mcp2 = zMatClone( m2 );
  idx = zIndexCreate( zMatRowSizeNC(m1) );
  s = zVecAlloc( zMatRowSizeNC(m1) );
  if( !mcp1 || !mcp2 || !idx || !s ||
      !_zMulInvMat( mcp1, mcp2, m, idx, s ) ) m = NULL;

  zMatFree( mcp1 );
  zMatFree( mcp2 );
  zIndexFree( idx );
  zVecFree( s );
  return m;
}

/* zMulInvMatMat
 * - multiplication of inverse matrix and matrix.
 */
zMat zMulInvMatMat(zMat m1, zMat m2, zMat m)
/* m = m1^-1 m2 */
{
  if( !zMatIsSqr(m1) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( zMatColSize(m1) != zMatRowSize(m2) || !zMatSizeIsEqual(m2,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return _zMulInvMatMatNC( m1, m2, m );
}

/* zMulMatInvMat
 * - multiplication of matrix and inverse matrix.
 */
zMat zMulMatInvMat(zMat m1, zMat m2, zMat m)
/* m = m1 m2^-1 */
{
  zMat mcp1, mcp2, mcp;
  zVec s;
  zIndex idx;

  if( !zMatIsSqr(m2) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( zMatRowSize(m1)!=zMatColSize(m2) || !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  mcp1 = zMatAlloc( zMatColSizeNC(m1), zMatRowSizeNC(m1) );
  mcp2 = zMatAlloc( zMatColSizeNC(m2), zMatRowSizeNC(m2) );
  mcp  = zMatAlloc( zMatColSizeNC(m), zMatRowSizeNC(m) );
  idx = zIndexCreate( zMatRowSizeNC(m2) );
  s = zVecAlloc( zMatRowSizeNC(m2) );
  if( !mcp1 || !mcp2 || !mcp || !idx || !s ) goto TERMINATE;
  zMatTNC( m1, mcp1 );
  zMatTNC( m2, mcp2 );
  if( _zMulInvMat( mcp2, mcp1, mcp, idx, s ) )
    zMatTNC( mcp, m );
  else
    m = NULL;

 TERMINATE:
  zMatFree( mcp1 );
  zMatFree( mcp2 );
  zMatFree( mcp );
  zIndexFree( idx );
  zVecFree( s );
  return m;
}

/* zMatInv
 * - inverse matrix.
 */
zMat zMatInv(zMat m, zMat im)
{
  if( !zMatIsSqr(m) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( zMatColSize(m) != zMatRowSize(im) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  if( zMatColSize(m) == 1 ){ /* scalar case */
    im->elem[0] = 1.0 / m->elem[0];
  } else{
    zMatIdentNC( im );
    _zMulInvMatMatNC( m, im, im );
  }
  return im;
}

/* zMatInvHotelling
 * - inverse matrix by Hotelling's method.
 */
zMat zMatInvHotelling(zMat m, zMat im, double tol, int iter)
{
  register int i, j;
  zMat im2, tmp, mc, mn;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( !zMatSizeIsEqual( m, im ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  im2 = zMatAllocSqr( zMatRowSizeNC(m) );
  tmp = zMatAllocSqr( zMatRowSizeNC(m) );
  if( !im2 || !tmp ){
    ZALLOCERROR();
    im = NULL;
    goto TERMINATE;
  }

  zMatT( m, im );
  zMatDivDRC( im, zMatSqrNorm(m) );
  ZITERINIT( iter );
  for( mc=im, mn=im2, i=0; i<iter; i++ ){
    /* update */
    zMulMatMatNC( m, mc, mn );
    zMulMatMatNC( mc, mn, tmp );
    zMatMul( mc, 2, mn );
    zMatSubNCDRC( mn, tmp );
    /* evaluate */
    zMulMatMatNC( m, mn, tmp );
    for( j=0; j<zMatRowSizeNC(tmp); j++ )
      zMatElem(tmp,j,j) -= 1.0;
    if( zMatIsTol( tmp, tol ) ){
      if( mc != im ) zMatCopy( mc, im );
      goto TERMINATE;
    }
    /* swap working memory */
    zSwap( zMat, mc, mn );
  }
  ZITERWARN( iter );

 TERMINATE:
  zMatFree( im2 );
  zMatFree( tmp );
  return im;
}
