/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_minv - linear equation: determinant and inverse matrix.
 */

#include <zm/zm_le.h>

/* determinant of matrix (destructive). */
double zMatDetDST(zMat m, zIndex idx)
{
  int i, j, k, p, q;
  double det = 1.0;

  zIndexOrder( idx, 0 );
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    p = zIndexElemNC( idx, i );
    if( p != zMatPivoting( m, idx, i, i ) ){
      det = -det;
      p = zIndexElemNC( idx, i );
    }
    if( zIsTiny( ( det *= zMatElemNC( m, p, i ) ) ) ) return 0;
    for( j=i+1; j<zMatRowSizeNC(m); j++ ){
      q = zIndexElemNC( idx, j );
      for( k=i+1; k<zMatRowSizeNC(m); k++ )
        zMatElemNC(m,q,k) -=
          zMatElemNC(m,p,k) / zMatElemNC(m,p,i) * zMatElemNC(m,q,i);
    }
  }
  return det;
}

/* determinant of matrix. */
double zMatDet(zMat m)
{
  double det = 0;
  zMat mcp;
  zIndex idx;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return 0;
  }
  mcp = zMatClone( m );
  idx = zIndexCreate( zMatRowSizeNC(m) );
  if( mcp && idx )
    det = zMatDetDST( mcp, idx );

  zMatFree( mcp );
  zIndexFree( idx );
  return det;
}

/* adjoint matrix. */
zMat zMatAdj(zMat m, zMat adj)
{
  int i, j, k, l, u, v;
  zMat smat;
  zIndex idx;

  if( !zMatIsSqr(m) || !zMatIsSqr(adj) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( !zMatSizeEqual(m,adj) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
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
          zMatSetElemNC( smat, u, v, zMatElemNC(m,k,l) );
          v++;
	}
        u++;
      }
      zMatSetElemNC( adj, j, i,
        (i+j)%2 ? -zMatDetDST( smat, idx ) : zMatDetDST( smat, idx ) );
    }
  }
 TERMINATE:
  zMatFree( smat );
  zIndexFree( idx );
  return adj;
}

/* directly make a matrix row-balanced and column-balanced. */
static void _zBalancingMatDST(zMat m1, zMat m2, zVec s)
{
  int i, j;
  double *mp1, *mp2, tmp;

  /* column balancing */
  if( s )
    for( i=0; i<zMatColSizeNC(m1); i++ ){
      zVecSetElemNC( s, i, fabs( zMatBuf(m1)[i] ) );
      for( j=1; j<zMatRowSizeNC(m1); j++ ){
        tmp = fabs( zMatElemNC(m1,j,i) );
        if( tmp > zVecElem(s,i) ) zVecSetElemNC( s, i, tmp );
      }
      if( zVecElem(s,i) == 0 ) continue;
      /* inverse column-balancing factor */
      zVecSetElemNC( s, i, 1.0 / zVecElem(s,i) );
      for( j=0; j<zMatRowSizeNC(m1); j++ )
        zMatElemNC(m1,j,i) *= zVecElemNC(s,i);
    }
  /* row balancing */
  for( mp1=zMatBuf(m1), mp2=zMatBuf(m2), i=0; i<zMatRowSizeNC(m1); mp1+=zMatColSizeNC(m1), mp2+=zMatColSizeNC(m2), i++ ){
    if( ( tmp = zDataAbsMax( mp1, zMatColSizeNC(m1), NULL ) ) == 0 )
      continue;
    zRawVecDivDRC( mp1, tmp, zMatColSizeNC(m1) );
    zRawVecDivDRC( mp2, tmp, zMatColSizeNC(m2) );
  }
}

/* inner operation of zMulInvMatMat and zMulMatInvMat. */
static zMat _zMulInvMat(zMat m1, zMat m2, zMat m, zIndex idx, zVec s)
{
  int i, j, k;
  int n, p, q;
  double head;
  double x;

  n = zMatRowSizeNC( m );
  _zBalancingMatDST( m1, m2, s );
  /* forward elimination */
  for( i=0; i<n; i++ ){
    p = zMatPivoting( m1, idx, i, i );
    if( ( head = zMatElemNC(m1,p,i) ) == 0 ){
      ZRUNERROR( ZM_ERR_MAT_SINGULAR );
      return NULL;
    }
    head = 1.0 / head;
    zMatSetElemNC( m1, p, i, 1 );
    for( j=i+1; j<n; j++ )
      zMatElemNC(m1,p,j) *= head;
    for( j=0; j<zMatColSizeNC(m2); j++ )
      zMatElemNC(m2,p,j) *= head;
    for( j=i+1; j<n; j++ ){
      q = zIndexElemNC( idx, j );
      if( !zIsTiny( head = zMatElemNC(m1,q,i) ) ){
        for( k=i+1; k<n; k++ )
          zMatElemNC(m1,q,k) -= zMatElemNC(m1,p,k) * head;
        for( k=0; k<zMatColSizeNC(m2); k++ )
          zMatElemNC(m2,q,k) -= zMatElemNC(m2,p,k) * head;
      }
      zMatSetElemNC( m1, q, i, 0 );
    }
  }
  /* backward elimination */
  for( i=n-1; i>=0; i-- ){
    p = zIndexElemNC( idx, i );
    for( j=0; j<zMatColSizeNC(m2); j++ ){
      x = zMatElemNC( m2, p, j );
      for( k=n-1; k>i; k-- )
        x -= zMatElemNC(m1,p,k)*zMatElemNC(m,k,j);
      zMatSetElemNC( m, i, j, x );
    }
  }
  if( s )
    for( i=0; i<n; i++ )
      zRawVecMulDRC( zMatRowBuf(m,i), zVecElem(s,i), zMatColSizeNC(m) );
  return m;
}

/* multiplication of inverse matrix and matrix without checking sizes. */
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

/* multiplication of inverse matrix and matrix. */
zMat zMulInvMatMat(zMat m1, zMat m2, zMat m)
/* m = m1^-1 m2 */
{
  if( !zMatIsSqr(m1) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( zMatColSize(m1) != zMatRowSize(m2) || !zMatSizeEqual(m2,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return _zMulInvMatMatNC( m1, m2, m );
}

/* multiplication of matrix and inverse matrix. */
zMat zMulMatInvMat(zMat m1, zMat m2, zMat m)
/* m = m1 m2^-1 */
{
  zMat mcp1, mcp2, mcp;
  zVec s;
  zIndex idx;

  if( !zMatIsSqr(m2) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( zMatRowSize(m1) != zMatColSize(m2) || !zMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
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

/* inverse matrix. */
zMat zMatInv(zMat m, zMat im)
{
  if( !zMatIsSqr(m) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( zMatColSize(m) != zMatRowSize(im) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  if( zMatColSize(m) == 1 ){ /* scalar case */
    zMatBuf(im)[0] = 1.0 / zMatBuf(m)[0];
  } else{
    zMatIdentNC( im );
    _zMulInvMatMatNC( m, im, im );
  }
  return im;
}

/* inverse matrix by Hotelling's method. */
zMat zMatInvHotelling(zMat m, zMat im, double tol, int iter)
{
  int i, j;
  zMat im2, tmp, mc, mn;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( !zMatSizeEqual( m, im ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
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
      zMatElemNC(tmp,j,j) -= 1.0;
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
