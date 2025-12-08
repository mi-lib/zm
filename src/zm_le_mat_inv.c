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
double zMatDet(const zMat m)
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
zMat zMatAdj(const zMat m, zMat adj)
{
  int i, j, k, l, u, v;
  zMat smat;
  zIndex idx;

  if( !zMatIsSqr( m ) || !zMatIsSqr( adj ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( !zMatSizeEqual( m, adj ) ){
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

/* inner operation of zMulInvMatMat and zMulMatInvMat. */
static zMat _zMulInvMat(const zMat m1, const zMat m2, zMat m, zIndex idx, zVec s)
{
  int i, j, k;
  int n, p, q;
  double head;
  double x;

  n = zMatRowSizeNC( m );
  zMatMatBalancingDST( m1, m2, s );
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

/* multiplication of an inverse matrix and another matrix without checking sizes. */
zMat _zMulInvMatMatNC(const zMat m1, const zMat m2, zMat m)
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

/* multiplication of an inverse matrix and another matrix. */
zMat zMulInvMatMat(const zMat m1, const zMat m2, zMat m)
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

/* multiplication of a matrix and an inverse matrix. */
zMat zMulMatInvMat(const zMat m1, const zMat m2, zMat m)
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
zMat zMatInv(const zMat m, zMat im)
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

/* multiplication of an inverse matrix and another matrix and a vector (destructive). */
zMat zMulInvMatMatAndVecDST(const zMat m1, const zMat m2, const zVec v, zMat mat_ans, zVec vec_ans, zIndex idx, zVec s)
{
  int i, j, k;
  int n, p, q;
  double head;
  double x;

  n = zVecSizeNC( v );
  zMatMatVecBalancingDST( m1, m2, v, s );
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
    zVecElemNC(v,p) *= head;
    for( j=i+1; j<n; j++ ){
      q = zIndexElemNC( idx, j );
      if( !zIsTiny( head = zMatElemNC(m1,q,i) ) ){
        for( k=i+1; k<n; k++ )
          zMatElemNC(m1,q,k) -= zMatElemNC(m1,p,k) * head;
        for( k=0; k<zMatColSizeNC(m2); k++ )
          zMatElemNC(m2,q,k) -= zMatElemNC(m2,p,k) * head;
        zVecElemNC(v,q) -= zVecElemNC(v,p) * head;
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
        x -= zMatElemNC(m1,p,k)*zMatElemNC(mat_ans,k,j);
      zMatSetElemNC( mat_ans, i, j, x );
    }
    x = zVecElemNC( v, p );
    for( j=n-1; j>i; j-- )
      x -= zMatElemNC(m1,p,j)*zVecElemNC(vec_ans,j);
    zVecSetElemNC( vec_ans, i, x );
  }
  if( s ){
    for( i=0; i<n; i++ )
      zRawVecMulDRC( zMatRowBuf(mat_ans,i), zVecElem(s,i), zMatColSizeNC(mat_ans) );
    zVecAmpDRC( vec_ans, s );
  }
  return m1;
}

/* multiplication of an inverse matrix and another matrix and a vector without checking sizes. */
zMat zMulInvMatMatAndVecNC(const zMat m1, const zMat m2, const zVec v, zMat mat_ans, zVec vec_ans)
{
  zMat mcp1, mcp2;
  zVec vcp, s;
  zIndex idx;

  mcp1 = zMatClone( m1 );
  mcp2 = zMatClone( m2 );
  vcp  = zVecClone( v );
  idx = zIndexCreate( zMatRowSizeNC(m1) );
  s = zVecAlloc( zMatRowSizeNC(m1) );
  if( !mcp1 || !mcp2 || !vcp || !idx || !s ||
      !zMulInvMatMatAndVecDST( mcp1, mcp2, vcp, mat_ans, vec_ans, idx, s ) ) mat_ans = NULL;

  zMatFree( mcp1 );
  zMatFree( mcp2 );
  zVecFree( vcp );
  zIndexFree( idx );
  zVecFree( s );
  return mat_ans;
}

/* multiplication of inverse matrix and matrix. */
zMat zMulInvMatMatAndVec(const zMat m1, const zMat m2, const zVec v, zMat mat_ans, zVec vec_ans)
{
  if( !zMatIsSqr(m1) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( zMatColSize(m1) != zMatRowSize(m2) || !zMatSizeEqual(m2,mat_ans) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  if( !zMatColVecSizeEqual( m1, v ) || !zVecSizeEqual( v, vec_ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  return zMulInvMatMatAndVecNC( m1, m2, v, mat_ans, vec_ans );
}

/* inverse matrix by Hotelling's method. */
zMat zMatInvHotelling(const zMat m, zMat im, double tol, int iter)
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
