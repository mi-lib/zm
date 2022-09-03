/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_gen - linear equation: generalized inverse matrix.
 */

#include <zm/zm_le.h>
#include <zm/zm_eig.h> /* for MP-inverse and SVD */

/* initialize workspace for generalized linear equation solver. */
void zLEInit(zLE *le)
{
  le->m = le->l = le->r = NULL;
  le->b = le->c = le->v1 = le->v2 = le->s = NULL;
  le->idx1 = le->idx2 = NULL;
}

/* allocate workspace for generalized linear equation solvers. */
bool zLEAlloc(zLE *le, zVec b, int size)
{
  le->m = zMatAllocSqr( size );
  le->b = b ? zVecClone( b ) : NULL;
  le->v1 = zVecAlloc( size );
  le->s = zVecAlloc( size );
  le->idx1 = zIndexCreate( size );
  return le->m && le->v1 && le->s && le->idx1;
}

/* free workspace for generalized linear equation solvers. */
void zLEFree(zLE *le)
{
  zMatFree( le->m );
  zVecFree( le->b );
  zVecFree( le->v1 );
  zVecFree( le->s );
  zIndexFree( le->idx1 );
}

/* allocate workspace for generalized linear equation solvers with reference. */
static bool _zLEAllocRef(zLE *le, zVec b, int size)
{
  le->v2 = zVecAlloc( size );
  return zLEAlloc( le, b, size ) && le->v2;
}

/* free workspace for generalized linear equation solvers with reference. */
static void _zLEFreeRef(zLE *le)
{
  zLEFree( le );
  zVecFree( le->v2 );
}

/* allocate workspace for generalized linear equation solvers based on MP inverse. */
static bool _zLEAllocMP(zLE *le, zVec b, int size)
{
  le->c = zVecAlloc( size );
  return zLEAlloc( le, b, size ) && le->c;
}

/* free workspace for generalized linear equation solvers based on MP inverse. */
static void _zLEFreeMP(zLE *le)
{
  zLEFree( le );
  zVecFree( le->c );
}

/* allocate workspace for a lienar equation solver with matrix decomposition. */
static bool _zLEAllocLR(zLE *le, zMat a)
{
  le->l = zMatAllocSqr( zMatRowSizeNC(a) );
  le->r = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a) );
  le->idx2 = zIndexCreate( zMatRowSizeNC(a) );
  return le->l && le->r && le->idx2;
}

/* free workspace for a lienar equation solver with matrix decomposition. */
static void _zLEFreeLR(zLE *le)
{
  zMatFreeAO( 2, le->l, le->r );
  zIndexFree( le->idx2 );
}

/* weighted-norm-minimizing redundant linear equation solver
 * without checking size consistency. */
zVec zLESolveNormMinDST(zMat a, zVec b, zVec w, zVec ans, zLE *le)
{
  w ? zMatQuadNC( a, w, le->m ) : zMulMatMatTNC( a, a, le->m );
  if( !zLESolveGaussDST( le->m, b, le->v1, le->idx1, le->s ) ) return NULL;
  zMulMatTVecNC( a, le->v1, ans );
  return w ? zVecAmpNCDRC( ans, w ) : ans;
}

/* norm-minimizing redundant linear equation solver. */
zVec zLESolveNormMin(zMat a, zVec b, zVec w, zVec ans)
{
  zLE le;

  if( !zMatRowVecSizeIsEqual( a, b ) ||
      !zMatColVecSizeIsEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( w && !zVecSizeIsEqual( ans, w ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  ans = zLEAlloc( &le, b, zMatRowSizeNC(a) ) ?
    zLESolveNormMinDST( a, le.b, w, ans, &le ) : NULL;
  zLEFree( &le );
  return ans;
}

/* error-minimizing inferior linear equation solver
 * without checking size consistency. */
zVec zLESolveErrorMinDST(zMat a, zVec b, zVec w, zVec ans, zLE *le)
{
  if( w ) zVecAmpNCDRC( b, w );
  zMulMatTVecNC( a, b, le->v1 );
  zMatTQuadNC( a, w, le->m );
  return zLESolveGaussDST( le->m, le->v1, ans, le->idx1, le->s );
}

/* error-minimizing inferior linear equation solver. */
zVec zLESolveErrorMin(zMat a, zVec b, zVec w, zVec ans)
{
  zLE le;

  if( !zMatRowVecSizeIsEqual( a, b ) ||
      !zMatColVecSizeIsEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( w && !zVecSizeIsEqual( b, w ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  ans = zLEAlloc( &le, b, zMatColSizeNC(a) ) ?
    zLESolveErrorMinDST( a, le.b, w, ans, &le ) : NULL;
  zLEFree( &le );
  return ans;
}

/* weighted-referred-norm-minimizing redundant linear
 * equation solver without checking size consistency. */
zVec zLESolveRefMinDST(zMat a, zVec b, zVec w, zVec ref, zVec ans, zLE *le)
{
  w ? zMatQuadNC( a, w, le->m ) : zMulMatMatTNC( a, a, le->m );
  zLEResidual( a, b, ref, le->v1 );
  if( !zLESolveGaussDST( le->m, le->v1, le->v2, le->idx1, le->s ) ) return NULL;
  zMulMatTVecNC( a, le->v2, ans );
  if( w ) zVecAmpNCDRC( ans, w );
  return zVecAddNCDRC( ans, ref );
}

/* referred-norm-minimizing redundant linear equation solver. */
zVec zLESolveRefMin(zMat a, zVec b, zVec w, zVec ref, zVec ans)
{
  zLE le;

  if( !zMatRowVecSizeIsEqual( a, b ) ||
      !zMatColVecSizeIsEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( !zVecSizeIsEqual( ref, ans ) ||
      ( w && !zVecSizeIsEqual( ans, w ) ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  ans = _zLEAllocRef( &le, b, zMatRowSizeNC(a) ) ?
    zLESolveRefMinDST( a, le.b, w, ref, ans, &le ) : NULL;
  _zLEFreeRef( &le );
  return ans;
}

/* compute left-lower part of the linear equation. */
static void _zLESolveMP1(zLE *le, zVec we, int rank)
{
  if( rank < zMatColSizeNC(le->l) ){
    zMatColReg( le->l, rank );
    zMatRowReg( le->r, rank );
    zLESolveErrorMinDST( le->l, le->b, we, le->c, le );
  } else
    zLESolveL( le->l, le->b, le->c, le->idx2 );
}

/* generalized linear equation solver using Moore-Penrose's
 * inverse (MP-inverse, pseudoinverse) based on LQ decomposition. */
zVec zLESolveMP(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  int rank;
  zLE le;

  if( !_zLEAllocLR( &le, a ) ) goto TERMINATE2;
  if( ( rank = zMatDecompLQ( a, le.l, le.r, le.idx2 ) ) == 0 )
    goto TERMINATE2; /* extremely irregular case */
  if( !_zLEAllocMP( &le, b, rank ) ) goto TERMINATE1;

  _zLESolveMP1( &le, we, rank );
  zMatIsSqr(le.r) ?
    zMulMatTVec( le.r, le.c, ans ) :
    zLESolveNormMinDST( le.r, le.c, wn, ans, &le );

 TERMINATE1:
  _zLEFreeMP( &le );
 TERMINATE2:
  _zLEFreeLR( &le );
  return ans;
}

/* generalized linear equation solver with MP-inverse
 * based on LU decomposition. */
zVec zLESolveMPLU(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  int rank;
  zLE le;

  if( !_zLEAllocLR( &le, a ) ) goto TERMINATE2;
  if( ( rank = zMatDecompLU( a, le.l, le.r, le.idx2 ) ) == 0 )
    goto TERMINATE2; /* extremely irregular case */
  if( !_zLEAllocMP( &le, b, rank ) ) goto TERMINATE1;

  _zLESolveMP1( &le, we, rank );
  zMatIsSqr(le.r) ?
    zLESolveU( le.r, le.c, ans ) :
    zLESolveNormMinDST( le.r, le.c, wn, ans, &le );

 TERMINATE1:
  _zLEFreeMP( &le );
 TERMINATE2:
  _zLEFreeLR( &le );
  return ans;
}

/* generalized linear equation solver with MP-inverse
 * based on singular value decomposition. */
zVec zLESolveMPSVD(zMat a, zVec b, zVec ans)
{
  int rank;
  zMat u, v;
  zVec sv, tmp;

  u = zMatAllocSqr( zMatRowSizeNC(a) );
  v = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a) );
  sv = zVecAlloc( zMatRowSizeNC(a) );
  tmp = zVecAlloc( zMatRowSizeNC(a) );
  if( !u || !v || !sv || !tmp ) goto TERMINATE;
  if( ( rank = zSVD( a, sv, u, v ) ) < zMatRowSizeNC(a) ){
    zMatColReg( u, rank );
    zMatRowReg( v, rank );
    zVecSetSize( sv, rank );
    zVecSetSize( tmp, rank );
  }
  zMulMatTVecNC( u, b, tmp );
  zVecDemNCDRC( tmp, sv );
  zMulMatTVecNC( v, tmp, ans );

 TERMINATE:
  zMatFree( u );
  zMatFree( v );
  zVecFree( sv );
  zVecFree( tmp );
  return ans;
}

/* generalized linear equation solver using Moore-Penrose's inverse
 * (MP-inverse, pseudoinverse) based on LQ decomposition with the null space. */
zVec zLESolveMPNull(zMat a, zVec b, zVec wn, zVec we, zVec ans, zMat mn)
{
  int i, rank;
  zLE le;

  if( !_zLEAllocLR( &le, a ) ) goto TERMINATE2;
  if( ( rank = zMatDecompLQ( a, le.l, le.r, le.idx2 ) ) == 0 )
    goto TERMINATE2; /* extremely irregular case */
  if( !_zLEAllocMP( &le, b, rank ) ) goto TERMINATE1;

  _zLESolveMP1( &le, we, rank );
  if( zMatIsSqr(le.r) ){
    zMulMatTVec( le.r, le.c, ans );
    zMatZero( mn );
  } else{
    zLESolveNormMinDST( le.r, le.c, wn, ans, &le );
    zMulMatTMat( le.r, le.r, mn );
    for( i=0; i<zMatRowSizeNC(mn); i++ )
      zMatElemNC(mn,i,i) -= 1.0;
  }

 TERMINATE1:
  _zLEFreeMP( &le );
 TERMINATE2:
  _zLEFreeLR( &le );
  return ans;
}

/* generalized linear equation solver using MP-inverse
 * biasing a vector in the null space. */
zVec zLESolveMPAux(zMat a, zVec b, zVec wn, zVec we, zVec ans, zVec aux)
{
  zVec bb;

  if( !( bb = zVecAlloc( zVecSizeNC(b) ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  zLEResidual( a, b, aux, bb );
  zLESolveMP( a, bb, wn, we, ans );
  zVecFree( bb );
  return zVecAddNCDRC( ans, aux );
}

/* check sizes of vectors and matrices for a linear equation solver with SR-inverse matrix. */
static bool _zLESolveSRSizeIsEqual(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  if( !zMatRowVecSizeIsEqual( a, b ) ||
      !zMatColVecSizeIsEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return false;
  }
  if( !zVecSizeIsEqual( we, b ) ||
      !zVecSizeIsEqual( wn, ans ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return false;
  }
  return true;
}

/* linear equation solver using singularity-robust inverse
 * (SR-inverse) matrix (destructive). */
zVec zLESolveSRDST(zMat a, zVec b, zVec wn, zVec we, zVec ans, zLE *le)
{
  int i;

  if( we ) zVecAmpNCDRC( b, we );
  zMulMatTVecNC( a, b, le->v1 );
  zMatTQuadNC( a, we, le->m );
  for( i=0; i<zMatRowSizeNC(le->m); i++ )
    zMatElemNC(le->m,i,i) += zVecElemNC(wn,i);
  return zLESolveGaussDST( le->m, le->v1, ans, le->idx1, le->s );
}

/* linear equation solver using singularity-robust inverse
 * (SR-inverse) matrix, proposed by Y. Nakamura(1991). */
zVec zLESolveSR(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  zLE le;

  if( !_zLESolveSRSizeIsEqual( a, b, wn, we, ans ) ) return NULL;
  ans = zLEAlloc( &le, b, zVecSizeNC(ans) ) ?
    zLESolveSRDST( a, le.b, wn, we, ans, &le ) : NULL;
  zLEFree( &le );
  return ans;
}

/* generalized linear equation solver using SR-inverse
 * biasing a vector in the null space (destructive). */
zVec zLESolveSRAuxDST(zMat a, zVec b, zVec wn, zVec we, zVec ans, zVec aux, zLE *le, zVec bb)
{
  zLEResidual( a, b, aux, bb );
  zLESolveSRDST( a, bb, wn, we, ans, le );
  return zVecAddNCDRC( ans, aux );
}

/* generalized linear equation solver using SR-inverse
 * biasing a vector in the null space. */
zVec zLESolveSRAux(zMat a, zVec b, zVec wn, zVec we, zVec ans, zVec aux)
{
  zVec bb;

  if( !( bb = zVecAlloc( zVecSizeNC(b) ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  zLEResidual( a, b, aux, bb );
  zLESolveSR( a, bb, wn, we, ans );
  zVecFree( bb );
  return zVecAddNCDRC( ans, aux );
}

/* linear equation solver using referred singularity-robust
 * inverse matrix (destructive). */
zVec zLESolveRSRDST(zMat a, zVec b, zVec wn, zVec we, zVec ref, zVec ans, zLE *le)
{
  int i;

  if( we ) zVecAmpNCDRC( b, we );
  zMulMatTVecNC( a, b, le->v1 );
  for( i=0; i<zVecSizeNC(ref); i++ )
    zVecElemNC(le->v1,i) += zVecElemNC(wn,i) * zVecElemNC(ref,i);
  zMatTQuadNC( a, we, le->m );
  for( i=0; i<zMatRowSizeNC(le->m); i++ )
    zMatElemNC(le->m,i,i) += zVecElemNC(wn,i);
  return zLESolveGaussDST( le->m, le->v1, ans, le->idx1, le->s );
}

/* linear equation solver using referred singularity robust inverse matrix. */
zVec zLESolveRSR(zMat a, zVec b, zVec wn, zVec we, zVec ref, zVec ans)
{
  zLE le;

  if( !_zLESolveSRSizeIsEqual( a, b, wn, we, ans ) ) return NULL;
  ans = zLEAlloc( &le, b, zVecSizeNC(ans) ) ?
    zLESolveRSRDST( a, le.b, wn, we, ref, ans, &le ) : NULL;
  zLEFree( &le );
  return ans;
}
