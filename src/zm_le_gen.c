/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_gen - linear equation: generalized inverse matrix.
 */

#include <zm/zm_le.h>
#include <zm/zm_eig.h> /* for MP-inverse and SVD */

/* zLEAllocWork
 * - allocate working memory for generalized linear equation solvers.
 */
bool zLEAllocWork(zMat *m, zVec *v, zVec *s, zIndex *idx, int size)
{
  *m = zMatAllocSqr( size );
  *v = zVecAlloc( size );
  *s = zVecAlloc( size );
  *idx = zIndexCreate( size );
  return *m && *v && *s && *idx;
}

/* zLEFreeWork
 * - free working memory for generalized linear equation solvers.
 */
void zLEFreeWork(zMat m, zVec v, zVec s, zIndex idx)
{
  zMatFree( m );
  zVecFree( v );
  zVecFree( s );
  zIndexFree( idx );
}

/* zLESolveNormMinDST
 * - weighted-norm-minimizing redundant linear equation solver
 *   without checking size consistency.
 */
zVec zLESolveNormMinDST(zMat a, zVec b, zVec w, zVec ans, zMat m, zVec v, zIndex idx, zVec s)
{
  w ? zMatQuadNC( a, w, m ) : zMulMatMatTNC( a, a, m );
  if( !zLESolveGaussDST( m, b, v, idx, s ) ) return NULL;
  zMulMatTVecNC( a, v, ans );
  return w ? zVecAmpNCDRC( ans, w ) : ans;
}

/* zLESolveNormMin
 * - norm-minimizing redundant linear equation solver.
 */
zVec zLESolveNormMin(zMat a, zVec b, zVec w, zVec ans)
{
  zMat m;
  zVec bcp, v, s;
  zIndex idx;

  if( !zMatRowVecSizeIsEqual( a, b ) ||
      !zMatColVecSizeIsEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( w && !zVecSizeIsEqual( ans, w ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  bcp = zVecClone( b );
  zLEAllocWork( &m, &v, &s, &idx, zMatRowSizeNC(a) );
  if( bcp && m && v && s && idx )
    ans = zLESolveNormMinDST( a, bcp, w, ans, m, v, idx, s );
  else{
    ZALLOCERROR();
    ans = NULL;
  }
  zVecFree( bcp );
  zLEFreeWork( m, v, s, idx );
  return ans;
}

/* zLESolveErrorMinDST
 * - error-minimizing inferior linear equation solver
 *   without checking size consistency.
 */
zVec zLESolveErrorMinDST(zMat a, zVec b, zVec w, zVec ans, zMat m, zVec v, zIndex idx, zVec s)
{
  if( w ) zVecAmpNCDRC( b, w );
  zMulMatTVecNC( a, b, v );
  zMatTQuadNC( a, w, m );
  return zLESolveGaussDST( m, v, ans, idx, s );
}

/* zLESolveErrorMin
 * - error-minimizing inferior linear equation solver.
 */
zVec zLESolveErrorMin(zMat a, zVec b, zVec w, zVec ans)
{
  zMat m;
  zVec bcp, v, s;
  zIndex idx;

  if( !zMatRowVecSizeIsEqual( a, b ) ||
      !zMatColVecSizeIsEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( w && !zVecSizeIsEqual( b, w ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  bcp = zVecClone( b );
  zLEAllocWork( &m, &v, &s, &idx, zMatColSizeNC(a) );
  if( bcp && m && v && s && idx )
    ans = zLESolveErrorMinDST( a, bcp, w, ans, m, v, idx, s );
  else{
    ZALLOCERROR();
    ans = NULL;
  }
  zVecFree( bcp );
  zLEFreeWork( m, v, s, idx );
  return ans;
}

/* zLESolveRefMinDST
 * - weighted-referred-norm-minimizing redundant linear
 *   equation solver without checking size consistency.
 */
zVec zLESolveRefMinDST(zMat a, zVec b, zVec w, zVec ref, zVec ans, zMat m, zVec v1, zVec v2, zIndex idx, zVec s)
{
  w ? zMatQuadNC( a, w, m ) : zMulMatMatTNC( a, a, m );
  zLEResidual( a, b, ref, v1 );
  if( !zLESolveGaussDST( m, v1, v2, idx, s ) ) return NULL;
  zMulMatTVecNC( a, v2, ans );
  if( w ) zVecAmpNCDRC( ans, w );
  return zVecAddNCDRC( ans, ref );
}

/* zLESolveRefMin
 * - referred-norm-minimizing redundant linear equation solver.
 */
zVec zLESolveRefMin(zMat a, zVec b, zVec w, zVec ref, zVec ans)
{
  zMat m;
  zVec bcp, v1, v2, s;
  zIndex idx;

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
  bcp = zVecClone( b );
  v2 = zVecAlloc( zMatRowSizeNC(a) );
  zLEAllocWork( &m, &v1, &s, &idx, zMatRowSizeNC(a) );
  if( bcp && m && v1 && v2 && s && idx )
    ans = zLESolveRefMinDST( a, bcp, w, ref, ans, m, v1, v2, idx, s );
  else{
    ZALLOCERROR();
    ans = NULL;
  }
  zVecFree( bcp );
  zVecFree( v2 );
  zLEFreeWork( m, v1, s, idx );
  return ans;
}

/* (static)
 * _zLEAllocWorkMP
 * - allocate working memory for a lienar equation solver with Moore=Penrose's inverse matrix.
 */
static bool _zLEAllocWorkMP(zMat a, zVec b, zMat *l, zMat *u, zVec *bcp, zIndex *idx);
bool _zLEAllocWorkMP(zMat a, zVec b, zMat *l, zMat *u, zVec *bcp, zIndex *idx)
{
  *bcp = zVecClone( b );
  *l = zMatAllocSqr( zMatRowSizeNC(a) );
  *u = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a) );
  *idx = zIndexCreate( zMatRowSizeNC(a) );
  return *bcp && *l && *u && *idx ? true : false;
}

/* (static)
 * _zLEFreeWorkMP
 * - free working memory for a lienar equation solver with Moore=Penrose's inverse matrix.
 */
static void _zLEFreeWorkMP(zMat *l, zMat *u, zVec *bcp, zIndex *idx);
void _zLEFreeWorkMP(zMat *l, zMat *u, zVec *bcp, zIndex *idx)
{
  zVecFree( *bcp );
  zMatFreeAO( 2, *l, *u );
  zIndexFree( *idx );
}

/* (static)
 * _zLESolveMP1
 * - compute left-lower part of the linear equation.
 */
static void _zLESolveMP1(zMat l, zMat u, zVec b, zVec c, zVec we, zMat m, zVec v, zVec s, zIndex idx1, zIndex idx2, int rank);
void _zLESolveMP1(zMat l, zMat u, zVec b, zVec c, zVec we, zMat m, zVec v, zVec s, zIndex idx1, zIndex idx2, int rank)
{
  if( rank < zMatColSizeNC(l) ){
    zMatColReg( l, rank );
    zMatRowReg( u, rank );
    zLESolveErrorMinDST( l, b, we, c, m, v, idx2, s );
  } else
    zLESolve_L( l, b, c, idx1 );
}

/* zLESolveMP
 * - generalized linear equation solver using Moore-Penrose's
 *   inverse (MP-inverse, pseudoinverse) based on LQ decomposition.
 */
zVec zLESolveMP(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  int rank;
  zMat l, q, m;
  zVec bcp, c, v, s;
  zIndex idx1, idx2;

  if( !_zLEAllocWorkMP( a, b, &l, &q, &bcp, &idx1 ) ) goto TERMINATE2;
  if( ( rank = zLQDecomp( a, l, q, idx1 ) ) == 0 )
    goto TERMINATE2; /* extremely irregular case */
  c = zVecAlloc( rank );
  zLEAllocWork( &m, &v, &s, &idx2, rank );
  if( !c || !m || !v || !s || !idx2 ) goto TERMINATE1;

  _zLESolveMP1( l, q, bcp, c, we, m, v, s, idx1, idx2, rank );
  zMatIsSqr(q) ?
    zMulMatTVec( q, c, ans ) :
    zLESolveNormMinDST( q, c, wn, ans, m, v, idx2, s );

 TERMINATE1:
  zVecFree( c );
  zLEFreeWork( m, v, s, idx2 );
 TERMINATE2:
  _zLEFreeWorkMP( &l, &q, &bcp, &idx1 );
  return ans;
}

/* zLESolveMP_LU
 * - generalized linear equation solver with MP-inverse
 *   based on LU decomposition.
 */
zVec zLESolveMP_LU(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  int rank;
  zMat l, u, m;
  zVec bcp, c, v, s;
  zIndex idx1, idx2;

  if( !_zLEAllocWorkMP( a, b, &l, &u, &bcp, &idx1 ) ) goto TERMINATE2;
  if( ( rank = zLUDecomp( a, l, u, idx1 ) ) == 0 )
    goto TERMINATE2; /* extremely irregular case */
  c = zVecAlloc( rank );
  zLEAllocWork( &m, &v, &s, &idx2, rank );
  if( !c || !m || !v || !s || !idx2 ) goto TERMINATE1;

  _zLESolveMP1( l, u, bcp, c, we, m, v, s, idx1, idx2, rank );
  zMatIsSqr(u) ?
    zLESolve_U( u, c, ans ) :
    zLESolveNormMinDST( u, c, wn, ans, m, v, idx2, s );

 TERMINATE1:
  zVecFree( c );
  zLEFreeWork( m, v, s, idx2 );
 TERMINATE2:
  _zLEFreeWorkMP( &l, &u, &bcp, &idx1 );
  return ans;
}

/* zLESolveMP_SVD
 * - generalized linear equation solver with MP-inverse
 *   based on singular value decomposition.
 */
zVec zLESolveMP_SVD(zMat a, zVec b, zVec ans)
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

/* zLESolveMPNull
 * - generalized linear equation solver using Moore-Penrose's
 *   inverse (MP-inverse, pseudoinverse) based on LQ decomposition
 *   with the null space.
 */
zVec zLESolveMPNull(zMat a, zVec b, zVec wn, zVec we, zVec ans, zMat mn)
{
  int rank;
  zMat l, q, m;
  zVec bcp, c, v, s;
  zIndex idx1, idx2;
  register int i;

  if( !_zLEAllocWorkMP( a, b, &l, &q, &bcp, &idx1 ) ) goto TERMINATE2;
  if( ( rank = zLQDecomp( a, l, q, idx1 ) ) == 0 )
    goto TERMINATE2; /* extremely irregular case */
  c = zVecAlloc( rank );
  zLEAllocWork( &m, &v, &s, &idx2, rank );
  if( !c || !m || !v || !s || !idx2 ) goto TERMINATE1;

  _zLESolveMP1( l, q, bcp, c, we, m, v, s, idx1, idx2, rank );
  if( zMatIsSqr(q) ){
    zMulMatTVec( q, c, ans );
    zMatClear( mn );
  } else{
    zLESolveNormMinDST( q, c, wn, ans, m, v, idx2, s );
    zMulMatTMat( q, q, mn );
    for( i=0; i<zMatRowSizeNC(mn); i++ )
      zMatElem(mn,i,i) -= 1.0;
  }

 TERMINATE1:
  zVecFree( c );
  zLEFreeWork( m, v, s, idx2 );
 TERMINATE2:
  _zLEFreeWorkMP( &l, &q, &bcp, &idx1 );
  return ans;
}

/* zLESolveMPAux
 * - generalized linear equation solver using MP-inverse
 *   biasing a vector in the null space.
 */
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

/* (static)
 * _zLESolveSRSizeIsEqual
 * - check sizes of vectors and matrices for a linear equation solver with SR-inverse matrix.
 */
static bool _zLESolveSRSizeIsEqual(zMat a, zVec b, zVec wn, zVec we, zVec ans);
bool _zLESolveSRSizeIsEqual(zMat a, zVec b, zVec wn, zVec we, zVec ans)
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

/* zLESolveSRDST
 * - linear equation solver using singularity-robust inverse
 *   (SR-inverse) matrix (destructive).
 */
zVec zLESolveSRDST(zMat a, zVec b, zVec wn, zVec we, zVec ans, zMat m, zVec v, zIndex index, zVec s)
{
  register int i;

  if( we ) zVecAmpNCDRC( b, we );
  zMulMatTVecNC( a, b, v );
  zMatTQuadNC( a, we, m );
  for( i=0; i<zMatRowSizeNC(m); i++ )
    zMatElem(m,i,i) += zVecElem(wn,i);
  return zLESolveGaussDST( m, v, ans, index, s );
}

/* zLESolveSR
 * - linear equation solver using singularity-robust inverse
 *   (SR-inverse) matrix, proposed by Y. Nakamura(1991).
 */
zVec zLESolveSR(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  zMat m;
  zVec v, bcp, s;
  zIndex index;

  if( !_zLESolveSRSizeIsEqual( a, b, wn, we, ans ) ) return NULL;
  bcp = zVecClone( b );
  zLEAllocWork( &m, &v, &s, &index, zVecSizeNC(ans) );
  if( bcp && m && v && s && index )
    ans = zLESolveSRDST( a, bcp, wn, we, ans, m, v, index, s );
  else{
    ZALLOCERROR();
    ans = NULL;
  }
  zVecFree( bcp );
  zLEFreeWork( m, v, s, index );
  return ans;
}

/* zLESolveMPAuxDST
 * - generalized linear equation solver using SR-inverse
 *   biasing a vector in the null space (destructive).
 */
zVec zLESolveSRAuxDST(zMat a, zVec b, zVec wn, zVec we, zVec ans, zVec aux, zMat m, zVec v, zIndex idx, zVec s, zVec bb)
{
  zLEResidual( a, b, aux, bb );
  zLESolveSRDST( a, bb, wn, we, ans, m, v, idx, s );
  return zVecAddNCDRC( ans, aux );
}

/* zLESolveMPAux
 * - generalized linear equation solver using SR-inverse
 *   biasing a vector in the null space.
 */
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

/* zLESolveRSRDST
 * - linear equation solver using referred singularity-robust
 *   inverse matrix (destructive).
 */
zVec zLESolveRSRDST(zMat a, zVec b, zVec wn, zVec we, zVec ref, zVec ans, zMat m, zVec v, zIndex index, zVec s)
{
  register int i;

  if( we ) zVecAmpNCDRC( b, we );
  zMulMatTVecNC( a, b, v );
  for( i=0; i<zVecSizeNC(ref); i++ )
    zVecElem(v,i) += zVecElem(wn,i) * zVecElem(ref,i);
  zMatTQuadNC( a, we, m );
  for( i=0; i<zMatRowSizeNC(m); i++ )
    zMatElem(m,i,i) += zVecElem(wn,i);
  return zLESolveGaussDST( m, v, ans, index, s );
}

/* zLESolveRSR
 * - linear equation solver using referred singularity robust
 *   inverse matrix.
 */
zVec zLESolveRSR(zMat a, zVec b, zVec wn, zVec we, zVec ref, zVec ans)
{
  zMat m;
  zVec v, bcp, s;
  zIndex index;

  if( !_zLESolveSRSizeIsEqual( a, b, wn, we, ans ) ) return NULL;
  bcp = zVecClone( b );
  zLEAllocWork( &m, &v, &s, &index, zVecSizeNC(ans) );
  if( bcp && m && v && s && index )
    ans = zLESolveRSRDST( a, bcp, wn, we, ref, ans, m, v, index, s );
  else{
    ZALLOCERROR();
    ans = NULL;
  }
  zVecFree( bcp );
  zLEFreeWork( m, v, s, index );
  return ans;
}
