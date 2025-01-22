/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_gen - linear equation: generalized inverse matrix.
 */

#include <zm/zm_le.h>
#include <zm/zm_eig.h> /* for MP-inverse and SVD */

/* initialize workspace for generalized linear equation solver. */
void zLEWorkspaceInit(zLEWorkspace *workspace)
{
  workspace->m = workspace->l = workspace->r = NULL;
  workspace->b = workspace->c = workspace->v1 = workspace->v2 = workspace->s = NULL;
  workspace->idx1 = workspace->idx2 = NULL;
}

/* allocate workspace for generalized linear equation solvers. */
bool zLEWorkspaceAlloc(zLEWorkspace *workspace, zVec b, int size)
{
  workspace->m = zMatAllocSqr( size );
  workspace->b = b ? zVecClone( b ) : NULL;
  workspace->v1 = zVecAlloc( size );
  workspace->s = zVecAlloc( size );
  workspace->idx1 = zIndexCreate( size );
  return workspace->m && workspace->v1 && workspace->s && workspace->idx1;
}

/* clone workspace for generalized linear equation solvers. */
bool zLEWorkspaceClone(zLEWorkspace *src, zLEWorkspace *cln)
{
  bool ret = true;

  if( src->m && !( cln->m = zMatClone( src->m ) ) ) ret = false;
  if( src->l && !( cln->l = zMatClone( src->l ) ) ) ret = false;
  if( src->r && !( cln->r = zMatClone( src->r ) ) ) ret = false;
  if( src->b && !( cln->b = zVecClone( src->b ) ) ) ret = false;
  if( src->c && !( cln->c = zVecClone( src->c ) ) ) ret = false;
  if( src->v1 && !( cln->v1 = zVecClone( src->v1 ) ) ) ret = false;
  if( src->v2 && !( cln->v2 = zVecClone( src->v2 ) ) ) ret = false;
  if( src->s && !( cln->s = zVecClone( src->s ) ) ) ret = false;
  if( src->idx1 && !( cln->idx1 = zIndexClone( src->idx1 ) ) ) ret = false;
  if( src->idx2 && !( cln->idx2 = zIndexClone( src->idx2 ) ) ) ret = false;
  if( !ret ) zLEWorkspaceFree( cln );
  return ret;
}

/* free workspace for generalized linear equation solvers. */
void zLEWorkspaceFree(zLEWorkspace *workspace)
{
  zMatFree( workspace->m );
  zVecFree( workspace->b );
  zVecFree( workspace->v1 );
  zVecFree( workspace->s );
  zIndexFree( workspace->idx1 );
}

/* allocate workspace for generalized linear equation solvers with reference. */
static bool _zLEWorkspaceAllocRef(zLEWorkspace *workspace, zVec b, int size)
{
  workspace->v2 = zVecAlloc( size );
  return zLEWorkspaceAlloc( workspace, b, size ) && workspace->v2;
}

/* free workspace for generalized linear equation solvers with reference. */
static void _zLEWorkspaceFreeRef(zLEWorkspace *workspace)
{
  zLEWorkspaceFree( workspace );
  zVecFree( workspace->v2 );
}

/* allocate workspace for generalized linear equation solvers based on MP inverse. */
static bool _zLEWorkspaceAllocMP(zLEWorkspace *workspace, zVec b, int size)
{
  workspace->c = zVecAlloc( size );
  return zLEWorkspaceAlloc( workspace, b, size ) && workspace->c;
}

/* free workspace for generalized linear equation solvers based on MP inverse. */
static void _zLEWorkspaceFreeMP(zLEWorkspace *workspace)
{
  zLEWorkspaceFree( workspace );
  zVecFree( workspace->c );
}

/* allocate workspace for a lienar equation solver with matrix decomposition. */
static bool _zLEWorkspaceAllocLR(zLEWorkspace *workspace, zMat a)
{
  workspace->l = zMatAllocSqr( zMatRowSizeNC(a) );
  workspace->r = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a) );
  workspace->idx2 = zIndexCreate( zMatRowSizeNC(a) );
  return workspace->l && workspace->r && workspace->idx2;
}

/* free workspace for a lienar equation solver with matrix decomposition. */
static void _zLEWorkspaceFreeLR(zLEWorkspace *workspace)
{
  zMatFreeAtOnce( 2, workspace->l, workspace->r );
  zIndexFree( workspace->idx2 );
}

/* weighted-norm-minimizing redundant linear equation solver without checking size consistency. */
zVec zLESolveNormMinDST(zMat a, zVec b, zVec w, zVec ans, zLEWorkspace *workspace)
{
  w ? zMatQuadNC( a, w, workspace->m ) : zMulMatMatTNC( a, a, workspace->m );
  if( !zLESolveGaussDST( workspace->m, b, workspace->v1, workspace->idx1, workspace->s ) ) return NULL;
  zMulMatTVecNC( a, workspace->v1, ans );
  return w ? zVecAmpNCDRC( ans, w ) : ans;
}

/* norm-minimizing redundant linear equation solver. */
zVec zLESolveNormMin(zMat a, zVec b, zVec w, zVec ans)
{
  zLEWorkspace workspace;

  if( !zMatRowVecSizeEqual( a, b ) ||
      !zMatColVecSizeEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( w && !zVecSizeEqual( ans, w ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  ans = zLEWorkspaceAlloc( &workspace, b, zMatRowSizeNC(a) ) ?
    zLESolveNormMinDST( a, workspace.b, w, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}

/* error-minimizing inferior linear equation solver without checking size consistency. */
zVec zLESolveErrorMinDST(zMat a, zVec b, zVec w, zVec ans, zLEWorkspace *workspace)
{
  if( w ) zVecAmpNCDRC( b, w );
  zMulMatTVecNC( a, b, workspace->v1 );
  zMatTQuadNC( a, w, workspace->m );
  return zLESolveGaussDST( workspace->m, workspace->v1, ans, workspace->idx1, workspace->s );
}

/* error-minimizing inferior linear equation solver. */
zVec zLESolveErrorMin(zMat a, zVec b, zVec w, zVec ans)
{
  zLEWorkspace workspace;

  if( !zMatRowVecSizeEqual( a, b ) ||
      !zMatColVecSizeEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( w && !zVecSizeEqual( b, w ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  ans = zLEWorkspaceAlloc( &workspace, b, zMatColSizeNC(a) ) ?
    zLESolveErrorMinDST( a, workspace.b, w, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}

/* weighted-referred-norm-minimizing redundant linear
 * equation solver without checking size consistency. */
zVec zLESolveRefMinDST(zMat a, zVec b, zVec w, zVec ref, zVec ans, zLEWorkspace *workspace)
{
  w ? zMatQuadNC( a, w, workspace->m ) : zMulMatMatTNC( a, a, workspace->m );
  zLEResidual( a, b, ref, workspace->v1 );
  if( !zLESolveGaussDST( workspace->m, workspace->v1, workspace->v2, workspace->idx1, workspace->s ) ) return NULL;
  zMulMatTVecNC( a, workspace->v2, ans );
  if( w ) zVecAmpNCDRC( ans, w );
  return zVecAddNCDRC( ans, ref );
}

/* referred-norm-minimizing redundant linear equation solver. */
zVec zLESolveRefMin(zMat a, zVec b, zVec w, zVec ref, zVec ans)
{
  zLEWorkspace workspace;

  if( !zMatRowVecSizeEqual( a, b ) ||
      !zMatColVecSizeEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( !zVecSizeEqual( ref, ans ) ||
      ( w && !zVecSizeEqual( ans, w ) ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  ans = _zLEWorkspaceAllocRef( &workspace, b, zMatRowSizeNC(a) ) ?
    zLESolveRefMinDST( a, workspace.b, w, ref, ans, &workspace ) : NULL;
  _zLEWorkspaceFreeRef( &workspace );
  return ans;
}

/* compute left-lower part of the linear equation. */
static void _zLESolveMP1(zLEWorkspace *workspace, zVec we, int rank)
{
  if( rank < (int)zMatColSizeNC(workspace->l) ){
    zMatColReg( workspace->l, rank );
    zMatRowReg( workspace->r, rank );
    zLESolveErrorMinDST( workspace->l, workspace->b, we, workspace->c, workspace );
  } else
    zLESolveL( workspace->l, workspace->b, workspace->c, workspace->idx2 );
}

/* generalized linear equation solver using Moore-Penrose's
 * inverse (MP-inverse, pseudoinverse) based on LQ decomposition. */
zVec zLESolveMP(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  int rank;
  zLEWorkspace workspace;

  if( !_zLEWorkspaceAllocLR( &workspace, a ) ) goto TERMINATE2;
  if( ( rank = zMatDecompLQ( a, workspace.l, workspace.r, workspace.idx2 ) ) == 0 )
    goto TERMINATE2; /* extremely irregular case */
  if( !_zLEWorkspaceAllocMP( &workspace, b, rank ) ) goto TERMINATE1;

  _zLESolveMP1( &workspace, we, rank );
  zMatIsSqr( workspace.r ) ?
    zMulMatTVec( workspace.r, workspace.c, ans ) :
    zLESolveNormMinDST( workspace.r, workspace.c, wn, ans, &workspace );

 TERMINATE1:
  _zLEWorkspaceFreeMP( &workspace );
 TERMINATE2:
  _zLEWorkspaceFreeLR( &workspace );
  return ans;
}

/* generalized linear equation solver with MP-inverse
 * based on LU decomposition. */
zVec zLESolveMPLU(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  int rank;
  zLEWorkspace workspace;

  if( !_zLEWorkspaceAllocLR( &workspace, a ) ) goto TERMINATE2;
  if( ( rank = zMatDecompLU( a, workspace.l, workspace.r, workspace.idx2 ) ) == 0 )
    goto TERMINATE2; /* extremely irregular case */
  if( !_zLEWorkspaceAllocMP( &workspace, b, rank ) ) goto TERMINATE1;

  _zLESolveMP1( &workspace, we, rank );
  zMatIsSqr( workspace.r ) ?
    zLESolveU( workspace.r, workspace.c, ans ) :
    zLESolveNormMinDST( workspace.r, workspace.c, wn, ans, &workspace );

 TERMINATE1:
  _zLEWorkspaceFreeMP( &workspace );
 TERMINATE2:
  _zLEWorkspaceFreeLR( &workspace );
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
  if( ( rank = zSVD( a, sv, u, v ) ) < (int)zMatRowSizeNC(a) ){
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
  zLEWorkspace workspace;

  if( !_zLEWorkspaceAllocLR( &workspace, a ) ) goto TERMINATE2;
  if( ( rank = zMatDecompLQ( a, workspace.l, workspace.r, workspace.idx2 ) ) <= 0 )
    goto TERMINATE2; /* extremely irregular case */
  if( !_zLEWorkspaceAllocMP( &workspace, b, rank ) ) goto TERMINATE1;

  _zLESolveMP1( &workspace, we, rank );
  if( zMatIsSqr( workspace.r ) ){
    zMulMatTVec( workspace.r, workspace.c, ans );
    zMatZero( mn );
  } else{
    zLESolveNormMinDST( workspace.r, workspace.c, wn, ans, &workspace );
    zMulMatTMat( workspace.r, workspace.r, mn );
    for( i=0; i<zMatRowSizeNC(mn); i++ )
      zMatElemNC(mn,i,i) -= 1.0;
  }

 TERMINATE1:
  _zLEWorkspaceFreeMP( &workspace );
 TERMINATE2:
  _zLEWorkspaceFreeLR( &workspace );
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
static bool _zLESolveSRSizeEqual(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  if( !zMatRowVecSizeEqual( a, b ) ||
      !zMatColVecSizeEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return false;
  }
  if( !zVecSizeEqual( we, b ) ||
      !zVecSizeEqual( wn, ans ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return false;
  }
  return true;
}

/* linear equation solver using singularity-robust inverse matrix (destructive). */
zVec zLESolveSRBiasDST(zMat a, zVec b, zVec wn, zVec we, double bias, zVec ans, zLEWorkspace *workspace)
{
  int i;

  if( we ) zVecAmpNCDRC( b, we );
  zMulMatTVecNC( a, b, workspace->v1 );
  zMatTQuadNC( a, we, workspace->m );
  for( i=0; i<zMatRowSizeNC(workspace->m); i++ )
    zMatElemNC(workspace->m,i,i) += zVecElemNC(wn,i) + bias;
  return zLESolveGaussDST( workspace->m, workspace->v1, ans, workspace->idx1, workspace->s );
}

/* linear equation solver using singularity-robust inverse matrix. */
zVec zLESolveSRBias(zMat a, zVec b, zVec wn, zVec we, double bias, zVec ans)
{
  zLEWorkspace workspace;

  if( !_zLESolveSRSizeEqual( a, b, wn, we, ans ) ) return NULL;
  ans = zLEWorkspaceAlloc( &workspace, b, zVecSizeNC(ans) ) ?
    zLESolveSRBiasDST( a, workspace.b, wn, we, bias, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}

/* generalized linear equation solver using SR-inverse biasing a vector in the null space (destructive). */
zVec zLESolveSRAuxDST(zMat a, zVec b, zVec wn, zVec we, zVec ans, zVec aux, zLEWorkspace *workspace, zVec bb)
{
  zLEResidual( a, b, aux, bb );
  zLESolveSRDST( a, bb, wn, we, ans, workspace );
  return zVecAddNCDRC( ans, aux );
}

/* generalized linear equation solver using SR-inverse biasing a vector in the null space. */
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

/* linear equation solver using referred singularity-robust inverse matrix (destructive). */
zVec zLESolveRSRDST(zMat a, zVec b, zVec wn, zVec we, zVec ref, zVec ans, zLEWorkspace *workspace)
{
  int i;

  if( we ) zVecAmpNCDRC( b, we );
  zMulMatTVecNC( a, b, workspace->v1 );
  for( i=0; i<zVecSizeNC(ref); i++ )
    zVecElemNC(workspace->v1,i) += zVecElemNC(wn,i) * zVecElemNC(ref,i);
  zMatTQuadNC( a, we, workspace->m );
  for( i=0; i<zMatRowSizeNC(workspace->m); i++ )
    zMatElemNC(workspace->m,i,i) += zVecElemNC(wn,i);
  return zLESolveGaussDST( workspace->m, workspace->v1, ans, workspace->idx1, workspace->s );
}

/* linear equation solver using referred singularity robust inverse matrix. */
zVec zLESolveRSR(zMat a, zVec b, zVec wn, zVec we, zVec ref, zVec ans)
{
  zLEWorkspace workspace;

  if( !_zLESolveSRSizeEqual( a, b, wn, we, ans ) ) return NULL;
  ans = zLEWorkspaceAlloc( &workspace, b, zVecSizeNC(ans) ) ?
    zLESolveRSRDST( a, workspace.b, wn, we, ref, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}
