/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_gen - linear equation: generalized inverse matrix.
 */

#include <zm/zm_le.h>
#include <zm/zm_mat_eig.h> /* for MP-inverse and SVD */

/* initialize workspace for generalized linear equation solver. */
void zLEWorkspaceInit(zLEWorkspace *workspace)
{
  /* for regular linear equation solver */
  workspace->m_regular = NULL;
  workspace->v_scale = NULL;
  workspace->index_le = NULL;
  workspace->b_copy = workspace->v_ispace = NULL;
  /* for referenced norm-minimization */
  workspace->v_ref = NULL;
  /* for MP-inverse */
  workspace->v_mp_ispace = NULL;
  /* for LU/LQ factorization */
  workspace->l_facto = workspace->r_facto = NULL;
  workspace->index_lu = NULL;
}

/* allocate workspace for generalized linear equation solvers. */
bool zLEWorkspaceAlloc(zLEWorkspace *workspace, int num_equation, int num_variable)
{
  workspace->m_regular = zMatAllocSqr( num_variable );
  workspace->b_copy = num_equation > 0 ? zVecAlloc( num_equation ) : NULL;
  workspace->v_ispace = zVecAlloc( num_variable );
  workspace->v_scale = zVecAlloc( num_variable );
  workspace->index_le = zIndexCreate( num_variable );
  return workspace->m_regular && workspace->v_ispace && workspace->v_scale && workspace->index_le && ( num_equation == 0 || workspace->b_copy );
}

/* allocate workspace for generalized linear equation solvers, and copy the right-hand-side vector. */
static bool _zLEWorkspaceAllocAndCopyB(zLEWorkspace *workspace, const zVec b, int num_variable)
{
  if( !zLEWorkspaceAlloc( workspace, zVecSize(b), num_variable ) ) return false;
  zVecCopy( b, workspace->b_copy );
  return true;
}

/* clone workspace for generalized linear equation solvers. */
bool zLEWorkspaceClone(zLEWorkspace *src, zLEWorkspace *cln)
{
  bool ret = true;

  if( src->m_regular && !( cln->m_regular = zMatClone( src->m_regular ) ) ) ret = false;
  if( src->v_scale && !( cln->v_scale = zVecClone( src->v_scale ) ) ) ret = false;
  if( src->index_le && !( cln->index_le = zIndexClone( src->index_le ) ) ) ret = false;

  if( src->b_copy && !( cln->b_copy = zVecClone( src->b_copy ) ) ) ret = false;
  if( src->v_ispace && !( cln->v_ispace = zVecClone( src->v_ispace ) ) ) ret = false;

  if( src->v_ref && !( cln->v_ref = zVecClone( src->v_ref ) ) ) ret = false;
  if( src->v_mp_ispace && !( cln->v_mp_ispace = zVecClone( src->v_mp_ispace ) ) ) ret = false;

  if( src->l_facto && !( cln->l_facto = zMatClone( src->l_facto ) ) ) ret = false;
  if( src->r_facto && !( cln->r_facto = zMatClone( src->r_facto ) ) ) ret = false;
  if( src->index_lu && !( cln->index_lu = zIndexClone( src->index_lu ) ) ) ret = false;
  if( !ret ) zLEWorkspaceFree( cln );
  return ret;
}

/* free workspace for generalized linear equation solvers. */
void zLEWorkspaceFree(zLEWorkspace *workspace)
{
  zMatFree( workspace->m_regular );
  zVecFree( workspace->v_scale );
  zIndexFree( workspace->index_le );

  zVecFree( workspace->b_copy );
  zVecFree( workspace->v_ispace );
}

/* allocate workspace for generalized linear equation solvers with reference. */
static bool _zLEWorkspaceAllocRef(zLEWorkspace *workspace, const zVec b, int num_variable)
{
  workspace->v_ref = zVecAlloc( num_variable );
  return _zLEWorkspaceAllocAndCopyB( workspace, b, num_variable ) && workspace->v_ref;
}

/* free workspace for generalized linear equation solvers with reference. */
static void _zLEWorkspaceFreeRef(zLEWorkspace *workspace)
{
  zLEWorkspaceFree( workspace );
  zVecFree( workspace->v_ref );
}

/* allocate workspace for generalized linear equation solvers based on MP inverse. */
static bool _zLEWorkspaceAllocMP(zLEWorkspace *workspace, int num_equation, int num_variable)
{
  workspace->v_mp_ispace = zVecAlloc( num_variable );
  return zLEWorkspaceAlloc( workspace, num_equation, num_variable ) && workspace->v_mp_ispace;
}

/* free workspace for generalized linear equation solvers based on MP inverse. */
static void _zLEWorkspaceFreeMP(zLEWorkspace *workspace)
{
  zLEWorkspaceFree( workspace );
  zVecFree( workspace->v_mp_ispace );
}

/* allocate workspace for a linear equation solver with matrix decomposition. */
static bool _zLEWorkspaceAllocLR(zLEWorkspace *workspace, int num_equation, int num_variable)
{
  workspace->l_facto = zMatAllocSqr( num_equation );
  workspace->r_facto = zMatAlloc( num_equation, num_variable );
  workspace->index_lu = zIndexCreate( num_equation );
  return workspace->l_facto && workspace->r_facto && workspace->index_lu;
}

/* free workspace for a linear equation solver with matrix decomposition. */
static void _zLEWorkspaceFreeLR(zLEWorkspace *workspace)
{
  zMatFreeAtOnce( 2, workspace->l_facto, workspace->r_facto );
  zIndexFree( workspace->index_lu );
}

/* allocate workspace for LQ/LU decomposition and generalized linear equation solvers based on MP inverse. */
bool zLEWorkspaceAllocMP(zLEWorkspace *workspace, int num_equation, int num_variable)
{
  zLEWorkspaceInit( workspace );
  if( !_zLEWorkspaceAllocLR( workspace, num_equation, num_variable ) ||
      !_zLEWorkspaceAllocMP( workspace, num_equation, num_variable ) ){
    ZALLOCERROR();
    _zLEWorkspaceFreeLR( workspace );
    _zLEWorkspaceFreeMP( workspace );
    return false;
  }
  return true;
}

/* allocate workspace for generalized linear equation solvers based on MP inverse, and copy the right-hand-side vector. */
static bool _zLEWorkspaceAllocMPAndCopyAB(zLEWorkspace *workspace, const zMat a, const zVec b)
{
  if( !zLEWorkspaceAllocMP( workspace, zMatRowSizeNC(a), zMatColSizeNC(a) ) ) return false;
  zVecCopy( b, workspace->b_copy );
  return true;
}

/* free workspace for LQ/LU decomposition and generalized linear equation solvers based on MP inverse. */
void zLEWorkspaceFreeMP(zLEWorkspace *workspace)
{
  _zLEWorkspaceFreeMP( workspace );
  _zLEWorkspaceFreeLR( workspace );
}

/* resize matrices and a vector for LQ/LU decomposition and generalized linear equation solvers based on MP inverse. */
void zLEWorkspaceResizeMP(zLEWorkspace *workspace, int rank)
{
  zMatColResize( workspace->l_facto, rank );
  zMatRowResize( workspace->r_facto, rank );
  zVecResize( workspace->v_mp_ispace, rank );
}

/* weighted-norm-minimizing redundant linear equation solver without checking size consistency. */
zVec zLESolveNormMinDST(const zMat a, const zVec b, const zVec w, zVec ans, zLEWorkspace *workspace)
{
  w ? zMatQuadNC( a, w, workspace->m_regular ) : zMulMatMatTNC( a, a, workspace->m_regular );
  if( !zLESolveGaussDST( workspace->m_regular, b, workspace->v_ispace, workspace->index_le, workspace->v_scale ) ) return NULL;
  zMulMatTVecNC( a, workspace->v_ispace, ans );
  return w ? zVecAmpNCDRC( ans, w ) : ans;
}

/* norm-minimizing redundant linear equation solver. */
zVec zLESolveNormMin(const zMat a, const zVec b, const zVec w, zVec ans)
{
  zLEWorkspace workspace;

  if( !zMatRowVecSizeEqual( a, b ) || !zMatColVecSizeEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( w && !zVecSizeEqual( ans, w ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  ans = _zLEWorkspaceAllocAndCopyB( &workspace, b, zMatRowSizeNC(a) ) ?
    zLESolveNormMinDST( a, workspace.b_copy, w, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}

/* error-minimizing inferior linear equation solver without checking size consistency. */
zVec zLESolveErrorMinDST(zMat a, zVec b, zVec w, zVec ans, zLEWorkspace *workspace)
{
  if( w ) zVecAmpNCDRC( b, w );
  zMulMatTVecNC( a, b, workspace->v_ispace );
  zMatTQuadNC( a, w, workspace->m_regular );
  return zLESolveGaussDST( workspace->m_regular, workspace->v_ispace, ans, workspace->index_le, workspace->v_scale );
}

/* error-minimizing inferior linear equation solver. */
zVec zLESolveErrorMin(const zMat a, const zVec b, const zVec w, zVec ans)
{
  zLEWorkspace workspace;

  if( !zMatRowVecSizeEqual( a, b ) || !zMatColVecSizeEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( w && !zVecSizeEqual( b, w ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  ans = _zLEWorkspaceAllocAndCopyB( &workspace, b, zMatColSizeNC(a) ) ?
    zLESolveErrorMinDST( a, workspace.b_copy, w, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}

/* weighted-referred-norm-minimizing redundant linear
 * equation solver without checking size consistency. */
zVec zLESolveRefNormMinDST(const zMat a, const zVec b, const zVec w, const zVec ref, zVec ans, zLEWorkspace *workspace)
{
  w ? zMatQuadNC( a, w, workspace->m_regular ) : zMulMatMatTNC( a, a, workspace->m_regular );
  zLEResidual( a, b, ref, workspace->v_ispace );
  if( !zLESolveGaussDST( workspace->m_regular, workspace->v_ispace, workspace->v_ref, workspace->index_le, workspace->v_scale ) ) return NULL;
  zMulMatTVecNC( a, workspace->v_ref, ans );
  if( w ) zVecAmpNCDRC( ans, w );
  return zVecAddNCDRC( ans, ref );
}

/* referred-norm-minimizing redundant linear equation solver. */
zVec zLESolveRefNormMin(const zMat a, const zVec b, const zVec w, const zVec ref, zVec ans)
{
  zLEWorkspace workspace;

  if( !zMatRowVecSizeEqual( a, b ) || !zMatColVecSizeEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( !zVecSizeEqual( ref, ans ) || ( w && !zVecSizeEqual( ans, w ) ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  ans = _zLEWorkspaceAllocRef( &workspace, b, zMatRowSizeNC(a) ) ?
    zLESolveRefNormMinDST( a, workspace.b_copy, w, ref, ans, &workspace ) : NULL;
  _zLEWorkspaceFreeRef( &workspace );
  return ans;
}

/* compute left-lower part of the linear equation. */
static void _zLESolveMP_L(zLEWorkspace *workspace, const zVec we, int rank)
{
  if( rank < zMatColSizeNC(workspace->l_facto) ){
    zMatColResize( workspace->l_facto, rank );
    zMatRowResize( workspace->r_facto, rank );
    zLESolveErrorMinDST( workspace->l_facto, workspace->b_copy, we, workspace->v_mp_ispace, workspace );
  } else
    zLESolveL( workspace->l_facto, workspace->b_copy, workspace->v_mp_ispace, workspace->index_lu );
}

/* generalized linear equation solver with MP-inverse based on LQ decomposition (destructive). */
zVec zLESolveMP_LQ_DST(zLEWorkspace *workspace, const zVec we, const zVec wn, int rank, zVec ans)
{
  _zLESolveMP_L( workspace, we, rank );
  zMatIsSqr( workspace->r_facto ) ?
    zMulMatTVec( workspace->r_facto, workspace->v_mp_ispace, ans ) :
    zLESolveNormMinDST( workspace->r_facto, workspace->v_mp_ispace, wn, ans, workspace );
  return ans;
}

/* generalized linear equation solver using MP-inverse based on LU decomposition (destructive). */
zVec zLESolveMP_LU_DST(zLEWorkspace *workspace, const zVec we, const zVec wn, int rank, zVec ans)
{
  _zLESolveMP_L( workspace, we, rank );
  zMatIsSqr( workspace->r_facto ) ?
    zLESolveU( workspace->r_facto, workspace->v_mp_ispace, ans ) :
    zLESolveNormMinDST( workspace->r_facto, workspace->v_mp_ispace, wn, ans, workspace );
  return ans;
}

/* generalized linear equation solver using MP-inverse with the null space (destructive). */
zVec zLESolveMPNullDST(zLEWorkspace *workspace, const zVec we, const zVec wn, int rank, zVec ans, zMat mn)
{
  int i;

  _zLESolveMP_L( workspace, we, rank );
  if( zMatIsSqr( workspace->r_facto ) ){
    zMulMatTVec( workspace->r_facto, workspace->v_mp_ispace, ans );
    zMatZero( mn );
  } else{
    zLESolveNormMinDST( workspace->r_facto, workspace->v_mp_ispace, wn, ans, workspace );
    zMulMatTMat( workspace->r_facto, workspace->r_facto, mn );
    for( i=0; i<zMatRowSizeNC(mn); i++ )
      zMatElemNC(mn,i,i) -= 1.0;
  }
  return ans;
}

/* generalized linear equation solver with MP-inverse based on LQ decomposition. */
zVec zLESolveMP_LQ(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans)
{
  int rank;
  zLEWorkspace workspace;

  if( !_zLEWorkspaceAllocMPAndCopyAB( &workspace, a, b ) ) goto TERMINATE;
  if( ( rank = zMatDecompLQ( a, workspace.l_facto, workspace.r_facto ) ) > 0 ){
    zLEWorkspaceResizeMP( &workspace, rank );
    zLESolveMP_LQ_DST( &workspace, we, wn, rank, ans );
  }
 TERMINATE:
  zLEWorkspaceFreeMP( &workspace );
  return ans;
}

/* generalized linear equation solver using MP-inverse based on LU decomposition. */
zVec zLESolveMP_LU(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans)
{
  int rank;
  zLEWorkspace workspace;

  if( !_zLEWorkspaceAllocMPAndCopyAB( &workspace, a, b ) ) goto TERMINATE;
  if( ( rank = zMatDecompLU( a, workspace.l_facto, workspace.r_facto, workspace.index_lu ) ) > 0 ){
    zLEWorkspaceResizeMP( &workspace, rank );
    zLESolveMP_LU_DST( &workspace, we, wn, rank, ans );
  }
 TERMINATE:
  zLEWorkspaceFreeMP( &workspace );
  return ans;
}

/* generalized linear equation solver using MP-inverse based on singular value decomposition. */
zVec zLESolveMP_SVD(const zMat a, const zVec b, zVec ans)
{
  int rank;
  zMat u, v;
  zVec sv, tmp;

  u = zMatAllocSqr( zMatRowSizeNC(a) );
  v = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a) );
  sv = zVecAlloc( zMatRowSizeNC(a) );
  tmp = zVecAlloc( zMatRowSizeNC(a) );
  if( !u || !v || !sv || !tmp ) goto TERMINATE;
  if( ( rank = zMatSVD( a, u, sv, v ) ) < zMatRowSizeNC(a) ){
    zMatColResize( u, rank );
    zMatRowResize( v, rank );
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

/* generalized linear equation solver using MP-inverse based on LQ decomposition with the null space. */
zVec zLESolveMPNull(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans, zMat mn)
{
  int rank;
  zLEWorkspace workspace;

  if( !_zLEWorkspaceAllocMPAndCopyAB( &workspace, a, b ) ) goto TERMINATE;
  if( ( rank = zMatDecompLQ( a, workspace.l_facto, workspace.r_facto ) ) > 0 ){
    zLEWorkspaceResizeMP( &workspace, rank );
    zLESolveMPNullDST( &workspace, we, wn, rank, ans, mn );
  }
 TERMINATE:
  zLEWorkspaceFreeMP( &workspace );
  return ans;
}

/* generalized linear equation solver using MP-inverse biasing a vector in the null space. */
zVec zLESolveMPAux(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans, const zVec aux)
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
static bool _zLESolveSRSizeEqual(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans)
{
  if( !zMatRowVecSizeEqual( a, b ) || !zMatColVecSizeEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return false;
  }
  if( !zVecSizeEqual( we, b ) || !zVecSizeEqual( wn, ans ) ){
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
  zMulMatTVecNC( a, b, workspace->v_ispace );
  zMatTQuadNC( a, we, workspace->m_regular );
  for( i=0; i<zMatRowSizeNC(workspace->m_regular); i++ )
    zMatElemNC(workspace->m_regular,i,i) += zVecElemNC(wn,i) + bias;
  return zLESolveGaussDST( workspace->m_regular, workspace->v_ispace, ans, workspace->index_le, workspace->v_scale );
}

/* linear equation solver using singularity-robust inverse matrix. */
zVec zLESolveSRBias(const zMat a, const zVec b, const zVec wn, const zVec we, double bias, zVec ans)
{
  zLEWorkspace workspace;

  if( !_zLESolveSRSizeEqual( a, b, wn, we, ans ) ) return NULL;
  ans = _zLEWorkspaceAllocAndCopyB( &workspace, b, zVecSizeNC(ans) ) ?
    zLESolveSRBiasDST( a, workspace.b_copy, wn, we, bias, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}

/* generalized linear equation solver using SR-inverse biasing a vector in the null space (destructive). */
zVec zLESolveSRAuxDST(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans, const zVec aux, zLEWorkspace *workspace, zVec bb)
{
  zLEResidual( a, b, aux, bb );
  zLESolveSRDST( a, bb, wn, we, ans, workspace );
  return zVecAddNCDRC( ans, aux );
}

/* generalized linear equation solver using SR-inverse biasing a vector in the null space. */
zVec zLESolveSRAux(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans, const zVec aux)
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
zVec zLESolveSRRefDST(zMat a, zVec b, zVec wn, zVec we, zVec ref, zVec ans, zLEWorkspace *workspace)
{
  int i;

  if( we ) zVecAmpNCDRC( b, we );
  zMulMatTVecNC( a, b, workspace->v_ispace );
  for( i=0; i<zVecSizeNC(ref); i++ )
    zVecElemNC(workspace->v_ispace,i) += zVecElemNC(wn,i) * zVecElemNC(ref,i);
  zMatTQuadNC( a, we, workspace->m_regular );
  for( i=0; i<zMatRowSizeNC(workspace->m_regular); i++ )
    zMatElemNC(workspace->m_regular,i,i) += zVecElemNC(wn,i);
  return zLESolveGaussDST( workspace->m_regular, workspace->v_ispace, ans, workspace->index_le, workspace->v_scale );
}

/* linear equation solver using referred singularity robust inverse matrix. */
zVec zLESolveSRRef(const zMat a, const zVec b, const zVec wn, const zVec we, const zVec ref, zVec ans)
{
  zLEWorkspace workspace;

  if( !_zLESolveSRSizeEqual( a, b, wn, we, ans ) ) return NULL;
  ans = _zLEWorkspaceAllocAndCopyB( &workspace, b, zVecSizeNC(ans) ) ?
    zLESolveSRRefDST( a, workspace.b_copy, wn, we, ref, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}
