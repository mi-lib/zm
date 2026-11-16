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
  /* for Householder transformations */
  workspace->m_house = NULL;
  /* for MP-inverse */
  workspace->v_mp_ispace = NULL;
  /* for LU/LQ factorization */
  workspace->l_facto = workspace->r_facto = NULL;
  workspace->index_row = workspace->index_col = NULL;
}

/* allocate workspace for generalized linear equation solvers. */
bool zLEWorkspaceAlloc(zLEWorkspace *workspace, int num_equation, int num_variable)
{
  workspace->m_regular = zMatAllocSqr( num_variable );
  workspace->b_copy = num_equation > 0 ? zVecAlloc( num_equation ) : NULL;
  workspace->v_ispace = zVecAlloc( num_variable );
  workspace->v_scale = zVecAlloc( num_variable );
  workspace->index_le = zIndexCreate( num_variable );
  /* the matrix for Householder transformations is large enough for both the coefficient matrix and its
   * transpose, and for those of the factors of LQ/LU decomposition. */
  workspace->m_house = num_equation > 0 ?
    zMatAlloc( zMin( num_equation, num_variable ), zMax( num_equation, num_variable ) ) : NULL;
  return workspace->m_regular && workspace->v_ispace && workspace->v_scale && workspace->index_le &&
    ( num_equation == 0 || ( workspace->b_copy && workspace->m_house ) );
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
  if( src->m_house && !( cln->m_house = zMatClone( src->m_house ) ) ) ret = false;

  if( src->v_mp_ispace && !( cln->v_mp_ispace = zVecClone( src->v_mp_ispace ) ) ) ret = false;

  if( src->l_facto && !( cln->l_facto = zMatClone( src->l_facto ) ) ) ret = false;
  if( src->r_facto && !( cln->r_facto = zMatClone( src->r_facto ) ) ) ret = false;
  if( src->index_row && !( cln->index_row = zIndexClone( src->index_row ) ) ) ret = false;
  if( src->index_col && !( cln->index_col = zIndexClone( src->index_col ) ) ) ret = false;
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
  zMatFree( workspace->m_house );
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
  int minsize;

  minsize = zMin( num_equation, num_variable );
  workspace->l_facto = zMatAlloc( num_equation, minsize );
  workspace->r_facto = zMatAlloc( minsize, num_variable );
  workspace->index_row = zIndexCreate( num_equation );
  workspace->index_col = zIndexCreate( num_variable );
  return workspace->l_facto && workspace->r_facto && workspace->index_row && workspace->index_col;
}

/* free workspace for a linear equation solver with matrix decomposition. */
static void _zLEWorkspaceFreeLR(zLEWorkspace *workspace)
{
  zMatFreeAtOnce( 2, workspace->l_facto, workspace->r_facto );
  zIndexFree( workspace->index_row );
  zIndexFree( workspace->index_col );
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
  zMatResize( workspace->m_regular, rank, rank );
  zVecResize( workspace->v_ispace, rank );
  zVecResize( workspace->v_scale, rank );
  zIndexResize( workspace->index_le, rank );
}

/* Householder LQ decomposition of a row-full-rank matrix m (r x c, r <= c) in place, namely,
 *  m = ( L O ) H_{r-1} ... H_1 H_0,
 * where L is a lower triangular matrix and H_k is the Householder reflector I - beta_k v_k v_k^T.
 * The strictly lower part of L and v_k are stored in the strictly lower part and the k-th row (from
 * the k-th column) of m, respectively, and the diagonal components of L are stored in d, where
 * beta_k = -1 / ( d_k v_kk ). Since it applies the reflectors to the rows, which are contiguous in
 * memory, it is applied to the transpose of a column-full-rank matrix for QR decomposition.
 * A row is regarded as dependent on the preceding ones if the norm of its components remaining after
 * the preceding reflections is not greater than zTOL times the maximum norm of the rows of the original
 * matrix, which is the same criterion with that of zMatDecompLQ(). It returns false in that case. */
static bool _zLEHouseholderLQ(zMat m, double d[])
{
  int i, j, k, r, c;
  double *mk, *mi, s, nu, beta, tol = 0;

  r = zMatRowSizeNC(m);
  c = zMatColSizeNC(m);
  for( k=0; k<r; k++ ) /* the maximum norm of the rows */
    if( ( s = zRawVecNorm( zMatRowBufNC(m,k), c ) ) > tol ) tol = s;
  tol *= zTOL;
  for( k=0; k<r; k++ ){
    mk = zMatRowBufNC(m,k);
    if( !( ( s = zRawVecNorm( mk+k, c-k ) ) > tol ) ) return false;
    d[k] = mk[k] >= 0 ? -s : s;
    mk[k] -= d[k]; /* v_k = m_k - d_k e_k, |v_kk| = |m_kk| + s */
    beta = -1.0 / ( d[k] * mk[k] );
    for( i=k+1; i<r; i++ ){
      mi = zMatRowBufNC(m,i);
      for( nu=0, j=k; j<c; j++ ) nu += mk[j] * mi[j];
      for( nu*=beta, j=k; j<c; j++ ) mi[j] -= nu * mk[j];
    }
  }
  return true;
}

/* apply the Householder reflectors stored by _zLEHouseholderLQ() to a vector v, namely, v is replaced
 * by H_{r-1} ... H_0 v if forward is true, or by H_0 ... H_{r-1} v otherwise. */
static void _zLEHouseholderApply(const zMat m, const double d[], double v[], bool forward)
{
  int i, j, k, r, c;
  double *mk, nu;

  r = zMatRowSizeNC(m);
  c = zMatColSizeNC(m);
  for( i=0; i<r; i++ ){
    mk = zMatRowBufNC(m,( k = forward ? i : r-1-i ));
    for( nu=0, j=k; j<c; j++ ) nu += mk[j] * v[j];
    for( nu/=-d[k]*mk[k], j=k; j<c; j++ ) v[j] -= nu * mk[j];
  }
}

/* set the size of the matrix for Householder transformations in a workspace to r x c, and check if
 * the workspace is large enough. */
static bool _zLEWorkspaceSetSizeHouseholder(zLEWorkspace *workspace, int r, int c)
{
  if( !workspace->m_house || !workspace->v_ispace ||
      zMatRowCapacity(workspace->m_house) < r || zMatColCapacity(workspace->m_house) < c ||
      zVecSizeNC(workspace->v_ispace) < r ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return false;
  }
  zMatSetSizeNC( workspace->m_house, r, c );
  return true;
}

/* minimum-weighted-norm solution of a redundant linear equation a x = b, which minimizes
 * ( x - ref )^T W^-1 ( x - ref ) for W = diag(w), where ref is the null pointer for ref = 0.
 * With z = W^-1/2 ( x - ref ), it is the minimum-norm solution of a W^1/2 z = b - a ref, which is
 * found by LQ decomposition a W^1/2 = ( L O ) Q as z = Q^T ( y O )^T, L y = b - a ref. */
static zVec _zLESolveNormMinHouseholder(const zMat a, const zVec b, const zVec w, const zVec ref, zVec ans, zLEWorkspace *workspace)
{
  int i, j, r, c;
  double *d, x;

  r = zMatRowSizeNC(a);
  c = zMatColSizeNC(a);
  if( r > c ){ /* a cannot be row-full-rank */
    ZRUNERROR( ZM_ERR_MAT_SINGULAR );
    return NULL;
  }
  if( !_zLEWorkspaceSetSizeHouseholder( workspace, r, c ) ) return NULL;
  for( i=0; i<r; i++ )
    for( j=0; j<c; j++ )
      zMatSetElemNC( workspace->m_house, i, j, w ? zMatElemNC(a,i,j) * sqrt( zVecElemNC(w,j) ) : zMatElemNC(a,i,j) );
  d = zVecBufNC(workspace->v_ispace);
  if( !_zLEHouseholderLQ( workspace->m_house, d ) ){
    ZRUNERROR( ZM_ERR_MAT_SINGULAR );
    return NULL;
  }
  for( i=0; i<r; i++ ){ /* L y = b - a ref */
    x = zVecElemNC(b,i);
    if( ref ) x -= zRawVecInnerProd( zMatRowBufNC(a,i), zVecBufNC(ref), c );
    for( j=0; j<i; j++ )
      x -= zMatElemNC(workspace->m_house,i,j) * zVecElemNC(ans,j);
    zVecSetElemNC( ans, i, x / d[i] );
  }
  for( ; i<c; i++ ) zVecSetElemNC( ans, i, 0 );
  _zLEHouseholderApply( workspace->m_house, d, zVecBufNC(ans), false );
  if( w )
    for( j=0; j<c; j++ ) zVecElemNC(ans,j) *= sqrt( zVecElemNC(w,j) );
  return ref ? zVecAddNCDRC( ans, ref ) : ans;
}

/* weighted-norm-minimizing redundant linear equation solver without checking size consistency. */
zVec zLESolveNormMinDST(const zMat a, const zVec b, const zVec w, zVec ans, zLEWorkspace *workspace)
{
  return _zLESolveNormMinHouseholder( a, b, w, NULL, ans, workspace );
}

/* allocate workspace for solvers based on Householder transformations, and copy the right-hand-side
 * vector, where the matrices for the normal equations are not allocated. */
static bool _zLEWorkspaceAllocHouseholder(zLEWorkspace *workspace, const zMat a, const zVec b)
{
  int minsize, maxsize;

  zLEWorkspaceInit( workspace );
  minsize = zMin( zMatRowSizeNC(a), zMatColSizeNC(a) );
  maxsize = zMax( zMatRowSizeNC(a), zMatColSizeNC(a) );
  workspace->b_copy = zVecClone( b );
  workspace->v_ispace = zVecAlloc( minsize );
  workspace->m_house = zMatAlloc( minsize, maxsize );
  if( !workspace->b_copy || !workspace->v_ispace || !workspace->m_house ){
    ZALLOCERROR();
    return false;
  }
  return true;
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
  ans = _zLEWorkspaceAllocHouseholder( &workspace, a, b ) ?
    zLESolveNormMinDST( a, workspace.b_copy, w, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}

/* error-minimizing inferior linear equation solver without checking size consistency.
 * it minimizes ( a x - b )^T W ( a x - b ) for W = diag(w), namely, |W^1/2 a x - W^1/2 b|^2, which is
 * solved by QR decomposition W^1/2 a = Q^T ( L^T O )^T obtained by LQ decomposition of a^T W^1/2 as
 * L^T x = ( Q W^1/2 b )_upper. b is destroyed. */
zVec zLESolveErrorMinDST(zMat a, zVec b, zVec w, zVec ans, zLEWorkspace *workspace)
{
  int i, j, r, c;
  double *d, x;

  r = zMatRowSizeNC(a);
  c = zMatColSizeNC(a);
  if( c > r ){ /* a cannot be column-full-rank */
    ZRUNERROR( ZM_ERR_MAT_SINGULAR );
    return NULL;
  }
  if( !_zLEWorkspaceSetSizeHouseholder( workspace, c, r ) ) return NULL;
  for( i=0; i<r; i++ )
    for( j=0; j<c; j++ )
      zMatSetElemNC( workspace->m_house, j, i, w ? zMatElemNC(a,i,j) * sqrt( zVecElemNC(w,i) ) : zMatElemNC(a,i,j) );
  if( w )
    for( i=0; i<r; i++ ) zVecElemNC(b,i) *= sqrt( zVecElemNC(w,i) );
  d = zVecBufNC(workspace->v_ispace);
  if( !_zLEHouseholderLQ( workspace->m_house, d ) ){
    ZRUNERROR( ZM_ERR_MAT_SINGULAR );
    return NULL;
  }
  _zLEHouseholderApply( workspace->m_house, d, zVecBufNC(b), true );
  for( j=c-1; j>=0; j-- ){ /* L^T x = ( Q W^1/2 b )_upper */
    x = zVecElemNC(b,j);
    for( i=j+1; i<c; i++ )
      x -= zMatElemNC(workspace->m_house,i,j) * zVecElemNC(ans,i);
    zVecSetElemNC( ans, j, x / d[j] );
  }
  return ans;
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
  ans = _zLEWorkspaceAllocHouseholder( &workspace, a, b ) ?
    zLESolveErrorMinDST( a, workspace.b_copy, w, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}

/* weighted-referred-norm-minimizing redundant linear
 * equation solver without checking size consistency. */
zVec zLESolveRefNormMinDST(const zMat a, const zVec b, const zVec w, const zVec ref, zVec ans, zLEWorkspace *workspace)
{
  return _zLESolveNormMinHouseholder( a, b, w, ref, ans, workspace );
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
  ans = _zLEWorkspaceAllocHouseholder( &workspace, a, b ) ?
    zLESolveRefNormMinDST( a, workspace.b_copy, w, ref, ans, &workspace ) : NULL;
  zLEWorkspaceFree( &workspace );
  return ans;
}

/* The weighted MP-inverse solution, which minimizes ( a x - b )^T We ( a x - b ) and then x^T Wn^-1 x
 * for We = diag(we) and Wn = diag(wn), is x = Wn^1/2 y, where y is the MP-inverse solution of
 * We^1/2 a Wn^1/2 y = We^1/2 b. Hence, the weighted coefficient matrix We^1/2 a Wn^1/2 is decomposed
 * into a column-full-rank matrix L and a row-full-rank matrix R as L R, and y = R# L# We^1/2 b, where
 * # denotes the MP-inverse. The rank is detected for the weighted coefficient matrix. */

/* weight a coefficient matrix as We^1/2 a Wn^1/2. */
static zMat _zLEWeightMat(const zMat a, const zVec wn, const zVec we, zMat aw)
{
  int i, j;
  double s;

  zMatCopyNC( a, aw );
  if( we )
    for( i=0; i<zMatRowSizeNC(aw); i++ )
      zRawVecMulDRC( zMatRowBufNC(aw,i), sqrt( zVecElemNC(we,i) ), zMatColSizeNC(aw) );
  if( wn )
    for( j=0; j<zMatColSizeNC(aw); j++ ){
      s = sqrt( zVecElemNC(wn,j) );
      for( i=0; i<zMatRowSizeNC(aw); i++ )
        zMatElemNC(aw,i,j) *= s;
    }
  return aw;
}

/* solve L z = We^1/2 b in the sense of least squares for the left factor L. */
static zVec _zLESolveMP_L(zLEWorkspace *workspace, const zVec we)
{
  int i;

  if( we )
    for( i=0; i<zVecSizeNC(workspace->b_copy); i++ )
      zVecElemNC(workspace->b_copy,i) *= sqrt( zVecElemNC(we,i) );
  return zMatIsSqr( workspace->l_facto ) ?
    zLESolveL( workspace->l_facto, workspace->b_copy, workspace->v_mp_ispace, workspace->index_row ) :
    zLESolveErrorMinDST( workspace->l_facto, workspace->b_copy, NULL, workspace->v_mp_ispace, workspace );
}

/* scale the answer as Wn^1/2 y. */
static zVec _zLESolveMPScaleAns(zVec ans, const zVec wn)
{
  int j;

  if( ans && wn )
    for( j=0; j<zVecSizeNC(ans); j++ )
      zVecElemNC(ans,j) *= sqrt( zVecElemNC(wn,j) );
  return ans;
}

/* generalized linear equation solver with MP-inverse based on LQ decomposition (destructive). */
zVec zLESolveMP_LQ_DST(zLEWorkspace *workspace, const zVec we, const zVec wn, int rank, zVec ans)
{
  if( !_zLESolveMP_L( workspace, we ) ) return NULL;
  /* Q^T z is not used for a row-long Q, since Q by Gram=Schmidt's method might lose orthogonality. */
  return _zLESolveMPScaleAns( zMatIsSqr( workspace->r_facto ) ?
    zMulMatTVec( workspace->r_facto, workspace->v_mp_ispace, ans ) :
    zLESolveNormMinDST( workspace->r_facto, workspace->v_mp_ispace, NULL, ans, workspace ), wn );
}

/* generalized linear equation solver using MP-inverse based on LU decomposition (destructive). */
zVec zLESolveMP_LU_DST(zLEWorkspace *workspace, const zVec we, const zVec wn, int rank, zVec ans)
{
  if( !_zLESolveMP_L( workspace, we ) ) return NULL;
  return _zLESolveMPScaleAns( zMatIsSqr( workspace->r_facto ) ?
    zLESolveUCP( workspace->r_facto, workspace->v_mp_ispace, ans, workspace->index_col ) :
    zLESolveNormMinDST( workspace->r_facto, workspace->v_mp_ispace, NULL, ans, workspace ), wn );
}

/* generalized linear equation solver using MP-inverse with the null space (destructive). */
zVec zLESolveMPNullDST(zLEWorkspace *workspace, const zVec we, const zVec wn, int rank, zVec ans, zMat mn)
{
  int i, j;

  if( wn )
    for( j=0; j<zVecSizeNC(wn); j++ )
      if( !( zVecElemNC(wn,j) > 0 ) ){
        ZRUNERROR( ZM_ERR_LE_NONPOSITIVEWEIGHT );
        return NULL;
      }
  if( !zLESolveMP_LQ_DST( workspace, we, wn, rank, ans ) ) return NULL;
  if( zMatIsSqr( workspace->r_facto ) ){
    zMatZero( mn );
    return ans;
  }
  /* mn = a# a - I = Wn^1/2 ( Q^T Q - I ) Wn^-1/2 */
  zMulMatTMat( workspace->r_facto, workspace->r_facto, mn );
  for( i=0; i<zMatRowSizeNC(mn); i++ ){
    zMatElemNC(mn,i,i) -= 1.0;
    if( wn )
      for( j=0; j<zMatColSizeNC(mn); j++ )
        if( i != j ) zMatElemNC(mn,i,j) *= sqrt( zVecElemNC(wn,i) / zVecElemNC(wn,j) );
  }
  return ans;
}

/* allocate workspace for generalized linear equation solvers based on MP inverse, copy the
 * right-hand-side vector, and weight the coefficient matrix. */
static zMat _zLEWorkspaceAllocMPWeighted(zLEWorkspace *workspace, const zMat a, const zVec b, const zVec wn, const zVec we)
{
  zMat aw;

  if( !_zLEWorkspaceAllocMPAndCopyAB( workspace, a, b ) ) return NULL;
  if( !( aw = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a) ) ) ) return NULL;
  return _zLEWeightMat( a, wn, we, aw );
}

/* generalized linear equation solver with MP-inverse based on LQ decomposition. */
zVec zLESolveMP_LQ(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans)
{
  int rank;
  zMat aw;
  zLEWorkspace workspace;
  zVec ret = NULL;

  if( !( aw = _zLEWorkspaceAllocMPWeighted( &workspace, a, b, wn, we ) ) ) goto TERMINATE;
  if( ( rank = zMatDecompLQ( aw, workspace.l_facto, workspace.r_facto, workspace.index_row ) ) > 0 ){
    zLEWorkspaceResizeMP( &workspace, rank );
    ret = zLESolveMP_LQ_DST( &workspace, we, wn, rank, ans );
  } else
  if( rank == 0 )
    ret = zVecZero( ans );
 TERMINATE:
  zMatFree( aw );
  zLEWorkspaceFreeMP( &workspace );
  return ret;
}

/* generalized linear equation solver using MP-inverse based on LU decomposition. */
zVec zLESolveMP_LU(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans)
{
  int rank;
  zMat aw;
  zLEWorkspace workspace;
  zVec ret = NULL;

  if( !( aw = _zLEWorkspaceAllocMPWeighted( &workspace, a, b, wn, we ) ) ) goto TERMINATE;
  if( ( rank = zMatDecompLUCP( aw, workspace.l_facto, workspace.r_facto, workspace.index_row, workspace.index_col ) ) > 0 ){
    zLEWorkspaceResizeMP( &workspace, rank );
    ret = zLESolveMP_LU_DST( &workspace, we, wn, rank, ans );
  } else
  if( rank == 0 )
    ret = zVecZero( ans );
 TERMINATE:
  zMatFree( aw );
  zLEWorkspaceFreeMP( &workspace );
  return ret;
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
  int i, rank;
  zMat aw;
  zLEWorkspace workspace;
  zVec ret = NULL;

  if( !( aw = _zLEWorkspaceAllocMPWeighted( &workspace, a, b, wn, we ) ) ) goto TERMINATE;
  if( ( rank = zMatDecompLQ( aw, workspace.l_facto, workspace.r_facto, workspace.index_row ) ) > 0 ){
    zLEWorkspaceResizeMP( &workspace, rank );
    ret = zLESolveMPNullDST( &workspace, we, wn, rank, ans, mn );
  } else
  if( rank == 0 ){ /* a# = O */
    zMatZero( mn );
    for( i=0; i<zMatRowSizeNC(mn); i++ )
      zMatSetElemNC( mn, i, i, -1.0 );
    ret = zVecZero( ans );
  }
 TERMINATE:
  zMatFree( aw );
  zLEWorkspaceFreeMP( &workspace );
  return ret;
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
  ans = zLESolveMP( a, bb, wn, we, ans );
  zVecFree( bb );
  return ans ? zVecAddNCDRC( ans, aux ) : NULL;
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
