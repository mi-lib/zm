/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_qp_asm - optimization tools: active set method for quadratic programming.
 */

#include <zm/zm_opt.h>

typedef struct{
  zMat qinv;
  zMat aqinvat;
  zVec qinvc;
  zVec aqinvc_plus_b;
  zVec xtmp;
  zIndex ia; /* active set */
  zIndex in; /* non-active set */
  /* workspace */
  zMat _m;
  zVec _v;
  zVec _lambda;
  zIndex _idx;
  zVec _s; /* for balancing */
} zQPASM;

/* initialize working space. */
static bool _zQPASMInit(zQPASM *qpasm, const zMat q, const zVec c, const zMat a, const zVec b)
{
  qpasm->qinv = zMatAllocSqr( zVecSizeNC(c) );
  qpasm->aqinvat = zMatAllocSqr( zMatRowSizeNC(a) );
  qpasm->qinvc = zVecAlloc( zVecSizeNC(c) );
  qpasm->aqinvc_plus_b = zVecAlloc( zVecSizeNC(b) );
  qpasm->xtmp = zVecAlloc( zVecSizeNC(c) );
  qpasm->ia = zIndexAlloc( zMatRowSizeNC(a) );
  qpasm->in = zIndexAlloc( zMatRowSizeNC(a) );

  qpasm->_m = zMatAllocSqr( zMatRowSizeNC(a) );
  qpasm->_v = zVecAlloc( zMatRowSizeNC(a) );
  qpasm->_lambda = zVecAlloc( zMatRowSizeNC(a) );
  qpasm->_idx = zIndexAlloc( zMatRowSizeNC(a) );
  qpasm->_s = zVecAlloc( _zMax( zVecSizeNC(c), zMatRowSizeNC(a) ) ); /* used for both active set and solution */
  if( !qpasm->qinv || !qpasm->aqinvat || !qpasm->qinvc || !qpasm->aqinvc_plus_b ||
      !qpasm->xtmp || !qpasm->ia || !qpasm->in ||
      !qpasm->_m || !qpasm->_v || !qpasm->_lambda || !qpasm->_idx || !qpasm->_s ) return false;

  zMatInv( q, qpasm->qinv );
  zMulMatMatMatTNC( a, qpasm->qinv, qpasm->aqinvat );
  zMulMatVecNC( qpasm->qinv, c, qpasm->qinvc );
  zMulMatVecNC( a, qpasm->qinvc, qpasm->aqinvc_plus_b );
  zVecAddNCDRC( qpasm->aqinvc_plus_b, b );
  return true;
}

/* destroy working space. */
static void _zQPASMDestroy(zQPASM *qpasm)
{
  zMatFree( qpasm->qinv );
  zMatFree( qpasm->aqinvat );
  zVecFree( qpasm->qinvc );
  zVecFree( qpasm->aqinvc_plus_b );
  zVecFree( qpasm->xtmp );
  zIndexFree( qpasm->ia );
  zIndexFree( qpasm->in );
  zMatFree( qpasm->_m );
  zVecFree( qpasm->_v );
  zVecFree( qpasm->_lambda );
  zIndexFree( qpasm->_idx );
  zVecFree( qpasm->_s );
}

/* get initial feasible solution and active set by simplex method. */
static bool _zQPASMInitBase(zQPASM *qpasm, const zMat a, const zVec b)
{
  zLPTableau tab;
  int i, n, i_to_activate;
  bool ret = true;

  tab.a = zMatAlloc( zMatRowSizeNC(a), 2 * ( zMatColSizeNC(a) + zMatRowSizeNC(a) ) );
  tab.b = zVecAlloc( zMatRowSizeNC(a) );
  tab.c = zVecAlloc( zMatColSizeNC(tab.a) );
  tab.ib = zIndexCreate( zMatRowSizeNC(a) );
  tab.in = zIndexCreate( zMatColSizeNC(tab.a) - zMatRowSizeNC(a) );
  tab.ir = zIndexCreate( zMatRowSizeNC(a) );
  if( !tab.a || !tab.b ||!tab.c || !tab.ib || !tab.in || !tab.ir ){
    zLPTableauDestroy( &tab );
    return false;
  }
  n = 2 * zMatColSizeNC(a);
  for( i=0; i<zVecSizeNC(b); i++ ){
    if( zVecElemNC(b,i) >= 0 ){
      zRawVecCopy( zMatRowBufNC(a,i), zMatRowBufNC(tab.a,i), zMatColSizeNC(a) );
      zRawVecRev( zMatRowBufNC(a,i), zMatRowBufNC(tab.a,i)+zMatColSize(a), zMatColSizeNC(a) );
      zMatSetElemNC( tab.a, i, n+i,-1.0 );
      zVecSetElemNC( tab.b, i, zVecElemNC(b,i) );
    } else{
      zRawVecRev( zMatRowBufNC(a,i), zMatRowBufNC(tab.a,i), zMatColSizeNC(a) );
      zRawVecCopy( zMatRowBufNC(a,i), zMatRowBufNC(tab.a,i)+zMatColSize(a), zMatColSizeNC(a) );
      zMatSetElemNC( tab.a, i, n+i, 1.0 );
      zVecSetElemNC( tab.b, i,-zVecElemNC(b,i) );
    }
    zMatSetElemNC( tab.a, i, zIndexSizeNC(tab.in)+i, 1.0 );
  }
  for( i=zIndexSizeNC(tab.in); i<zVecSizeNC(tab.c); i++ )
    zVecSetElemNC( tab.c, i, 1.0 );
  tab.d = 0;
  zIndexOrder( tab.ib, zIndexSizeNC(tab.in) );
  zIndexOrder( tab.in, 0 );
  zIndexOrder( tab.ir, 0 );
  if( !zLPTableauSimplex( &tab ) || !zIsTiny(tab.d) ){
    ZRUNWARN( ZM_ERR_OPT_UNSOLVABLE );
    ret = false;
    goto TERMINATE;
  }
  zLPTableauFindBase( &tab ); /* remove artificial variables from the bases */
  /* rearrange initial feasible solution and active set */
  zIndexSetSize( qpasm->ia, zIndexSizeNC(tab.ir) );
  zIndexCopyNC( tab.ir, qpasm->ia );
  zIndexSetSize( qpasm->in, 0 );
  zVecZero( qpasm->xtmp );
  for( i=0; i<zIndexSizeNC(tab.ib); i++ ){
    if( zIndexElemNC(tab.ib,i) < zMatColSizeNC(a) ){
      zVecSetElemNC( qpasm->xtmp, zIndexElemNC(tab.ib,i), zVecElemNC(tab.b,i) );
    } else
    if( ( i_to_activate = zIndexElemNC(tab.ib,i) - n ) < 0 ){
      zVecSetElemNC( qpasm->xtmp, zIndexElemNC(tab.ib,i)-zMatColSizeNC(a),-zVecElemNC(tab.b,i) );
    } else
    if( i_to_activate < zMatRowSizeNC(a) ){
      zIndexRemoveVal( qpasm->ia, i_to_activate );
      zIndexInsertVal( qpasm->in, i_to_activate );
    } else{
      ZRUNERROR( ZM_ERR_FATAL );
      ret = false;
      goto TERMINATE;
    }
  }
 TERMINATE:
  zLPTableauDestroy( &tab );
  return ret;
}

/* solve an equation A_a Q^-1 A_a^T lambda = v with respect to lambda for active set a,
 * and find dx = Q^-1 A_a^T lambda, where v has to be preset in the workspace. */
static void _zQPASMSolveActive(zQPASM *qpasm, const zMat a, zVec dx)
{
  int i, j;

  zMatSetSizeNC( qpasm->_m, zIndexSizeNC(qpasm->ia), zIndexSizeNC(qpasm->ia) );
  zVecSetSizeNC( qpasm->_lambda, zIndexSizeNC(qpasm->ia) );
  zVecSetSizeNC( qpasm->_s, zIndexSizeNC(qpasm->ia) );
  zIndexSetSize( qpasm->_idx, zIndexSizeNC(qpasm->ia) );
  zIndexOrder( qpasm->_idx, 0 );
  for( i=0; i<zIndexSizeNC(qpasm->ia); i++ )
    for( j=0; j<zIndexSizeNC(qpasm->ia); j++ )
      zMatElemNC(qpasm->_m,i,j) = zMatElemNC(qpasm->aqinvat,zIndexElemNC(qpasm->ia,i),zIndexElemNC(qpasm->ia,j));
  zLESolveGaussDST( qpasm->_m, qpasm->_v, qpasm->_lambda, qpasm->_idx, qpasm->_s );
  zVecSetSizeNC( qpasm->_s, zVecSizeNC(dx) );
  for( i=0; i<zVecSizeNC(qpasm->_s); i++ ){
    zVecSetElemNC( qpasm->_s, i, 0 );
    for( j=0; j<zIndexSizeNC(qpasm->ia); j++ )
      zVecElemNC(qpasm->_s,i) += zMatElemNC(a,zIndexElemNC(qpasm->ia,j),i) * zVecElemNC(qpasm->_lambda,j);
  }
  zMulMatVecNC( qpasm->qinv, qpasm->_s, dx );
}

/* solve an equation to find temporary solution and adjoint variables. */
static bool _zQPASMSolveEq(zQPASM *qpasm, const zMat a, const zVec b, zVec x)
{
  int i;

  zVecSetSizeNC( qpasm->_v, zIndexSizeNC(qpasm->ia) );
  for( i=0; i<zIndexSizeNC(qpasm->ia); i++ )
    zVecElemNC(qpasm->_v,i) = zVecElemNC(qpasm->aqinvc_plus_b,zIndexElemNC(qpasm->ia,i));
  _zQPASMSolveActive( qpasm, a, x );
  zVecSubNCDRC( x, qpasm->qinvc );
  return true;
}

/* refine the solution by one step if it does not satisfy the active constraints accurately
 * due to the round-off error in an ill-conditioned case. Note that the adjoint variables
 * and the current solution in the workspace are overwritten. */
static void _zQPASMRefine(zQPASM *qpasm, const zMat a, const zVec b, zVec x)
{
  int i;

  if( zIndexSizeNC(qpasm->ia) == 0 ) return;
  zVecSetSizeNC( qpasm->_v, zIndexSizeNC(qpasm->ia) );
  for( i=0; i<zIndexSizeNC(qpasm->ia); i++ ) /* residual of active constraints */
    zVecElemNC(qpasm->_v,i) = zVecElemNC(b,zIndexElemNC(qpasm->ia,i))
      - zRawVecInnerProd( zMatRowBufNC(a,zIndexElemNC(qpasm->ia,i)), zVecBufNC(x), zVecSizeNC(x) );
  if( zVecNorm( qpasm->_v ) <= zTOL ) return;
  _zQPASMSolveActive( qpasm, a, qpasm->xtmp );
  zVecAddNCDRC( x, qpasm->xtmp );
}

/* check if a blocking constraint i is more transversal to the step than another constraint j,
 * where api and apj are the inner products of the constraint vectors and the step. */
static bool _zQPASMIsMoreTransversal(const zMat a, int i, double api, int j, double apj)
{
  return api / zRawVecNorm( zMatRowBufNC(a,i), zMatColSizeNC(a) )
       < apj / zRawVecNorm( zMatRowBufNC(a,j), zMatColSizeNC(a) );
}

/* move the current feasible solution toward the temporary solution, and add a blocking
 * constraint to active set if any. If no constraint blocks, the current solution is
 * replaced with the temporary solution and the false value is returned. */
static bool _zQPASMAddEq(zQPASM *qpasm, const zMat a, const zVec b, const zVec x)
{
  int i, imin;
  double ax, axtmp, d, dmin, axmin = 0;
  bool is_violated = false;

  /* the temporary solution coincides with the current solution if active set
   * determines the solution uniquely. */
  if( zIndexSizeNC(qpasm->ia) >= zVecSizeNC(x) ) return false;
  for( dmin=HUGE_VAL, imin=-1, i=0; i<zIndexSizeNC(qpasm->in); i++ ){
    ax = zRawVecInnerProd( zMatRowBufNC(a,zIndexElemNC(qpasm->in,i)), zVecBufNC(x), zVecSizeNC(x) );
    if( zVecElemNC(b,zIndexElemNC(qpasm->in,i)) <= ax ) continue;
    ax -= ( axtmp = zRawVecInnerProd( zMatRowBufNC(a,zIndexElemNC(qpasm->in,i)), zVecBufNC(qpasm->xtmp), zVecSizeNC(qpasm->xtmp) ) );
    if( ax > -zTOL ) continue; /* not a blocking constraint */
    is_violated = true;
    if( ( d = ( zVecElemNC(b,zIndexElemNC(qpasm->in,i)) - axtmp ) / ax ) < 0 ) d = 0;
    if( d < dmin - zTOL ){
      dmin = d;
      imin = zIndexElemNC(qpasm->in,i);
      axmin = ax;
    } else
    if( d < dmin + zTOL ){ /* tie (degenerate case): prefer a constraint more transversal to the step */
      if( _zQPASMIsMoreTransversal( a, zIndexElemNC(qpasm->in,i), ax, imin, axmin ) ){
        imin = zIndexElemNC(qpasm->in,i);
        axmin = ax;
      }
      if( d < dmin ) dmin = d;
    }
  }
  if( !is_violated ){ /* full step */
    zVecCopyNC( x, qpasm->xtmp );
    return false;
  }
  zIndexRemoveVal( qpasm->in, imin );
  zIndexInsertVal( qpasm->ia, imin );
  zVecInterDivDRC( qpasm->xtmp, x, dmin );
  return true;
}

/* delete an equality constraint with the most negative Lagrange multiplier from active set. */
static bool _zQPASMDelEq(zQPASM *qpasm, const zMat a)
{
  int i, imin;
  double lambda_min;

  for( lambda_min=0, imin=-1, i=0; i<zVecSizeNC(qpasm->_lambda); i++ )
    if( zVecElemNC(qpasm->_lambda,i) < lambda_min ){
      lambda_min = zVecElemNC(qpasm->_lambda,i);
      imin = i;
    }
  if( imin < 0 ) return false;
  zIndexInsertVal( qpasm->in, zIndexElemNC(qpasm->ia,imin) );
  zIndexRemove( qpasm->ia, imin );
  return true;
}

/* check if vector/matrix sizes are consistent and the quadratic term matrix is positive-definite. */
static bool _zQPASMCheck(const zMat q, const zVec c, const zMat a, const zVec b, zVec ans)
{
  zMat l;
  zIndex ic;
  bool ret;

  if( !zMatIsSymmetric( q ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSYMMETRIC );
    return false;
  }
  if( !zMatRowVecSizeEqual( q, c ) || !zVecSizeEqual( c, ans ) ||
      !zMatRowVecSizeEqual( a, b ) || !zMatColVecSizeEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return false;
  }
  if( !( ret = !( zMatDecompCholeskyAlloc( q, &l, &ic ) < zMatRowSizeNC(q) ) ) )
    ZRUNERROR( ZM_ERR_OPT_NONCONVEX );
  zMatFree( l );
  zIndexFree( ic );
  return ret;
}

/* solve quadratic programming by active set method. */
bool zQPSolveASM(const zMat q, const zVec c, const zMat a, const zVec b, zVec ans, double *cost)
{
  zQPASM qpasm;
  bool ret = false;

  if( !_zQPASMCheck( q, c, a, b, ans ) ) return false;
  if( !_zQPASMInit( &qpasm, q, c, a, b ) ) goto TERMINATE;
  if( !_zQPASMInitBase( &qpasm, a, b ) ) goto TERMINATE;
  do{
    do{
      _zQPASMSolveEq( &qpasm, a, b, ans );
    } while( _zQPASMAddEq( &qpasm, a, b, ans ) );
  } while( _zQPASMDelEq( &qpasm, a ) );
  _zQPASMRefine( &qpasm, a, b, ans );
  if( cost ) *cost = zQuadraticValue( q, c, ans );
  ret = true;
 TERMINATE:
  _zQPASMDestroy( &qpasm );
  return ret;
}
