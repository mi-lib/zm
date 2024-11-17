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
static bool _zQPASMInit(zQPASM *qpasm, zMat q, zVec c, zMat a, zVec b)
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
  qpasm->_s = zVecAlloc( zMatRowSizeNC(a) );
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
static bool _zQPASMInitBase(zQPASM *qpasm, zMat a, zVec b)
{
  zLPTableau tab;
  int i, n;
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
    zMatSetElemNC( tab.a, i, zArraySize(tab.in)+i, 1.0 );
  }
  for( i=zArraySize(tab.in); i<zVecSizeNC(tab.c); i++ )
    zVecSetElemNC( tab.c, i, 1.0 );
  tab.d = 0;
  zIndexOrder( tab.ib, zArraySize(tab.in) );
  zIndexOrder( tab.in, 0 );
  zIndexOrder( tab.ir, 0 );
  if( !zLPTableauSimplex( &tab ) || !zIsTiny(tab.d) ){
    ZRUNWARN( ZM_ERR_OPT_UNSOLVABLE );
    ret = false;
    goto TERMINATE;
  }
  zLPTableauFindBase( &tab ); /* remove artificial variables from the bases */
  /* rearrange initial feasible solution and active set */
  zIndexSizeNC(qpasm->ia) = zIndexSizeNC(tab.ir);
  zIndexCopyNC( tab.ir, qpasm->ia );
  zIndexSizeNC(qpasm->in) = 0;
  zVecZero( qpasm->xtmp );
  for( i=0; i<zArraySize(tab.ib); i++ ){
    if( zIndexElemNC(tab.ib,i) < zMatColSizeNC(a) ){
      zVecSetElemNC( qpasm->xtmp, zIndexElemNC(tab.ib,i), zVecElemNC(tab.b,i) );
    } else
    if( zIndexElemNC(tab.ib,i) < n ){
      zVecSetElemNC( qpasm->xtmp, zIndexElemNC(tab.ib,i)-zMatColSizeNC(a),-zVecElemNC(tab.b,i) );
    } else
    if( zIndexElemNC(tab.ib,i) < n + zMatRowSizeNC(a) ){
      zIndexRemoveVal( qpasm->ia, zIndexElemNC(tab.ib,i)-n );
      zIndexInsertVal( qpasm->in, zMatRowSizeNC(a), zIndexElemNC(tab.ib,i)-n );
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

/* solve an equation to find temporary solution and adjoint variables. */
static bool _zQPASMSolveEq(zQPASM *qpasm, zMat a, zVec b, zVec x)
{
  int i, j;

  zMatSetSize( qpasm->_m, zIndexSizeNC(qpasm->ia), zIndexSizeNC(qpasm->ia) );
  zVecSetSize( qpasm->_v, zIndexSizeNC(qpasm->ia) );
  zVecSetSize( qpasm->_lambda, zIndexSizeNC(qpasm->ia) );
  zArraySize( qpasm->_idx ) = zIndexSizeNC(qpasm->ia);
  zVecSetSize( qpasm->_s, zIndexSizeNC(qpasm->ia) );
  zIndexOrder( qpasm->_idx, 0 );
  for( i=0; i<zIndexSizeNC(qpasm->ia); i++ ){
    for( j=0; j<zIndexSizeNC(qpasm->ia); j++ ){
      zMatElemNC(qpasm->_m,i,j) = zMatElemNC(qpasm->aqinvat,zIndexElemNC(qpasm->ia,i),zIndexElemNC(qpasm->ia,j));
    }
    zVecElemNC(qpasm->_v,i) = zVecElemNC(qpasm->aqinvc_plus_b,zIndexElemNC(qpasm->ia,i));
  }
  zLESolveGaussDST( qpasm->_m, qpasm->_v, qpasm->_lambda, qpasm->_idx, qpasm->_s );
  zVecSetSize( qpasm->_s, zVecSizeNC(x) );
  for( i=0; i<zVecSizeNC(qpasm->_s); i++ ){
    zVecElemNC(qpasm->_s,i) = 0;
    for( j=0; j<zIndexSizeNC(qpasm->ia); j++ )
      zVecElemNC(qpasm->_s,i) += zMatElemNC(a,zIndexElemNC(qpasm->ia,j),i) * zVecElemNC(qpasm->_lambda,j);
  }
  zMulMatVecNC( qpasm->qinv, qpasm->_s, x );
  zVecSubNCDRC( x, qpasm->qinvc );
  return true;
}

/* add an equality constraint to active set. */
static bool _zQPASMAddEq(zQPASM *qpasm, zMat a, zVec b, zVec x)
{
  int i, imin;
  double ax, axtmp, d, dmin;
  bool is_violated = false;

  for( dmin=HUGE_VAL, imin=-1, i=0; i<zIndexSizeNC(qpasm->in); i++ ){
    ax = zRawVecInnerProd( zMatRowBufNC(a,zIndexElemNC(qpasm->in,i)), zVecBufNC(x), zVecSizeNC(x) );
    if( zVecElemNC(b,zIndexElemNC(qpasm->in,i)) <= ax ) continue;
    ax -= ( axtmp = zRawVecInnerProd( zMatRowBufNC(a,zIndexElemNC(qpasm->in,i)), zVecBufNC(qpasm->xtmp), zVecSizeNC(qpasm->xtmp) ) );
    if( zIsTiny( ax ) ) continue;
    is_violated = true;
    if( ( d = ( zVecElemNC(b,zIndexElemNC(qpasm->in,i)) - axtmp ) / ax ) < dmin ){
      dmin = d;
      imin = zIndexElemNC(qpasm->in,i);
    }
  }
  if( !is_violated ) return false;
  zIndexRemoveVal( qpasm->in, imin );
  zIndexInsertVal( qpasm->ia, zMatRowSizeNC(a), imin );
  zVecInterDivDRC( qpasm->xtmp, x, dmin );
  return true;
}

/* delete an equality constraint from active set. */
static bool _zQPASMDelEq(zQPASM *qpasm, zMat a)
{
  int i, ia;
  bool is_violated = false;

  for( i=0; i<zVecSizeNC(qpasm->_lambda); i++ ){
    if( zVecElemNC(qpasm->_lambda,i) < 0 ){
      is_violated = true;
      ia = zIndexElemNC(qpasm->ia,i);
      zIndexRemove( qpasm->ia, i );
      zIndexInsertVal( qpasm->in, zMatRowSizeNC(a), ia );
    }
  }
  return is_violated;
}

/* check if vector/matrix sizes are consistent and the quadratic term matrix is positive-definite. */
static bool _zQPASMCheck(zMat q, zVec c, zMat a, zVec b, zVec ans)
{
  zMat l;
  zIndex ic;
  bool ret;

  if( !zMatIsSymmetric( q ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSYMMETRIC );
    return false;
  }
  if( !zMatRowVecSizeIsEqual( q, c ) || !zVecSizeIsEqual( c, ans ) ||
      !zMatRowVecSizeIsEqual( a, b ) || !zMatColVecSizeIsEqual( a, ans ) ){
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
bool zQPSolveASM(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost)
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
  if( cost ) *cost = zQuadraticValue( q, c, ans );
  ret = true;
 TERMINATE:
  _zQPASMDestroy( &qpasm );
  return ret;
}
