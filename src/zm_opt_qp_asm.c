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
static bool _zQPASMInit(zQPASM *asm, zMat q, zVec c, zMat a, zVec b)
{
  asm->qinv = zMatAllocSqr( zVecSizeNC(c) );
  asm->aqinvat = zMatAllocSqr( zMatRowSizeNC(a) );
  asm->qinvc = zVecAlloc( zVecSizeNC(c) );
  asm->aqinvc_plus_b = zVecAlloc( zVecSizeNC(b) );
  asm->xtmp = zVecAlloc( zVecSizeNC(c) );
  asm->ia = zIndexAlloc( zMatRowSizeNC(a) );
  asm->in = zIndexAlloc( zMatRowSizeNC(a) );

  asm->_m = zMatAllocSqr( zMatRowSizeNC(a) );
  asm->_v = zVecAlloc( zMatRowSizeNC(a) );
  asm->_lambda = zVecAlloc( zMatRowSizeNC(a) );
  asm->_idx = zIndexAlloc( zMatRowSizeNC(a) );
  asm->_s = zVecAlloc( zMatRowSizeNC(a) );
  if( !asm->qinv || !asm->aqinvat || !asm->qinvc || !asm->aqinvc_plus_b ||
      !asm->xtmp || !asm->ia || !asm->in ||
      !asm->_m || !asm->_v || !asm->_lambda || !asm->_idx || !asm->_s ) return false;

  zMatInv( q, asm->qinv );
  zMulMatMatMatTNC( a, asm->qinv, asm->aqinvat );
  zMulMatVecNC( asm->qinv, c, asm->qinvc );
  zMulMatVecNC( a, asm->qinvc, asm->aqinvc_plus_b );
  zVecAddNCDRC( asm->aqinvc_plus_b, b );
  return true;
}

/* destroy working space. */
static void _zQPASMDestroy(zQPASM *asm)
{
  zMatFree( asm->qinv );
  zMatFree( asm->aqinvat );
  zVecFree( asm->qinvc );
  zVecFree( asm->aqinvc_plus_b );
  zVecFree( asm->xtmp );
  zIndexFree( asm->ia );
  zIndexFree( asm->in );
  zMatFree( asm->_m );
  zVecFree( asm->_v );
  zVecFree( asm->_lambda );
  zIndexFree( asm->_idx );
  zVecFree( asm->_s );
}

/* get initial feasible solution and active set by simplex method. */
static bool _zQPASMInitBase(zQPASM *asm, zMat a, zVec b)
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
    ZRUNWARN( ZM_ERR_OPT_UNSOLVE );
    ret = false;
    goto TERMINATE;
  }
  zLPTableauFindBase( &tab ); /* remove artificial variables from the bases */
  /* rearrange initial feasible solution and active set */
  zIndexSizeNC(asm->ia) = zIndexSizeNC(tab.ir);
  zIndexCopyNC( tab.ir, asm->ia );
  zIndexSizeNC(asm->in) = 0;
  zVecZero( asm->xtmp );
  for( i=0; i<zArraySize(tab.ib); i++ ){
    if( zIndexElemNC(tab.ib,i) < zMatColSizeNC(a) ){
      zVecSetElemNC( asm->xtmp, zIndexElemNC(tab.ib,i), zVecElemNC(tab.b,i) );
    } else
    if( zIndexElemNC(tab.ib,i) < n ){
      zVecSetElemNC( asm->xtmp, zIndexElemNC(tab.ib,i)-zMatColSizeNC(a),-zVecElemNC(tab.b,i) );
    } else
    if( zIndexElemNC(tab.ib,i) < n + zMatRowSizeNC(a) ){
      zIndexRemoveVal( asm->ia, zIndexElemNC(tab.ib,i)-n );
      zIndexInsertVal( asm->in, zMatRowSizeNC(a), zIndexElemNC(tab.ib,i)-n );
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
static bool _zQPASMSolveEq(zQPASM *asm, zMat a, zVec b, zVec x)
{
  int i, j;

  zMatSetSize( asm->_m, zIndexSizeNC(asm->ia), zIndexSizeNC(asm->ia) );
  zVecSetSize( asm->_v, zIndexSizeNC(asm->ia) );
  zVecSetSize( asm->_lambda, zIndexSizeNC(asm->ia) );
  zArraySize( asm->_idx ) = zIndexSizeNC(asm->ia);
  zVecSetSize( asm->_s, zIndexSizeNC(asm->ia) );
  zIndexOrder( asm->_idx, 0 );
  for( i=0; i<zIndexSizeNC(asm->ia); i++ ){
    for( j=0; j<zIndexSizeNC(asm->ia); j++ ){
      zMatElemNC(asm->_m,i,j) = zMatElemNC(asm->aqinvat,zIndexElemNC(asm->ia,i),zIndexElemNC(asm->ia,j));
    }
    zVecElemNC(asm->_v,i) = zVecElemNC(asm->aqinvc_plus_b,zIndexElemNC(asm->ia,i));
  }
  zLESolveGaussDST( asm->_m, asm->_v, asm->_lambda, asm->_idx, asm->_s );
  zVecSetSize( asm->_s, zVecSizeNC(x) );
  for( i=0; i<zVecSizeNC(asm->_s); i++ ){
    zVecElemNC(asm->_s,i) = 0;
    for( j=0; j<zIndexSizeNC(asm->ia); j++ )
      zVecElemNC(asm->_s,i) += zMatElemNC(a,zIndexElemNC(asm->ia,j),i) * zVecElemNC(asm->_lambda,j);
  }
  zMulMatVecNC( asm->qinv, asm->_s, x );
  zVecSubNCDRC( x, asm->qinvc );
  return true;
}

/* add an equality constraint to active set. */
static bool _zQPASMAddEq(zQPASM *asm, zMat a, zVec b, zVec x)
{
  int i, imin;
  double ax, axtmp, d, dmin;
  bool is_violated = false;

  for( dmin=HUGE_VAL, imin=-1, i=0; i<zIndexSizeNC(asm->in); i++ ){
    ax = zRawVecInnerProd( zMatRowBufNC(a,zIndexElemNC(asm->in,i)), zVecBufNC(x), zVecSizeNC(x) );
    if( zVecElemNC(b,zIndexElemNC(asm->in,i)) <= ax ) continue;
    ax -= ( axtmp = zRawVecInnerProd( zMatRowBufNC(a,zIndexElemNC(asm->in,i)), zVecBufNC(asm->xtmp), zVecSizeNC(asm->xtmp) ) );
    if( zIsTiny( ax ) ) continue;
    is_violated = true;
    if( ( d = ( zVecElemNC(b,zIndexElemNC(asm->in,i)) - axtmp ) / ax ) < dmin ){
      dmin = d;
      imin = zIndexElemNC(asm->in,i);
    }
  }
  if( !is_violated ) return false;
  zIndexRemoveVal( asm->in, imin );
  zIndexInsertVal( asm->ia, zMatRowSizeNC(a), imin );
  zVecInterDivDRC( asm->xtmp, x, dmin );
  return true;
}

/* delete an equality constraint from active set. */
static bool _zQPASMDelEq(zQPASM *asm, zMat a)
{
  int i, ia;
  bool is_violated = false;

  for( i=0; i<zVecSizeNC(asm->_lambda); i++ ){
    if( zVecElemNC(asm->_lambda,i) < 0 ){
      is_violated = true;
      ia = zIndexElemNC(asm->ia,i);
      zIndexRemove( asm->ia, i );
      zIndexInsertVal( asm->in, zMatRowSizeNC(a), ia );
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
    ZRUNERROR( ZM_ERR_NONSYMMETRIC_MAT );
    return false;
  }
  if( !zMatRowVecSizeIsEqual( q, c ) || !zVecSizeIsEqual( c, ans ) ||
      !zMatRowVecSizeIsEqual( a, b ) || !zMatColVecSizeIsEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
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
  zQPASM asm;
  bool ret = false;

  if( !_zQPASMCheck( q, c, a, b, ans ) ) return false;
  if( !_zQPASMInit( &asm, q, c, a, b ) ) goto TERMINATE;
  if( !_zQPASMInitBase( &asm, a, b ) ) goto TERMINATE;
  do{
    do{
      _zQPASMSolveEq( &asm, a, b, ans );
    } while( _zQPASMAddEq( &asm, a, b, ans ) );
  } while( _zQPASMDelEq( &asm, a ) );
  if( cost ) *cost = zQuadraticValue( q, c, ans );
  ret = true;
 TERMINATE:
  _zQPASMDestroy( &asm );
  return ret;
}
