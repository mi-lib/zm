/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_qp_asm - optimization tools: active set method for quadratic programming.
 */

#include <zm/zm_opt.h>

/* The problem is transformed with the Cholesky factor L of Q (Q = L L^T) and z = L^T x into
 *   0.5 |z|^2 + cbar^T z -> min. subject to abar z >= b,
 * where abar = A L^-T and cbar = L^-1 c. The coefficient matrix of the active constraints is decomposed as
 *   abar_W = [ L_W O ] U,
 * where L_W is lower triangular and U is orthogonal, which is updated when a constraint is added
 * (Householder's reflection) or deleted (Givens's rotations) instead of decomposed from scratch.
 * (P. E. Gill, G. H. Golub, W. Murray and M. A. Saunders, Methods for Modifying Matrix Factorizations,
 *  Mathematics of Computation, 28(126), pp. 505-535, 1974.)
 */
typedef struct{
  zMat l;         /* Cholesky factor of Q (rows permuted by il) */
  zIndex il;
  zMat abar;      /* A L^-T */
  zVec cbar;      /* L^-1 c */
  zVec z;         /* current feasible solution in the transformed space */
  zMat lw;        /* lower triangular factor of the active constraints */
  zMat u;         /* orthogonal factor of the active constraints */
  zIndex ia;      /* active constraints corresponding to the rows of lw */
  bool *is_active;
  /* workspace */
  zVec _t;
  zVec _y;
  zVec _lambda;
  zVec _v;
} zQPASM;

/* out = L^-1 v. */
static void _zQPASMCholeskyFwdSubst(const zQPASM *qpasm, const double *v, double *out)
{
  int i, j;
  double s;

  for( i=0; i<zIndexSizeNC(qpasm->il); i++ ){
    s = v[zIndexElemNC(qpasm->il,i)];
    for( j=0; j<i; j++ ) s -= zMatElemNC(qpasm->l,zIndexElemNC(qpasm->il,i),j) * out[j];
    out[i] = s / zMatElemNC(qpasm->l,zIndexElemNC(qpasm->il,i),i);
  }
}

/* x = L^-T z, where w is a workspace. */
static void _zQPASMCholeskyBwdSubst(const zQPASM *qpasm, const zVec z, zVec x, double *w)
{
  int i, j;
  double s;

  for( i=zIndexSizeNC(qpasm->il)-1; i>=0; i-- ){
    s = zVecElemNC(z,i);
    for( j=i+1; j<zIndexSizeNC(qpasm->il); j++ ) s -= zMatElemNC(qpasm->l,zIndexElemNC(qpasm->il,j),i) * w[j];
    w[i] = s / zMatElemNC(qpasm->l,zIndexElemNC(qpasm->il,i),i);
  }
  for( i=0; i<zIndexSizeNC(qpasm->il); i++ )
    zVecSetElemNC( x, zIndexElemNC(qpasm->il,i), w[i] );
}

/* initialize working space, and transform the problem. */
static bool _zQPASMInit(zQPASM *qpasm, const zMat q, const zVec c, const zMat a)
{
  int i, n, m;

  n = zVecSizeNC(c);
  m = zMatRowSizeNC(a);
  qpasm->l = zMatAllocSqr( n );
  qpasm->il = zIndexCreate( n );
  qpasm->abar = zMatAlloc( m, n );
  qpasm->cbar = zVecAlloc( n );
  qpasm->z = zVecAlloc( n );
  qpasm->lw = zMatAllocSqr( n );
  qpasm->u = zMatAllocSqr( n );
  qpasm->ia = zIndexAlloc( n );
  qpasm->is_active = zAlloc( bool, m );
  qpasm->_t = zVecAlloc( n );
  qpasm->_y = zVecAlloc( n );
  qpasm->_lambda = zVecAlloc( n );
  qpasm->_v = zVecAlloc( n );
  if( !qpasm->l || !qpasm->il || !qpasm->abar || !qpasm->cbar || !qpasm->z ||
      !qpasm->lw || !qpasm->u || !qpasm->ia || ( m > 0 && !qpasm->is_active ) ||
      !qpasm->_t || !qpasm->_y || !qpasm->_lambda || !qpasm->_v ) return false;
  if( zMatDecompCholesky( q, qpasm->l, qpasm->il ) < n ){
    ZRUNERROR( ZM_ERR_OPT_NONCONVEX );
    return false;
  }
  for( i=0; i<m; i++ ){
    _zQPASMCholeskyFwdSubst( qpasm, zMatRowBufNC(a,i), zMatRowBufNC(qpasm->abar,i) );
    qpasm->is_active[i] = false;
  }
  _zQPASMCholeskyFwdSubst( qpasm, zVecBufNC(c), zVecBufNC(qpasm->cbar) );
  zMatZero( qpasm->lw );
  zMatIdent( qpasm->u );
  zIndexSetSize( qpasm->ia, 0 );
  return true;
}

/* destroy working space. */
static void _zQPASMDestroy(zQPASM *qpasm)
{
  zMatFree( qpasm->l );
  zIndexFree( qpasm->il );
  zMatFree( qpasm->abar );
  zVecFree( qpasm->cbar );
  zVecFree( qpasm->z );
  zMatFree( qpasm->lw );
  zMatFree( qpasm->u );
  zIndexFree( qpasm->ia );
  zFree( qpasm->is_active );
  zVecFree( qpasm->_t );
  zVecFree( qpasm->_y );
  zVecFree( qpasm->_lambda );
  zVecFree( qpasm->_v );
}

/* add a constraint to active set, where the factors are updated by Householder's reflection.
 * If the constraint is linearly dependent on the active constraints, it is not added. */
static bool _zQPASMAddCon(zQPASM *qpasm, int i)
{
  int j, r, n, k;
  double *ai, *w, *v, s, sigma, vv;

  n = zVecSizeNC(qpasm->z);
  if( ( k = zIndexSizeNC(qpasm->ia) ) >= n ) return false;
  ai = zMatRowBufNC(qpasm->abar,i);
  w = zVecBufNC(qpasm->_t);
  v = zVecBufNC(qpasm->_v);
  for( r=0; r<n; r++ ) /* w = U ai */
    w[r] = zRawVecInnerProd( zMatRowBufNC(qpasm->u,r), ai, n );
  if( ( s = zRawVecNorm( w+k, n-k ) ) <= zTOL * zRawVecNorm( ai, n ) )
    return false; /* linearly dependent */
  /* the reflection to map the lower part of w onto the k-th axis is applied to the lower rows of U. */
  sigma = w[k] > 0 ? -s : s;
  zRawVecCopy( w+k, v+k, n-k );
  v[k] -= sigma;
  vv = zRawVecSqrNorm( v+k, n-k );
  for( j=0; j<n; j++ ){
    for( s=0, r=k; r<n; r++ ) s += v[r] * zMatElemNC(qpasm->u,r,j);
    s *= 2 / vv;
    for( r=k; r<n; r++ ) zMatElemNC(qpasm->u,r,j) -= s * v[r];
  }
  zRawVecCopy( w, zMatRowBufNC(qpasm->lw,k), k );
  zMatSetElemNC( qpasm->lw, k, k, sigma );
  zIndexSetSize( qpasm->ia, k+1 );
  zIndexSetElemNC( qpasm->ia, k, i );
  qpasm->is_active[i] = true;
  return true;
}

/* delete the j-th active constraint, where the factors are updated by Givens's rotations. */
static void _zQPASMDelCon(zQPASM *qpasm, int j)
{
  int i, r, n, k;
  double a, b, h, cs, sn, t1, t2;

  n = zVecSizeNC(qpasm->z);
  k = zIndexSizeNC(qpasm->ia);
  qpasm->is_active[zIndexElemNC(qpasm->ia,j)] = false;
  for( i=j; i<k-1; i++ ){
    zRawVecCopy( zMatRowBufNC(qpasm->lw,i+1), zMatRowBufNC(qpasm->lw,i), i+2 );
    zIndexSetElemNC( qpasm->ia, i, zIndexElemNC(qpasm->ia,i+1) );
  }
  zRawVecZero( zMatRowBufNC(qpasm->lw,k-1), n );
  zIndexSetSize( qpasm->ia, k-1 );
  /* each of the rows from the j-th has a non-zero component next to the diagonal, which is
   * eliminated by a rotation of the columns applied to the same rows of U. */
  for( i=j; i<k-1; i++ ){
    a = zMatElemNC(qpasm->lw,i,i);
    if( ( b = zMatElemNC(qpasm->lw,i,i+1) ) == 0 ) continue;
    h = sqrt( a*a + b*b );
    cs = a / h;
    sn = b / h;
    for( r=i; r<k-1; r++ ){
      t1 = zMatElemNC(qpasm->lw,r,i);
      t2 = zMatElemNC(qpasm->lw,r,i+1);
      zMatSetElemNC( qpasm->lw, r, i,   cs*t1 + sn*t2 );
      zMatSetElemNC( qpasm->lw, r, i+1,-sn*t1 + cs*t2 );
    }
    zMatSetElemNC( qpasm->lw, i, i+1, 0 );
    for( r=0; r<n; r++ ){
      t1 = zMatElemNC(qpasm->u,i,r);
      t2 = zMatElemNC(qpasm->u,i+1,r);
      zMatSetElemNC( qpasm->u, i,   r, cs*t1 + sn*t2 );
      zMatSetElemNC( qpasm->u, i+1, r,-sn*t1 + cs*t2 );
    }
  }
}

/* get initial feasible solution and active set by simplex method. */
static bool _zQPASMInitBase(zQPASM *qpasm, const zMat a, const zVec b, zVec x)
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
  if( !zLPTableauFindFeasibleBase( &tab ) ){
    ret = false;
    goto TERMINATE;
  }
  /* rearrange initial feasible solution and active set, where the constraints whose slack variables
   * are not in the bases are active. */
  for( i=0; i<zIndexSizeNC(tab.ir); i++ )
    qpasm->is_active[zIndexElemNC(tab.ir,i)] = true;
  zVecZero( x );
  for( i=0; i<zIndexSizeNC(tab.ib); i++ ){
    if( zIndexElemNC(tab.ib,i) < zMatColSizeNC(a) ){
      zVecSetElemNC( x, zIndexElemNC(tab.ib,i), zVecElemNC(tab.b,zIndexElemNC(tab.ir,i)) );
    } else
    if( ( i_to_activate = zIndexElemNC(tab.ib,i) - n ) < 0 ){
      zVecSetElemNC( x, zIndexElemNC(tab.ib,i)-zMatColSizeNC(a),-zVecElemNC(tab.b,zIndexElemNC(tab.ir,i)) );
    } else
    if( i_to_activate < zMatRowSizeNC(a) ){
      qpasm->is_active[i_to_activate] = false;
    } else{
      ZRUNERROR( ZM_ERR_FATAL );
      ret = false;
      goto TERMINATE;
    }
  }
  zMulMatTVecNC( qpasm->l, x, qpasm->z ); /* z = L^T x */
  /* the active constraints are added one by one, where linearly dependent ones are skipped. */
  for( i=0; i<zMatRowSizeNC(a); i++ )
    if( qpasm->is_active[i] ){
      qpasm->is_active[i] = false;
      _zQPASMAddCon( qpasm, i );
    }
 TERMINATE:
  zLPTableauDestroy( &tab );
  return ret;
}

/* solve the equality-constrained problem with the active constraints to find the temporary solution
 * zt and the Lagrange multipliers, namely,
 *  zt = U^T [ L_W^-1 b_W; -U2 cbar ], L_W^T lambda = L_W^-1 b_W + U1 cbar,
 * where U1 and U2 are the upper rows of U corresponding to the active constraints and the rest. */
static void _zQPASMSolveEq(zQPASM *qpasm, const zVec b, zVec zt)
{
  int i, j, n, k;
  double s, *t, *y, *lambda;

  n = zVecSizeNC(qpasm->z);
  k = zIndexSizeNC(qpasm->ia);
  t = zVecBufNC(qpasm->_t);
  y = zVecBufNC(qpasm->_y);
  lambda = zVecBufNC(qpasm->_lambda);
  for( i=0; i<n; i++ ) /* t = U cbar */
    t[i] = zRawVecInnerProd( zMatRowBufNC(qpasm->u,i), zVecBufNC(qpasm->cbar), n );
  for( i=0; i<k; i++ ){
    s = zVecElemNC(b,zIndexElemNC(qpasm->ia,i));
    for( j=0; j<i; j++ ) s -= zMatElemNC(qpasm->lw,i,j) * y[j];
    y[i] = s / zMatElemNC(qpasm->lw,i,i);
  }
  for( i=k-1; i>=0; i-- ){
    s = y[i] + t[i];
    for( j=i+1; j<k; j++ ) s -= zMatElemNC(qpasm->lw,j,i) * lambda[j];
    lambda[i] = s / zMatElemNC(qpasm->lw,i,i);
  }
  for( i=k; i<n; i++ ) y[i] = -t[i];
  zVecZero( zt );
  for( i=0; i<n; i++ )
    zRawVecCatDRC( zVecBufNC(zt), y[i], zMatRowBufNC(qpasm->u,i), n );
}

/* refine the solution by one step so as to satisfy the active constraints accurately in spite of
 * the round-off error in an ill-conditioned case, namely, zt += U1^T L_W^-1 ( b_W - abar_W zt ). */
static void _zQPASMRefine(zQPASM *qpasm, const zVec b, zVec zt)
{
  int i, j, n, k;
  double s, *r, *y;

  if( ( k = zIndexSizeNC(qpasm->ia) ) == 0 ) return;
  n = zVecSizeNC(zt);
  r = zVecBufNC(qpasm->_t);
  y = zVecBufNC(qpasm->_y);
  for( i=0; i<k; i++ ) /* residual of active constraints */
    r[i] = zVecElemNC(b,zIndexElemNC(qpasm->ia,i))
      - zRawVecInnerProd( zMatRowBufNC(qpasm->abar,zIndexElemNC(qpasm->ia,i)), zVecBufNC(zt), n );
  for( i=0; i<k; i++ ){
    s = r[i];
    for( j=0; j<i; j++ ) s -= zMatElemNC(qpasm->lw,i,j) * y[j];
    y[i] = s / zMatElemNC(qpasm->lw,i,i);
  }
  for( i=0; i<k; i++ )
    zRawVecCatDRC( zVecBufNC(zt), y[i], zMatRowBufNC(qpasm->u,i), n );
}

/* check if a blocking constraint i is more transversal to the step than another constraint j,
 * where api and apj are the inner products of the constraint vectors and the step. */
static bool _zQPASMIsMoreTransversal(const zMat a, int i, double api, int j, double apj)
{
  return api / zRawVecNorm( zMatRowBufNC(a,i), zMatColSizeNC(a) )
       < apj / zRawVecNorm( zMatRowBufNC(a,j), zMatColSizeNC(a) );
}

/* move the current feasible solution toward the temporary solution zt, and add a blocking
 * constraint to active set if any. If no constraint blocks, the current solution is
 * replaced with the temporary solution and the false value is returned.
 * Whether the current solution moves or not is stored where moved points. */
static bool _zQPASMAddEq(zQPASM *qpasm, const zVec b, const zVec zt, bool *moved)
{
  int i, imin;
  double ap, apmin = 0, az, d, dmin, scale;
  bool is_violated = false;

  *moved = false;
  /* the temporary solution coincides with the current solution if active set
   * determines the solution uniquely. */
  if( zIndexSizeNC(qpasm->ia) >= zVecSizeNC(zt) ) return false;
  /* a constraint is regarded as blocking if the inner product of the constraint vector and the step
   * is negative beyond the round-off errors, which are proportional to the norms of the current and
   * temporary solutions. A constraint linearly dependent on the active constraints is not blocking. */
  scale = zTOL * zMax( zVecNorm(zt), zVecNorm(qpasm->z) );
  for( dmin=HUGE_VAL, imin=-1, i=0; i<zMatRowSizeNC(qpasm->abar); i++ ){
    if( qpasm->is_active[i] ) continue;
    ap = zRawVecInnerProd( zMatRowBufNC(qpasm->abar,i), zVecBufNC(zt), zVecSizeNC(zt) );
    if( zVecElemNC(b,i) <= ap ) continue;
    ap -= ( az = zRawVecInnerProd( zMatRowBufNC(qpasm->abar,i), zVecBufNC(qpasm->z), zVecSizeNC(zt) ) );
    if( ap >= -scale * zRawVecNorm( zMatRowBufNC(qpasm->abar,i), zVecSizeNC(zt) ) ) continue;
    is_violated = true;
    if( ( d = ( zVecElemNC(b,i) - az ) / ap ) < 0 ) d = 0;
    if( d < dmin - zTOL ){
      dmin = d;
      imin = i;
      apmin = ap;
    } else
    if( d < dmin + zTOL ){ /* tie (degenerate case): prefer a constraint more transversal to the step */
      if( _zQPASMIsMoreTransversal( qpasm->abar, i, ap, imin, apmin ) ){
        imin = i;
        apmin = ap;
      }
      if( d < dmin ) dmin = d;
    }
  }
  if( !is_violated ){ /* full step */
    *moved = !zVecEqual( zt, qpasm->z, scale );
    zVecCopyNC( zt, qpasm->z );
    return false;
  }
  *moved = dmin > zTOL;
  zVecInterDivDRC( qpasm->z, zt, dmin );
  _zQPASMAddCon( qpasm, imin );
  return true;
}

/* delete an active constraint with the most negative Lagrange multiplier. */
static bool _zQPASMDelEq(zQPASM *qpasm)
{
  int i, imin;
  double lambda_min;

  for( lambda_min=0, imin=-1, i=0; i<zIndexSizeNC(qpasm->ia); i++ )
    if( zVecElemNC(qpasm->_lambda,i) < lambda_min ){
      lambda_min = zVecElemNC(qpasm->_lambda,i);
      imin = i;
    }
  if( imin < 0 ) return false;
  _zQPASMDelCon( qpasm, imin );
  return true;
}

/* check if vector/matrix sizes are consistent and the quadratic term matrix is symmetric. */
static bool _zQPASMCheck(const zMat q, const zVec c, const zMat a, const zVec b, zVec ans)
{
  if( !zMatIsSymmetric( q ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSYMMETRIC );
    return false;
  }
  if( !zMatRowVecSizeEqual( q, c ) || !zVecSizeEqual( c, ans ) ||
      !zMatRowVecSizeEqual( a, b ) || !zMatColVecSizeEqual( a, ans ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return false;
  }
  return true;
}

/* solve quadratic programming by active set method. */
bool zQPSolveASM(const zMat q, const zVec c, const zMat a, const zVec b, zVec ans, double *cost)
{
  zQPASM qpasm;
  zVec zt;
  int stall = 0;
  bool moved, ret = false;

  if( !_zQPASMCheck( q, c, a, b, ans ) ) return false;
  zt = zVecAlloc( zVecSizeNC(c) );
  if( !_zQPASMInit( &qpasm, q, c, a ) || !zt ) goto TERMINATE;
  if( !_zQPASMInitBase( &qpasm, a, b, ans ) ) goto TERMINATE;
  do{
    _zQPASMSolveEq( &qpasm, b, zt );
    if( !_zQPASMAddEq( &qpasm, b, zt, &moved ) && !_zQPASMDelEq( &qpasm ) ) break; /* optimal */
    /* the cost strictly decreases whenever the current solution moves, so that active set never
     * returns to the previous ones. the iteration is regarded as stalled if the current solution
     * does not move more times consecutively than the number of constraints. */
    stall = moved ? 0 : stall + 1;
  } while( stall <= zMatRowSizeNC(a) );
  if( stall > zMatRowSizeNC(a) ){
    ZRUNWARN( ZM_ERR_OPT_TRAPPED );
    goto TERMINATE;
  }
  _zQPASMRefine( &qpasm, b, zt );
  _zQPASMCholeskyBwdSubst( &qpasm, zt, ans, zVecBufNC(qpasm._t) );
  if( cost ) *cost = zQuadraticValue( q, c, ans );
  ret = true;
 TERMINATE:
  _zQPASMDestroy( &qpasm );
  zVecFree( zt );
  return ret;
}
