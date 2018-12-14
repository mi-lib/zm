/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_qp - optimization tools: quadratic programming.
 */

#include <zm/zm_opt.h>

/* zQuadraticValue
 * - calculation of a quadratic value.
 */
double zQuadraticValue(zMat q, zVec c, zVec x)
{
  register uint i, j;
  uint n;
  double val;

  n = zVecSize( x );
  if( !zMatIsSqr(q) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return 0;
  }
  if( zMatRowSize(q) != n ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return 0;
  }
  if( zVecSize(c) != n ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return 0;
  }
  for( val=0, i=0; i<n; i++ ){
    for( j=0; j<n; j++ )
      val += 0.5*zMatElem(q,i,j)*zVecElem(x,i)*zVecElem(x,j);
    val += zVecElem(c,i)*zVecElem(x,i);
  }
  return val;
}

/* zQPSolve
 * - quadratic programming solver.
 */
/* NOTE: this function is to be deleted due to mathematical illegality. */
bool zQPSolve(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost)
{
  register uint i, j;
  uint m, n;
  zMat d;
  zVec f, y;

  n = zVecSize( ans );
  m = b ?  zVecSize( b ) : 0;
  if( !zMatIsSqr(q) || zMatRowSize(q)!=n ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return false;
  }
  if( c && zVecSize(c)!=n ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return false;
  }
  if( a && zMatColSize(a) == 0 )
    a = NULL;
  if( a && ( zMatColSize(a)!=n || zMatRowSize(a)!=m ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return false;
  }

  d = zMatAllocSqr( n+m );
  f = zVecAlloc( n+m );
  y = zVecAlloc( n+m );
  for( i=0; i<n; i++ ){
    for( j=0; j<n; j++ )
      zMatSetElem( d, i, j,
        0.5*( zMatElem( q, i, j )+zMatElem( q, j, i ) ) );
    if( c ) zVecSetElem( f, i,-zVecElem( c, i ) );
  }
  for( i=0; i<m; i++ ){
    for( j=0; j<n; j++ ){
      zMatSetElem( d, n+i, j, zMatElem( a, i, j ) );
      zMatSetElem( d, j, n+i, zMatElem( a, i, j ) );
    }
    zVecSetElem( f, n+i, zVecElem( b, i ) );
  }
  zLESolveGauss( d, f, y );
  for( i=0; i<n; i++ )
    zVecSetElem( ans, i, zVecElem( y, i ) );

  zMatFree( d );
  zVecFreeAO( 2, f, y );
  if( cost ) *cost = zQuadraticValue( q, c, ans );
  return true;
}

/* (static)
 * _zQP2LCP
 * - transform quadratic programming problem to linear complementary problem.
 */
static bool _zQP2LCP(zMat q, zVec c, zMat a, zVec b, zMat *lm, zVec *lq, zVec *z);
bool _zQP2LCP(zMat q, zVec c, zMat a, zVec b, zMat *lm, zVec *lq, zVec *z)
{
  register int i, j;

  *lm = zMatAllocSqr( zMatRowSizeNC(q) + zMatRowSize(a) );
  *lq = zVecAlloc( zVecSizeNC(c) + zVecSize(b) );
  *z  = zVecAlloc( zVecSizeNC(c) + zVecSize(b) );
  if( !*lm || !*lq || !*z ){
    ZALLOCERROR();
    zMatFree( *lm );
    zVecFree( *lq );
    zVecFree( *z );
    return false;
  }
  zMatPut( *lm, 0, 0, q );
  zMatPut( *lm, zMatRowSizeNC(q), 0, a );
  for( i=0; i<zMatRowSize(a); i++ )
    for( j=0; j<zMatColSize(a); j++ )
      zMatSetElem( *lm, j, i+zMatColSizeNC(q), -zMatElem(a,i,j) );
  zVecPut( *lq, 0, c );
  for( i=0; i<zVecSize(b); i++ )
    zVecSetElem( *lq, zVecSizeNC(c)+i, -zVecElem(b,i) );
  return true;
}

/* zQPSolverDef
 * - define quadratic programming problem solver
 */
#define zQPSolverDef( name ) \
bool zQPSolve##name(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost)\
{\
  zMat lm;\
  zVec lq, z;\
  bool ret = false;\
\
  if( !_zQP2LCP( q, c, a, b, &lm, &lq, &z ) ) return false;\
\
  if( zLCPSolve##name( lm, lq, NULL, z ) ){\
    zVecGet( z, 0, ans );\
    if( cost ) *cost = zQuadraticValue( q, c, ans );\
    ret = true;\
  }\
  zMatFree( lm );\
  zVecFree( lq );\
  zVecFree( z );\
  return ret;\
}

/* zQPSolveLemke
 * - solve quadratic programming problem with Lemke's method.
 */
zQPSolverDef( Lemke )

/* zQPSolveIP
 * - solve quadratic programming problem with interior-point method.
 */
zQPSolverDef( IP )


/* zQPSolveASM
 * - solve a quadratic programming by active set method
 *   implemented by N. Wakisaka
 */
typedef struct{ /* list of active set indices */
  zIndex idx;
  double min;
} _zQPASMIndexData;
zListClass( _zQPASMIndexList, _zQPASMIndex, _zQPASMIndexData );

zVec _zQPSolveASMInitDefault(zMat a, zVec b, zVec ans, void *util)
{
  return zVecClear( ans );
}

double _zQPSolveASMConditionDefault(zMat a, zVec ans, int i, void *util)
{
  return zRawVecInnerProd(zMatRowBuf(a,i),zVecBuf(ans),zVecSizeNC(ans));
}

int _zQPSolveASMInitIndex(zIndex idx, zMat a, zVec b, zVec ans, void *util, double cond(zMat,zVec,int,void*))
{
  int m, i;

  /* initialize the active set of constraints */
  for( m=0, i=0; i<zVecSizeNC(b); i++ ){
    if( zIsTiny( cond( a, ans, i, util ) - zVecElem(b,i) ) ){
      zIndexElem(idx,i) = 1;
      m++;
    } else
      zIndexElem(idx,i) = 0;
  }
  return m;
}

#define ZM_OPT_QP_ASM_TOL 1.0e-8
bool _zQPSolveASM(zMat q, zVec c, zMat a, zVec b, zVec ans, zIndex idx, void *util, double cond(zMat,zVec,int,void*))
{
  zMat qa;
  zVec xy, cb, d;
  double tempd, tempd2, objv;
  _zQPASMIndexList ilist;
  _zQPASMIndex *idata;
  bool endflag, ret = true;
  int tempi;
  int n, nm, m;
  register int i, j, k;

  n = zVecSizeNC(ans);
  zListInit( &ilist );
  d = zVecAlloc(n);
  qa = zMatAllocSqr( n + zVecSizeNC(b) );
  xy = zVecAlloc( n + zVecSizeNC(b) );
  cb = zVecAlloc( n + zVecSizeNC(b) );
  if( !d || !qa || !xy || !cb ){
    ret = false;
    goto RET;
  }
  if( cond == NULL ) cond = _zQPSolveASMConditionDefault;
  m = _zQPSolveASMInitIndex( idx, a, b, ans, util, cond );
  do{
    nm = n + m;
    zMatSetSize( qa, nm, nm );
    zVecSetSize( xy, nm );
    zVecSetSize( cb, nm );

    for( i=0; i<n; i++ )
      for( j=0; j<n; j++ )
        zMatElem(qa,i,j) = -zMatElem(q,i,j);
    for( k=0, j=n; j<nm; j++ ){
      while( !zIndexElem(idx,k) && k < zVecSizeNC(b) ) k++;
      for( i=0; i<n; i++ )
        zMatElem(qa,i,j) = zMatElem(a,k,i);
      k++;
    }
    for( k=0, i=n; i<nm; i++ ){
      while( !zIndexElem(idx,k) && k < zVecSizeNC(b) ) k++;
      for( j=0; j<n; j++ )
        zMatElem(qa,i,j) = zMatElem(a,k,j);
      k++;
    }
    for( i=n; i<nm; i++ )
      for( j=n; j<nm; j++ )
        zMatElem(qa,i,j) = 0.0;
    for( i=0; i<n; i++ )
      zVecElem(cb,i) = zVecElem(c,i);
    for( k=0; i<nm; i++ ){
      while( !zIndexElem(idx,k) && k < zVecSizeNC(b) ) k++;
      zVecElem(cb,i) = zVecElem(b,k);
      k++;
    }
    zLESolveMP( qa, cb, NULL, NULL, xy );

    for( i=0; i<n; i++ )
      if( !zIsTiny( zVecElem(xy,i) - zVecElem(ans,i) ) ) goto STEP2;

    for( i=0; i<n; i++ )
      zVecElem(ans,i) = zVecElem(xy,i);
    for( i=0; i<m; i++ )
      if( zVecElem(xy,n+i) < 0 ) goto NEXT;
    goto RET; /* found the optimal solution */

   NEXT:
    tempd = zVecElem(xy,n);
    for( i=1; i<m; i++ )
      if( zVecElem(xy,n+i) < tempd )
        tempd = zVecElem(xy,n+i);
    tempi = 0;
    for( i=0; i<zVecSizeNC(b); i++ )
      if( zIndexElem(idx,i) != 0 ){
        if( fabs( zVecElem(xy,tempi+n) - tempd ) < ZM_OPT_QP_ASM_TOL ){
          zIndexElem(idx,i) = 0;
          m--;
        }
        tempi++;
      }
    continue;

   STEP2:
    /* find a new feasible direction */
    for( i=0; i<n; i++ )
      zVecElem(d,i) = zVecElem(xy,i) - zVecElem(ans,i);
    tempd = 1.0;
    for( i=0; i<zVecSizeNC(b); i++ ){
      tempd2 = zRawVecInnerProd( zMatRowBuf(a,i), zVecBuf(d), n );
      if( zIndexElem(idx,i) == 0 && tempd2 < 0 ){
        tempd2 = ( zVecElem(b,i) - cond( a, ans, i, util ) ) / tempd2;
        if( tempd2 < tempd ) tempd = tempd2;
      }
    }
    for( i=0; i<n; i++ )
      zVecElem(ans,i) += tempd * zVecElem(d,i);
    for( i=0; i<zVecSizeNC(b); i++ )
      if( zIndexElem(idx,i) == 0 && zIsTiny( cond( a, ans, i, util ) - zVecElem(b,i) ) ){
        zIndexElem(idx,i) = 1;
        m++;
      }
    /* check if circulation happens due to degeneracy */
    endflag = false;
    objv = zQuadraticValue( q, c, ans );
    zListForEach( &ilist, idata ){
      for( i=0; i<zVecSizeNC(b);i++ )
        if( zIndexElem(idx,i) != zIndexElem(idata->data.idx,i) ) goto CONTINUE;
      if( fabs( idata->data.min - objv ) > ZM_OPT_QP_ASM_TOL ) goto CONTINUE;
      endflag = true;
      break;
     CONTINUE: ;
    }
    if( endflag == true ) goto RET;

    /* register to the list of bases */
    idata = zAlloc( _zQPASMIndex, 1 );
    idata->data.idx = zIndexCreate( zVecSizeNC(b) );
    for( i=0; i<zVecSizeNC(b); i++ )
      zIndexElem(idata->data.idx,i) = zIndexElem(idx,i);
    idata->data.min = objv;
    zListInsertTail( &ilist ,idata );
  } while( 1 );

 RET:
  zListForEach( &ilist, idata )
    zIndexFree( idata->data.idx );
  zListDestroy( _zQPASMIndex, &ilist );
  zMatFree( qa );
  zVecFreeAO( 3, xy, cb, d );
  return ret;
}

bool zQPSolveASM(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost, zVec init(zMat,zVec,zVec,void*), void *util)
{
  zIndex idx;
  bool ret;

  if( !init ) init = _zQPSolveASMInitDefault;
  init( a, b, ans, util );
  idx = zIndexCreate( zVecSizeNC(b) );
  ret = _zQPSolveASM( q, c, a, b, ans, idx, util, _zQPSolveASMConditionDefault );
  zIndexFree( idx );
  return ret;
}

/* zCGSolve
 * - quadratic programming solver by conjugate gradient method.
 */
double zCGSolve(zMat q, zVec c, zVec ans, int iter)
{
  zVec d, g, qd;
  double a, b, s, result;
  register int i=0;

  d = zVecAlloc( zVecSizeNC(ans) );
  g = zVecAlloc( zVecSizeNC(ans) );
  qd= zVecAlloc( zVecSizeNC(ans) );
  if( !d || !g || !qd ) goto TERMINATE;

  zMulMatVec( q, ans, g );
  zVecAddDRC( g, c );
  if( zVecIsTiny( g ) ) goto TERMINATE;
  zVecCopy( g, d );
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zMulMatVecNC( q, d, qd );
    s = zVecInnerProd( d, qd );
    a = zVecInnerProd( d, g ) / s;
    zVecCatNCDRC( ans, -a, d );
    zVecCatNCDRC( g, -a, qd );
    if( zVecIsTiny( g ) ) goto TERMINATE;
    b = zVecInnerProd( g, qd ) / s;
    zVecCat( g, b, d, d );
  }
  ZITERWARN( iter );

 TERMINATE:
  result = zQuadraticValue( q, c, ans );
  zVecFreeAO( 3, d, g, qd );
  return result;
}
