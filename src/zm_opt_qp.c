/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_qp - optimization tools: quadratic programming.
 */

#include <zm/zm_opt.h>

/* calculation of a quadratic value. */
double zQuadraticValue(const zMat q, const zVec c, const zVec x)
{
  int i, j, n;
  double val;

  n = zVecSize(x);
  if( !zMatIsSqr( q ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return 0;
  }
  if( zMatRowSize(q) != n ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return 0;
  }
  if( zVecSize(c) != n ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return 0;
  }
  for( val=0, i=0; i<n; i++ ){
    for( j=0; j<n; j++ )
      val += 0.5*zMatElemNC(q,i,j)*zVecElemNC(x,i)*zVecElemNC(x,j);
    val += zVecElemNC(c,i)*zVecElemNC(x,i);
  }
  return val;
}

#if 0
/* quadratic programming solver. */
/* NOTE: this function is to be deleted due to mathematical illegality. */
bool zQPSolve(const zMat q, const zVec c, const zMat a, const zVec b, zVec ans, double *cost)
{
  int i, j, m, n;
  zMat d;
  zVec f, y;

  n = zVecSize( ans );
  m = b ?  zVecSize( b ) : 0;
  if( !zMatIsSqr(q) || zMatRowSize(q)!=n ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return false;
  }
  if( c && zVecSize(c)!=n ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return false;
  }
  if( a && zMatColSize(a) == 0 )
    a = NULL;
  if( a && ( zMatColSize(a)!=n || zMatRowSize(a)!=m ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return false;
  }

  d = zMatAllocSqr( n+m );
  f = zVecAlloc( n+m );
  y = zVecAlloc( n+m );
  for( i=0; i<n; i++ ){
    for( j=0; j<n; j++ )
      zMatSetElemNC( d, i, j,
        0.5*( zMatElemNC( q, i, j )+zMatElemNC( q, j, i ) ) );
    if( c ) zVecSetElemNC( f, i,-zVecElemNC( c, i ) );
  }
  for( i=0; i<m; i++ ){
    for( j=0; j<n; j++ ){
      zMatSetElemNC( d, n+i, j, zMatElemNC( a, i, j ) );
      zMatSetElemNC( d, j, n+i, zMatElemNC( a, i, j ) );
    }
    zVecSetElemNC( f, n+i, zVecElemNC( b, i ) );
  }
  zLESolveGauss( d, f, y );
  for( i=0; i<n; i++ )
    zVecSetElemNC( ans, i, zVecElemNC( y, i ) );

  zMatFree( d );
  zVecFreeAtOnce( 2, f, y );
  if( cost ) *cost = zQuadraticValue( q, c, ans );
  return true;
}
#endif

/* transform a quadratic programming problem to a linear complementary problem. */
static bool _zQP2LCP(const zMat q, const zVec c, const zMat a, const zVec b, zMat *lm, zVec *lq, zVec *z)
{
  int i, j;

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
      zMatSetElemNC( *lm, j, i+zMatColSizeNC(q), -zMatElemNC(a,i,j) );
  zVecPut( *lq, 0, c );
  for( i=0; i<zVecSize(b); i++ )
    zVecSetElemNC( *lq, zVecSizeNC(c)+i, -zVecElemNC(b,i) );
  return true;
}

/* define quadratic programming problem solver */
#define zQPSolverDef( name ) \
bool zQPSolve##name(const zMat q, const zVec c, const zMat a, const zVec b, zVec ans, double *cost)\
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

/* solve quadratic programming problem with Lemke's method. */
zQPSolverDef( Lemke )

/* solve quadratic programming problem with interior-point method. */
zQPSolverDef( IP )

/* quadratic programming solver by conjugate gradient method. */
double zCGSolve(const zMat q, const zVec c, zVec ans, int iter)
{
  zVec d, g, qd;
  double a, b, s, result;
  int i = 0;

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
  zVecFreeAtOnce( 3, d, g, qd );
  return result;
}
