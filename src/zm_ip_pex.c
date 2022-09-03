/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_pex - interpolation: polynomial curve.
 */

#include <zm/zm_ip.h>

/* ********************************************************** */
/* CLASS: zPexIP
 * polynomial curve class with designable coefficients.
 * ********************************************************** */

/* allocate a polynomial curve. */
bool zPexIPAlloc(zPexIP *pc, double term, int dim)
{
  if( term == 0 ){
    ZRUNERROR( ZM_ERR_IP_INVTERM );
    return false;
  }
  if( !( pc->c = zPexAlloc( dim ) ) ) return false;
  zPexIPSetTerm( pc, term );
  return true;
}

/* free a polynomial curve. */
void zPexIPFree(zPexIP *pc)
{
  zPexFree( pc->c );
  zPexIPSetTerm( pc, 1 ); /* dummy */
}

/* set the boundary condition of a polynomial curve. */
bool zPexIPBoundary(zPexIP *pc, double x1, double v1, double a1, double x2, double v2, double a2, zVec v)
{
  double k1, k2, k3; /* internal working variables */
  uint i, n1, n2, n3, n;

  if( ( n = zPexIPDim(pc) ) != zVecSize(v) + 5 ) return false;
  n1 = n - 1;
  n2 = n - 2;
  n3 = n - 3;

  k1 = x2 - x1 - ( v1 + 0.5*a1*zPexIPTerm(pc) ) * zPexIPTerm(pc);
  k2 = ( v2 - v1 - a1*zPexIPTerm(pc) ) * zPexIPTerm(pc);
  k3 = ( a2 - a1 )*zSqr(zPexIPTerm(pc));

  zPexIPSetCoeff( pc, 0, x1 );
  zPexIPSetCoeff( pc, 1, v1*zPexIPTerm(pc) );
  zPexIPSetCoeff( pc, 2, 0.5*a1*zSqr(zPexIPTerm(pc)) );
  zPexIPSetCoeff( pc,n2, 0.5*k1*n*n1  -k2*n1     +0.5*k3 );
  zPexIPSetCoeff( pc,n1,    -k1*n*n2  +k2*(2*n-3)-    k3 );
  zPexIPSetCoeff( pc, n, 0.5*k1*n1*n2 -k2*n2     +0.5*k3 );
  for( i=3; i<=n3; i++ ){
    zPexIPSetCoeff( pc, i, ( k1 = zVecElemNC( v, i-3 ) ) );
    zVecElemNC(pc->c,n ) += -0.5*(n1-i)*(n2-i)*k1;
    zVecElemNC(pc->c,n1) +=      (n-i) *(n2-i)*k1;
    zVecElemNC(pc->c,n2) += -0.5*(n-i) *(n1-i)*k1;
  }
  return true;
}

/* allocate a polynomial curve from the boundary condition. */
bool zPexIPCreateBoundary(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, zVec v)
{
  if( !zPexIPAlloc( pc, term, zVecSize(v) + 5 ) ) return false;
  return zPexIPBoundary( pc, x1, v1, a1, x2, v2, a2, v );
}

/* fit a polynomial curve to a sequence of points based on the least square method. */
bool zPexIPLSM(zPexIP *pc, zVec t, zVec x)
{
  uint i, j, k, n, m;
  zMat a;
  zVec b, v;
  bool result = true;

  if( !zVecSizeIsEqual( t, x ) ){
    ZRUNERROR( ZM_ERR_IP_SIZMIS );
    return false;
  }
  n = zVecSizeNC( t );
  m = zPexIPDim(pc) + 1;
  a = zMatAllocSqr( m );
  b = zVecAlloc( m );
  v = zVecAlloc( 2 * m );
  if( !a || !b || !v ){
    result = false;
    goto TERMINATE;
  }
  for( i=0; i<n; i++ ){
    for( j=0; j<zVecSizeNC(v); j++ )
      zVecSetElemNC( v, j, pow( zVecElemNC(t,i)/zPexIPTerm(pc), j ) );
    for( j=0; j<m; j++ ){
      for( k=0; k<m; k++ )
        zMatElemNC(a,j,k) += zVecElemNC(v,j+k);
      zVecElemNC(b,j) += zVecElemNC(x,i) * zVecElemNC(v,j);
    }
  }
  zLESolveGauss( a, b, pc->c );
 TERMINATE:
  zMatFree( a );
  zVecFreeAO( 2, b, v );
  return result;
}

/* create a polynomial curve that fits a sequence of points based on the least square method. */
bool zPexIPCreateLSM(zPexIP *pc, double term, int dim, zVec t, zVec x)
{
  if( !zPexIPAlloc( pc, term, dim ) ) return false;
  return zPexIPLSM( pc, t, x );
}

/* create a polynomial curve from the boundary condition and a sequence of points to fit. */
bool zPexIPCreateBounderyLSM(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, int dim, zVec t, zVec x)
{
  uint i, j, k, n, m;
  zMat a;
  zVec b, v;
  bool result = true;

  if( !zVecSizeIsEqual( t, x ) ){
    ZRUNERROR( ZM_ERR_IP_SIZMIS );
    return false;
  }

  n = zVecSizeNC(t);
  if( !zPexIPAlloc( pc, term, dim ) ) return false;
  m = dim + 1;
  a = zMatAllocSqr( m );
  b = zVecAlloc( m );
  v = zVecAlloc( 2 * m );
  if( !a || !b || !v ){
    result = false;
    goto TERMINATE;
  }

  for( i=0; i<n; i++ ){
    for( j=0; j<zVecSizeNC(v); j++ )
      zVecSetElemNC( v, j, pow( zVecElemNC(t,i)/term, j ) );
    for( j=3; j<m-3; j++ ){
      for( k=0; k<m; k++ )
        zMatElemNC(a,j,k) += zVecElemNC(v,j+k);
      zVecElemNC(b,j) += zVecElemNC(x,i) * zVecElemNC(v,j);
    }
  }
  zMatSetElemNC( a, 0, 0, 1 );
  zMatSetElemNC( a, 1, 1, 1/term );
  zMatSetElemNC( a, 2, 2, 2/zSqr(term) );
  for( i=0; i<m; i++ )
    zMatSetElemNC( a, m-3, i, 1 );
  for( i=1; i<m; i++ )
    zMatSetElemNC( a, m-2, i, i );
  for( i=2; i<m; i++ )
    zMatSetElem( a, m-1, i, i*(i-1) );
  zVecSetElemNC( b, 0, x1 );
  zVecSetElemNC( b, 1, v1 );
  zVecSetElemNC( b, 2, a1 );
  zVecSetElemNC( b, m-3, x2 );
  zVecSetElemNC( b, m-2, v2 );
  zVecSetElemNC( b, m-1, a2 );

  zLESolveGauss( a, b, pc->c );

 TERMINATE:
  zMatFree( a );
  zVecFree( b );
  zVecFree( v );
  return result;
}

/* value of a polynomial curve. */
double zPexIPVal(zPexIP *pc, double t)
{
  return zPexVal( pc->c, t/zPexIPTerm(pc) );
}

/* velocity of a polynomial curve. */
double zPexIPVel(zPexIP *pc, double t)
{
  return zPexDifVal( pc->c, 1, t/zPexIPTerm(pc) ) / zPexIPTerm(pc);
}

/* acceleration of a polynomial curve. */
double zPexIPAcc(zPexIP *pc, double t)
{
  return zPexDifVal( pc->c, 2, t/zPexIPTerm(pc) ) / zSqr( zPexIPTerm(pc) );
}

/* print expression of a polynomial curve. */
void zPexIPFPrint(FILE *fp, zPexIP *pc)
{
  int i;

  for( i=zPexIPDim(pc); i>=0; i-- ){
    if( i < zPexIPDim(pc) && zPexIPCoeff(pc,i) > 0 )
      fprintf( fp, " + " );
    if( zPexIPCoeff(pc,i) != 0 ){
      fprintf( fp, "%.10g ", zPexIPCoeff(pc,i) );
      if( i > 0 ) fprintf( fp, "(t/%.10g)", zPexIPTerm(pc) );
      if( i > 1 ) fprintf( fp, "^%d", i );
    }
  }
  fprintf( fp, "\n" );
}
