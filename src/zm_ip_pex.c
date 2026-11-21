/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_pex - interpolation: polynomial function.
 */

#include <zm/zm_ip.h>

/* allocate a polynomial function. */
bool zPexIPAlloc(zPexIP *pc, double term, int dim)
{
  if( term == 0 ){
    ZRUNERROR( ZM_ERR_IP_INVALID_TERM );
    return false;
  }
  if( !( pc->c = zPexAlloc( dim ) ) ) return false;
  zPexIPSetTerm( pc, term );
  return true;
}

/* free a polynomial function. */
void zPexIPFree(zPexIP *pc)
{
  zPexFree( pc->c );
  zPexIPSetTerm( pc, 1 ); /* dummy */
}

/* allocate a polynomial function from the boundary condition. */
bool zPexIPCreateBoundary(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, const zVec v)
{
  double k1, k2, k3; /* internal working variables */
  int i, n1, n2, n3, n;

  if( !zPexIPAlloc( pc, term, ( n = zVecSize(v) + 5 ) ) ) return false;
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

/* create a polynomial function that fits a sequence of points based on the least square method,
 * where the coefficients are found by zLESolveErrorMin() without forming the normal equation. */
bool zPexIPCreateLSM(zPexIP *pc, double term, int dim, const zVec t, const zVec x)
{
  int i, j;
  double s, p;
  zMat a;
  zVec c;
  bool result = false;

  if( !zVecSizeEqual( t, x ) ){
    ZRUNERROR( ZM_ERR_IP_SIZEMISMATCH );
    return false;
  }
  if( term == 0 ){
    ZRUNERROR( ZM_ERR_IP_INVALID_TERM );
    return false;
  }
  if( dim < 0 ){
    ZRUNERROR( ZM_ERR_PEX_INVALID_DIM );
    return false;
  }
  a = zMatAlloc( zVecSizeNC(t), dim + 1 );
  c = zVecAlloc( dim + 1 );
  if( !a || !c ) goto TERMINATE;
  for( i=0; i<zVecSizeNC(t); i++ ){ /* Vandermonde matrix */
    s = zVecElemNC(t,i) / term;
    for( p=1, j=0; j<=dim; j++, p*=s ) zMatSetElemNC( a, i, j, p );
  }
  if( !zLESolveErrorMin( a, x, NULL, c ) || !zPexIPAlloc( pc, term, dim ) ) goto TERMINATE;
  zVecCopyNC( c, pc->c );
  result = true;
 TERMINATE:
  zMatFree( a );
  zVecFree( c );
  return result;
}

/* create a polynomial function from the boundary condition and a sequence of points to fit.
 * A polynomial function that satisfies the boundary condition is expressed as
 *   f(t) = f0(t) + sum_k v_k phi_k(t),
 * where f0 is created by zPexIPCreateBoundary() for v = 0, and phi_k is created by
 * zPexIPCreateBoundary() for the null boundary condition and v = e_k (the k-th unit vector).
 * The coefficients v_k are found by the least square method without any constraints. */
bool zPexIPCreateBoundaryLSM(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, int dim, const zVec t, const zVec x)
{
  zPexIP phi;
  int i, k, n;
  zMat a = NULL;
  zVec r = NULL, v = NULL;
  bool result = false;

  if( !zVecSizeEqual( t, x ) ){
    ZRUNERROR( ZM_ERR_IP_SIZEMISMATCH );
    return false;
  }
  if( dim < 5 ){
    ZRUNERROR( ZM_ERR_IP_TOOSMALLDIM, dim, 5 );
    return false;
  }
  if( dim == 5 ) /* no free coefficient */
    return zPexIPCreateBoundary( pc, term, x1, v1, a1, x2, v2, a2, NULL );
  n = zVecSizeNC(t);
  a = zMatAlloc( n, dim - 5 );
  r = zVecAlloc( n );
  v = zVecAlloc( dim - 5 );
  if( !a || !r || !v ) goto TERMINATE;
  for( k=0; k<zVecSizeNC(v); k++ ){ /* basis functions */
    zVecZero( v );
    zVecSetElemNC( v, k, 1 );
    if( !zPexIPCreateBoundary( &phi, term, 0, 0, 0, 0, 0, 0, v ) ) goto TERMINATE;
    for( i=0; i<n; i++ )
      zMatSetElemNC( a, i, k, zPexIPVal( &phi, zVecElemNC(t,i) ) );
    zPexIPFree( &phi );
  }
  zVecZero( v );
  if( !zPexIPCreateBoundary( &phi, term, x1, v1, a1, x2, v2, a2, v ) ) goto TERMINATE;
  for( i=0; i<n; i++ ) /* residual from f0 */
    zVecSetElemNC( r, i, zVecElemNC(x,i) - zPexIPVal( &phi, zVecElemNC(t,i) ) );
  zPexIPFree( &phi );
  if( !zLESolveErrorMin( a, r, NULL, v ) ) goto TERMINATE;
  result = zPexIPCreateBoundary( pc, term, x1, v1, a1, x2, v2, a2, v );
 TERMINATE:
  zMatFree( a );
  zVecFree( r );
  zVecFree( v );
  return result;
}

/* value of a polynomial function. */
double zPexIPVal(const zPexIP *pc, double t)
{
  return zPexVal( pc->c, t/zPexIPTerm(pc) );
}

/* velocity of a polynomial function. */
double zPexIPVel(const zPexIP *pc, double t)
{
  return zPexDifVal( pc->c, 1, t/zPexIPTerm(pc) ) / zPexIPTerm(pc);
}

/* acceleration of a polynomial function. */
double zPexIPAcc(const zPexIP *pc, double t)
{
  return zPexDifVal( pc->c, 2, t/zPexIPTerm(pc) ) / zSqr( zPexIPTerm(pc) );
}

/* print expression of a polynomial function. */
void zPexIPFPrint(FILE *fp, const zPexIP *pc)
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
