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

/* zPexIPAlloc
 * - allocate polynomial curve.
 */
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

/* zPexIPAllocBoundary
 * - allocate a polynomial curve from boundary condition.
 */
bool zPexIPAllocBoundary(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, zVec v)
{
  double k1, k2, k3; /* internal working variables */
  register uint i, n1, n2, n3;
  uint n;

  if( !zPexIPAlloc( pc, term, zVecSize(v) + 5 ) )
    return false;
  n = zPexIPDim( pc );
  n1 = n - 1;
  n2 = n - 2;
  n3 = n - 3;

  k1 = x2 - x1 - ( v1 + 0.5*a1*term ) * term;
  k2 = ( v2 - v1 - a1*term ) * term;
  k3 = ( a2 - a1 )*zSqr(term);

  zPexIPSetCoeff( pc, 0, x1 );
  zPexIPSetCoeff( pc, 1, v1*term );
  zPexIPSetCoeff( pc, 2, 0.5*a1*zSqr(term) );
  zPexIPSetCoeff( pc,n2, 0.5*k1*n*n1  -k2*n1     +0.5*k3 );
  zPexIPSetCoeff( pc,n1,    -k1*n*n2  +k2*(2*n-3)-    k3 );
  zPexIPSetCoeff( pc, n, 0.5*k1*n1*n2 -k2*n2     +0.5*k3 );
  for( i=3; i<=n3; i++ ){
    zPexIPSetCoeff( pc, i, ( k1 = zVecElem( v, i-3 ) ) );
    zVecElem(pc->c,n ) += -0.5*(n1-i)*(n2-i)*k1;
    zVecElem(pc->c,n1) +=      (n-i) *(n2-i)*k1;
    zVecElem(pc->c,n2) += -0.5*(n-i) *(n1-i)*k1;
  }
  return true;
}

/* zPexIPCreateLSM
 * - creation of polynomial curve based on the least square method.
 */
bool zPexIPCreateLSM(zPexIP *pc, double term, int dim, zVec t, zVec x)
{
  register uint i, j, k, n, m;
  zMat a;
  zVec b, v;
  bool result = true;

  if( !zVecSizeIsEqual( t, x ) ){
    ZRUNERROR( ZM_ERR_IP_SIZMIS );
    return false;
  }

  n = zVecSizeNC( t );
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
      zVecSetElem( v, j, pow( zVecElem(t,i)/term, j ) );
    for( j=0; j<m; j++ ){
      for( k=0; k<m; k++ )
        zMatElem(a,j,k) += zVecElem(v,j+k);
      zVecElem(b,j) += zVecElem(x,i) * zVecElem(v,j);
    }
  }
  zLESolveGauss( a, b, pc->c );

 TERMINATE:
  zMatFree( a );
  zVecFreeAO( 2, b, v );
  return result;
}

/* zPexIPFree
 * - destruction of polynomial curve.
 */
void zPexIPFree(zPexIP *pc)
{
  zPexFree( pc->c );
  zPexIPSetTerm( pc, 1 ); /* dummy */
}

/* zPexIPVal
 * - value of polynomial curve.
 */
double zPexIPVal(zPexIP *pc, double t)
{
  register int i, n;
  double value = 0, term;

  n = zPexIPDim( pc );
  if( ( t /= zPexIPTerm( pc ) ) == 0 )
    value = zPexIPCoeff( pc, 0 );
  else if( fabs(t) < 1.0 )
    for( i=n, term=pow(t,n); i>=0; i--, term/=t )
      value += zPexIPCoeff( pc, i ) * term;
  else
    for( i=0, term=1.0; i<=n; i++, term*=t )
      value += zPexIPCoeff( pc, i ) * term;
  return value;
}

/* zPexIPFWrite
 * - expression of polynomial curve.
 */
void zPexIPFWrite(FILE *fp, zPexIP *pc)
{
  register int i;

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
