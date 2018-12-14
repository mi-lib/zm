/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_pex_eq - polynomial expression class: polynomial equation solver.
 */

#include <zm/zm_pex.h>

/* Bairstow-Hitchcock's method */
static bool _zPexBH(double *p,double *q,zPex a, int n, double tol, int iter);

/* (static)
 * _zPexBH
 * - iterative calculation of quadratic factor.
 */
bool _zPexBH(double *p,double *q,zPex a, int n, double tol, int iter)
{
  double *b, *c;
  double dp, dq, d;
  int i, j;
  bool result = true;

  b = zAlloc( double, n + 1 );
  c = zAlloc( double, n );
  if( !b || !c ){
    ZALLOCERROR();
    result = false;
    goto TERMINATE;
  }
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    /* coefficients of quontient */
    b[n]   = zPexCoeff(a,n);
    b[n-1] = zPexCoeff(a,n-1) - *p * b[n];
    for( j=2; j<=n; j++ )
      b[n-j] = zPexCoeff(a,n-j) - *p * b[n-j+1] - *q * b[n-j+2];
    /* gradient of divider coefficients */
    c[n-1] = b[n];
    c[n-2] = b[n-1] - *p * c[n-1];
    for( j=n-3; j>=0; j-- )
      c[j] = b[j+1] - *p * c[j+1] - *q * c[j+2];
    /* delta of divider coefficients
       (based on Newton=Raphson's method) */
    d  = c[1] * c[1] +( b[1] - c[0] ) * c[2];
    dp = ( b[1] * c[1] - b[0] * c[2] ) / d;
    dq = ( b[0] * c[1] + b[1] * ( b[1] - c[0] ) ) / d;
    /* termination judgement */
    if( zIsTol(dp,tol) && zIsTol(dq,tol) ) goto SUCCESS;
    *p += dp;
    *q += dq;
  }
  ZITERWARN( iter );

 SUCCESS: /* shift of coefficients */
  for( i=2; i<=n; i++ )
    zPexSetCoeff( a, i-2, b[i] );
 TERMINATE:
  zFree( b );
  zFree( c );
  return result;
}

/* zPexBHDST
 * - numerical solution of polynomial equation
 *   based on Bairstow-Hitchcock's method.
 *   (destructive)
 */
zComplex *zPexBHDST(zPex a, zComplex *ans, double tol, int iter)
{
  int n;
  double p, q;

  if( tol == 0 ) tol = ZM_PEX_EQ_TOL;
  for( n=zPexDim(a); n>0; n-=2 )
    switch( n ){
    case 1: /* linear equation */
      zComplexCreate( &ans[0], -zPexCoeff(a,0)/zPexCoeff(a,1), 0 );
      return ans;
    case 2: /* quadratic equation */
      zQESolve( zPexCoeff(a,2), zPexCoeff(a,1), zPexCoeff(a,0), &ans[0] );
      return ans;
    default:
      /* initial values of p and q for iteration */
      p = zPexCoeff(a,1) + 1.0;
      q = zPexCoeff(a,0) + 1.0;
      /* division by quadratic expression
         based on Bairstow's method */
      if( !_zPexBH( &p, &q, a, n, tol, iter ) ){
        ZALLOCERROR();
        return NULL;
      }
      /* append solution of 'x^2 + p x + q = 0' */
      zQESolve( 1, p, q, &ans[n-2] );
    }
  return ans;
}

/* zPexBH
 * - numerical solution of polynomial equation
 *   based on Bairstow-Hitchcock's method.
 */
zComplex *zPexBH(zPex a, zComplex *ans, double tol, int iter)
{
  zPex ac;

  if( !( ac = zPexClone( a ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  ans = zPexBHDST( ac, ans, tol, iter );
  zPexFree( ac );
  return ans;
}

/* Durand-Kerner-Aberth's method */

/* zPexDKA
 * - numerical solution of polynomial equation
 *   based on Durand-Kerner-Aberth's method.
 */
zComplex *zPexDKA(zPex a, zComplex *ans, double tol, int iter)
{
  register int i, j, k, n;
  int count;
  double s, r, r0;
  zComplex *p, tmp1, tmp2;
  zPex b;

  n = zPexDim(a);
  p = zAlloc( zComplex, n );
  b = zPexAlloc( n );
  if( !p || !b ){
    ZALLOCERROR();
    return NULL;
  }
  if( tol == 0 ) tol = ZM_PEX_EQ_TOL;
  /* Aberth-Iri's initialization */
  s = zPexCoeff(a,n-1) / ( zPexCoeff(a,n) * n );
  zPexModulo( a, s, b );
  for( r0=0, i=1; i<=n; i++ ){
    r = pow( n*fabs(zPexCoeffHigh(b,i)/zPexCoeff(b,n)), 1.0/(double)i );
    if( r0 < r ) r0 = r;
  }
  for( i=0; i<n; i++ ){
    zComplexPolar( &ans[i], r0, ( zPIx2*i + 1.5 )/n );
    ans[i].re -= s;
  }
  /* Durand's method with Kerner's interpretation */
  ZITERINIT( iter );
  for( k=0; k<iter; k++ ){
    for( i=0; i<n; i++ ){
      zComplexCreate( &p[i], zPexCoeff(a,n), 0.0 );
      for( j=0; j<n; j++ ){
        if( j == i ) continue;
        zComplexSub( &ans[i], &ans[j], &tmp1 );
        zComplexCMul( &p[i], &tmp1, &tmp2 );
        zComplexCopy( &tmp2, &p[i] );
      }
    }
    for( count=0, i=0; i<n; i++ ){
      zPexCVal( a, &ans[i], &tmp1 );
      if( zComplexIsTol( &tmp1, tol ) ) count++;
      zComplexCDiv( &tmp1, &p[i], &tmp2 );
      zComplexSub( &ans[i], &tmp2, &ans[i] );
    }
    if( count == n ) goto TERMINATE;
  }
  ZITERWARN( iter );
 TERMINATE:
  zFree( p );
  zPexFree( b );
  return ans;
}
