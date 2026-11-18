/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_pex_eq - polynomial expression class: polynomial equation solver.
 */

#include <zm/zm_pex.h>

/* common procedures */

/* the exponent e such that 2^e is close to an upper bound of the absolute values of the roots of a
 * polynomial, max_k |a_{n-k}/a_n|^(1/k). */
static int _zPexEqScaleExp(const zPex pex)
{
  int k, n, e;
  double r, rmax = 0;

  n = zPexDim(pex);
  for( k=1; k<=n; k++ )
    if( ( r = pow( fabs( zPexCoeff(pex,n-k) / zPexCoeff(pex,n) ), 1.0/k ) ) > rmax ) rmax = r;
  if( !( rmax > 0 ) || !zIsFinite( rmax ) ) return 0;
  frexp( rmax, &e );
  return e;
}

/* scale a polynomial pex(x) to b(y) = pex(2^e y) / 2^f, where 2^f is close to the leading coefficient of
 * pex(2^e y). Since the coefficients are multiplied by powers of two, no rounding error is introduced. */
static zPex _zPexEqScale(const zPex pex, int e, zPex b)
{
  int k, n, f;

  n = zPexDim(pex);
  frexp( zPexCoeff(pex,n), &f );
  for( k=0; k<=n; k++ )
    zPexSetCoeff( b, k, ldexp( zPexCoeff(pex,k), ( k - n ) * e - f ) );
  return b;
}

/* scale the roots of a scaled polynomial back by 2^e. */
static zCVec _zPexEqScaleRoots(zCVec ans, int n, int e)
{
  int i;

  for( i=0; i<n; i++ ){
    zCVecElemNC(ans,i)->re = ldexp( zCVecElemNC(ans,i)->re, e );
    zCVecElemNC(ans,i)->im = ldexp( zCVecElemNC(ans,i)->im, e );
  }
  return ans;
}

/* the value of a polynomial at a complex number z, its derivative (if dval is not the null pointer),
 * and the sum of the absolute values of the terms sum_k |a_k| |z|^k, which is the scale of the rounding
 * error of the value. */
static void _zPexEqCVal(const zPex pex, const zComplex *z, zComplex *val, zComplex *dval, double *bound)
{
  int k;
  double az;
  zComplex tmp;

  az = zComplexAbs( z );
  zComplexCreate( val, zPexCoeff(pex,zPexDim(pex)), 0 );
  if( dval ) zComplexZero( dval );
  *bound = fabs( zPexCoeff(pex,zPexDim(pex)) );
  for( k=zPexDim(pex)-1; k>=0; k-- ){
    if( dval ){
      zComplexCMul( dval, z, &tmp );
      zComplexAdd( &tmp, val, dval );
    }
    zComplexCMul( val, z, &tmp );
    zComplexCreate( val, tmp.re + zPexCoeff(pex,k), tmp.im );
    *bound = *bound * az + fabs( zPexCoeff(pex,k) );
  }
}

/* Bairstow-Hitchcock's method. */

/* coefficients of the quotient of a polynomial pex of degree n divided by x^2 + p x + q, which are
 * stored in b[2], ..., b[n], while b[0] and b[1] are the coefficients of the remainder. */
static void _zPexBHQuotient(const zPex pex, int n, double p, double q, double *b)
{
  int j;

  b[n]   = zPexCoeff(pex,n);
  b[n-1] = zPexCoeff(pex,n-1) - p * b[n];
  for( j=2; j<=n; j++ )
    b[n-j] = zPexCoeff(pex,n-j) - p * b[n-j+1] - q * b[n-j+2];
}

/* iterative calculation of quadratic factor x^2 + p x + q of a polynomial pex of degree n based on
 * Bairstow-Hitchcock's method. b and c are workspaces with the sizes n+1 and n, respectively.
 * The iteration is regarded as converged if the updates of p and q are not greater than tol times the
 * scale of the roots of the factor m = |p| + sqrt(|q|) and its square, respectively.
 * It returns the true value if it converges within iter steps. */
static bool _zPexBH(double *p, double *q, const zPex pex, int n, double tol, int iter, double *b, double *c)
{
  double dp, dq, d, m;
  int i, j;

  for( i=0; i<iter; i++ ){
    /* coefficients of quontient */
    _zPexBHQuotient( pex, n, *p, *q, b );
    /* gradient of divider coefficients */
    c[n-1] = b[n];
    c[n-2] = b[n-1] - *p * c[n-1];
    for( j=n-3; j>=0; j-- )
      c[j] = b[j+1] - *p * c[j+1] - *q * c[j+2];
    /* delta of divider coefficients (based on Newton=Raphson's method) */
    d  = 1.0 / ( c[1] * c[1] +( b[1] - c[0] ) * c[2] );
    dp = ( b[1] * c[1] - b[0] * c[2] ) * d;
    dq = ( b[0] * c[1] + b[1] * ( b[1] - c[0] ) ) * d;
    /* termination judgement relative to the scale of the roots */
    m = fabs( *p ) + sqrt( fabs( *q ) );
    if( fabs( dp ) <= tol * m && fabs( dq ) <= tol * m * m ) return true;
    if( !zIsFinite( dp ) || !zIsFinite( dq ) ) return false; /* breakdown */
    *p += dp;
    *q += dq;
  }
  return false;
}

/* deflate a polynomial pex of degree n by a quadratic factor x^2 + p x + q. The quotient is computed
 * from the highest-order coefficient (forward deflation) if the roots of the factor are smaller than
 * the geometric mean of the roots of pex, or from the lowest-order one (backward deflation) otherwise,
 * since the forward and backward deflations are stable for small and large roots, respectively. */
static void _zPexBHDeflate(zPex pex, int n, double p, double q, double *b)
{
  int j;

  if( fabs( q ) <= pow( fabs( zPexCoeff(pex,0) / zPexCoeff(pex,n) ), 2.0/n ) ){ /* forward deflation */
    _zPexBHQuotient( pex, n, p, q, b );
    for( j=2; j<=n; j++ ) zPexSetCoeff( pex, j-2, b[j] );
  } else{ /* backward deflation, in which the k-th coefficient of the quotient c_k is computed from
           * a_k = q c_k + p c_{k-1} + c_{k-2}, and overwrites a_k after a_k is used. */
    zPexSetCoeff( pex, 0, zPexCoeff(pex,0) / q );
    zPexSetCoeff( pex, 1, ( zPexCoeff(pex,1) - p * zPexCoeff(pex,0) ) / q );
    for( j=2; j<=n-2; j++ )
      zPexSetCoeff( pex, j, ( zPexCoeff(pex,j) - p * zPexCoeff(pex,j-1) - zPexCoeff(pex,j-2) ) / q );
  }
}

/* numerical solution of a scaled polynomial equation pex=0 based on Bairstow-Hitchcock's method.
 * (destructive) If the original polynomial org is given, each quadratic factor found for the deflated
 * polynomial is refined by Bairstow-Hitchcock's method for org before deflation, in order to prevent
 * the accumulation of errors through the deflations. org must not have zero roots, since the factor
 * could converge to x^2, which has already been removed from pex. */
#define ZM_PEX_EQ_BH_TRIAL 8
static zCVec _zPexBHDST(zPex pex, const zPex org, zCVec ans, double tol, int iter)
{
  int n, k;
  double p, q, pr, qr, *b, *c;
  bool converged = false;

  ZITERINIT( iter );
  b = zAlloc( double, zPexDim(pex) + 1 );
  c = zAlloc( double, zPexDim(pex) + 1 );
  if( !b || !c ){
    ans = NULL;
    goto TERMINATE;
  }
  for( n=zPexDim(pex); n>0; n-=2 ){
    while( n > 0 && zPexCoeff(pex,0) == 0 ){ /* a zero root */
      zComplexZero( zCVecElemNC(ans,--n) );
      for( k=0; k<=n; k++ ) zPexSetCoeff( pex, k, zPexCoeff(pex,k+1) );
    }
    switch( n ){
    case 0: /* no more root */
      goto TERMINATE;
    case 1: /* linear equation */
      zComplexCreate( zCVecElemNC(ans,0), -zPexCoeff(pex,0)/zPexCoeff(pex,1), 0 );
      goto TERMINATE;
    case 2: /* quadratic equation */
      p = zPexCoeff(pex,1) / zPexCoeff(pex,2);
      q = zPexCoeff(pex,0) / zPexCoeff(pex,2);
      break;
    default:
      /* division by quadratic expression based on Bairstow's method, where the initial factor is
       * chosen as one with a pair of roots on the unit circle at different angles until it converges,
       * since the roots of the scaled polynomial are within a circle of radius about one. */
      for( k=0; k<ZM_PEX_EQ_BH_TRIAL; k++ ){
        p = -2 * cos( zPI * ( k + 0.5 ) / ZM_PEX_EQ_BH_TRIAL );
        q = 1;
        if( ( converged = _zPexBH( &p, &q, pex, n, tol, iter, b, c ) ) ) break;
      }
      if( !converged ){
        ZITERWARN( iter );
        ans = NULL;
        goto TERMINATE;
      }
    }
    if( org && zPexDim(org) > n ){ /* refinement for the original polynomial */
      pr = p;
      qr = q;
      if( _zPexBH( &pr, &qr, org, zPexDim(org), tol, iter, b, c ) ){
        p = pr;
        q = qr;
      }
    }
    /* append solution of 'x^2 + p x + q = 0' */
    zQESolve( 1, p, q, zCVecElemNC(ans,n-2) );
    if( n > 2 ) _zPexBHDeflate( pex, n, p, q, b );
  }
 TERMINATE:
  zFree( b );
  zFree( c );
  return ans;
}

/* polish a root of a polynomial by Newton=Raphson's method on the polynomial, in order to remove the
 * errors accumulated through the deflations. A step is adopted only if it decreases the residual, so
 * that the iteration terminates when the residual reaches the level of rounding errors. */
#define ZM_PEX_EQ_POLISH_ITER 5
static void _zPexEqPolish(const zPex pex, zComplex *z)
{
  zComplex val, dval, dz, z_new, val_new;
  double bound;
  int i;

  _zPexEqCVal( pex, z, &val, &dval, &bound );
  for( i=0; i<ZM_PEX_EQ_POLISH_ITER; i++ ){
    if( zComplexAbs( &val ) == 0 ) return; /* an exact root */
    if( zComplexAbs( &dval ) == 0 ) return; /* a multiple root */
    zComplexCDiv( &val, &dval, &dz );
    zComplexSub( z, &dz, &z_new );
    _zPexEqCVal( pex, &z_new, &val_new, &dval, &bound );
    if( !( zComplexAbs( &val_new ) < zComplexAbs( &val ) ) ) return;
    zComplexCopy( &z_new, z );
    zComplexCopy( &val_new, &val );
  }
}

/* polish the roots of a polynomial with real coefficients, where a complex root with a positive
 * imaginary part is polished and its conjugate is replaced by the conjugate of the result. */
static void _zPexEqPolishRoots(const zPex pex, zCVec ans, int n)
{
  int i, j, jmin;
  double d, dmin;
  zComplex c;

  for( i=0; i<n; i++ ) /* an exact zero root is not polished, since it is removed from pex */
    if( zCVecElemNC(ans,i)->im >= 0 && !( zCVecElemNC(ans,i)->re == 0 && zCVecElemNC(ans,i)->im == 0 ) ) _zPexEqPolish( pex, zCVecElemNC(ans,i) );
  for( i=0; i<n; i++ ){ /* the conjugate of a complex root with a negative imaginary part */
    if( !( zCVecElemNC(ans,i)->im < 0 ) ) continue;
    for( dmin=HUGE_VAL, jmin=-1, j=0; j<n; j++ ){
      if( !( zCVecElemNC(ans,j)->im > 0 ) ) continue;
      zComplexConj( zCVecElemNC(ans,j), &c );
      zComplexSub( &c, zCVecElemNC(ans,i), &c );
      if( ( d = zComplexAbs( &c ) ) < dmin ){
        dmin = d;
        jmin = j;
      }
    }
    if( jmin >= 0 ) zComplexConj( zCVecElemNC(ans,jmin), zCVecElemNC(ans,i) );
    else _zPexEqPolish( pex, zCVecElemNC(ans,i) );
  }
}

/* check the size of the vector to store the roots of a polynomial. */
static bool _zPexEqCheckSize(const zPex pex, const zCVec ans)
{
  char buf[BUFSIZ];

  if( zPexDim(pex) > (int)zCVecSizeNC(ans) ){
    ZRUNERROR( ZM_ERR_PEX_EQ_SIZEMISMATCH, zCVecSizeNC(ans), zI2AOrdinal(zPexDim(pex),buf,BUFSIZ) );
    return false;
  }
  return true;
}

/* numerical solution of polynomial equation based on Bairstow-Hitchcock's method.
 * (destructive) */
zCVec zPexBHDST(zPex pex, zCVec ans, double tol, int iter)
{
  int e;

  if( !_zPexEqCheckSize( pex, ans ) ) return NULL;
  if( tol <= 0 ) tol = zTOL;
  e = _zPexEqScaleExp( pex );
  _zPexEqScale( pex, e, pex );
  if( !_zPexBHDST( pex, NULL, ans, tol, iter ) ) return NULL;
  return _zPexEqScaleRoots( ans, zVecSizeNC(pex)-1, e );
}

/* numerical solution of polynomial equation based on Bairstow-Hitchcock's method. */
zCVec zPexBH(zPex pex, zCVec ans, double tol, int iter)
{
  zPex as, ac;
  int e, k;

  if( !_zPexEqCheckSize( pex, ans ) ) return NULL;
  if( tol <= 0 ) tol = zTOL;
  as = zPexAlloc( zPexDim(pex) );
  ac = zPexAlloc( zPexDim(pex) );
  if( !as || !ac ){
    ans = NULL;
    goto TERMINATE;
  }
  e = _zPexEqScaleExp( pex );
  _zPexEqScale( pex, e, as );
  zVecCopyNC( as, ac );
  while( zPexDim(as) > 0 && zPexCoeff(as,0) == 0 ){ /* zero roots are removed for refinement */
    for( k=0; k<zPexDim(as); k++ ) zPexSetCoeff( as, k, zPexCoeff(as,k+1) );
    zVecSetSizeNC( as, zPexDim(as) );
  }
  if( ( ans = _zPexBHDST( ac, as, ans, tol, iter ) ) ){
    /* the roots are polished on the scaled polynomial before deflation */
    _zPexEqPolishRoots( as, ans, zPexDim(pex) );
    _zPexEqScaleRoots( ans, zPexDim(pex), e );
  }
 TERMINATE:
  zPexFree( as );
  zPexFree( ac );
  return ans;
}

/* Durand-Kerner-Aberth's method */

/* numerical solution of polynomial equation based on Durand-Kerner-Aberth's method.
 * The iteration is regarded as converged if the residuals at all the roots are not greater than tol
 * times the scales of the rounding errors of the values of the polynomial, after which the correction
 * continues while the residuals decrease, in order to improve the accuracy of ill-conditioned roots. */
#define ZM_PEX_EQ_DKA_EXTRA 10
zCVec zPexDKA(zPex pex, zCVec ans, double tol, int iter)
{
  int k, e, extra = 0;
  uint i, j, n;
  bool converged = false;
  double s, r, r0, res, res_max, res_prev = HUGE_VAL;
  zComplex *p, tmp1, tmp2;
  zPex as, b;

  if( !_zPexEqCheckSize( pex, ans ) ) return NULL;
  n = zPexDim(pex);
  p = zAlloc( zComplex, n );
  as = zPexAlloc( n );
  b = zPexAlloc( n );
  if( !p || !as || !b ){
    ans = NULL;
    goto TERMINATE;
  }
  if( tol <= 0 ) tol = zTOL;
  e = _zPexEqScaleExp( pex );
  _zPexEqScale( pex, e, as );
  /* Aberth-Iri's initialization */
  s = zPexCoeff(as,n-1) / ( zPexCoeff(as,n) * n );
  zPexModulo( as, s, b );
  for( r0=0, i=1; i<=n; i++ ){
    r = pow( n*fabs(zPexCoeffHigh(b,i)/zPexCoeff(b,n)), 1.0/(double)i );
    if( r0 < r ) r0 = r;
  }
  for( i=0; i<n; i++ ){
    zComplexCreatePolar( zCVecElemNC(ans,i), r0, ( zPIx2*i + 1.5 )/n );
    zCVecElemNC(ans,i)->re -= s;
  }
  if( r0 == 0 ) goto SUCCESS; /* pex(x) = a_n (x+s)^n, the roots of which are all -s */
  /* Durand's method with Kerner's interpretation */
  ZITERINIT( iter );
  for( k=0; k<iter; k++ ){
    for( res_max=0, i=0; i<n; i++ ){
      _zPexEqCVal( as, zCVecElemNC(ans,i), &p[i], NULL, &r );
      if( ( res = r > 0 ? zComplexAbs( &p[i] ) / r : zComplexAbs( &p[i] ) ) > res_max || zIsNan( res ) ) res_max = res;
    }
    if( !converged && res_max <= tol ) converged = true;
    if( converged ){ /* the correction continues while the residual decreases to the level of rounding errors */
      if( !( res_max < res_prev ) || ++extra > ZM_PEX_EQ_DKA_EXTRA ) goto SUCCESS;
      res_prev = res_max;
    }
    for( i=0; i<n; i++ ){
      for( j=0; j<n; j++ ){
        if( j == i ) continue;
        zComplexSub( zCVecElemNC(ans,i), zCVecElemNC(ans,j), &tmp1 );
        zComplexCDiv( &p[i], &tmp1, &tmp2 );
        zComplexCopy( &tmp2, &p[i] );
      }
      zComplexDiv( &p[i], zPexCoeff(as,n), &p[i] );
    }
    for( i=0; i<n; i++ )
      zComplexSub( zCVecElemNC(ans,i), &p[i], zCVecElemNC(ans,i) );
  }
  ZITERWARN( iter );
 SUCCESS:
  _zPexEqScaleRoots( ans, n, e );
 TERMINATE:
  zFree( p );
  zPexFree( as );
  zPexFree( b );
  return ans;
}
