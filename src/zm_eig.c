/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_eig - eigenvalue analysis.
 */

#include <zm/zm_eig.h>

/* Householder transformation to Hessian matrices */

/* Householder transformation (destructive). */
void zHouseholder(zMat m, zMat p, int from, int to, zVec u, zVec v, zVec w)
{
  int i, j;
  double s, a;

  /* v = 2Mu, w = 2M^Tu */
  for( i=0; i<zVecSizeNC(v); i++ ){
    zVecSetElemNC( v, i, 0 );
    zVecSetElemNC( w, i, 0 );
    for( j=from; j<to; j++ ){
      zVecElemNC(v,i) += zMatElemNC(m,i,j)*zVecElemNC(u,j);
      zVecElemNC(w,i) += zMatElemNC(m,j,i)*zVecElemNC(u,j);
    }
    zVecElemNC(v,i) *= 2;
    zVecElemNC(w,i) *= 2;
  }
  /* s = v^T u */
  a = zRawVecInnerProd( &zVecElemNC(v,from), &zVecElemNC(u,from), to-from );
  /* v -= su, w -= su */
  for( i=from; i<to; i++ ){
    s = a * zVecElemNC(u,i);
    zVecElemNC(v,i) -= s;
    zVecElemNC(w,i) -= s;
  }
  /* M -= vu^T+uw^T */
  for( i=0; i<zVecSizeNC(v); i++ )
    for( j=from; j<to; j++ ){
      zMatElemNC(m,i,j) -= zVecElemNC(v,i)*zVecElemNC(u,j);
      zMatElemNC(m,j,i) -= zVecElemNC(u,j)*zVecElemNC(w,i);
    }
  /* transportation matrix */
  if( !p ) return;
  zVecZero( v );
  for( i=0; i<zVecSizeNC(v); i++ ){
    for( j=from; j<to; j++ )
      zVecElemNC(v,i) += zMatElemNC(p,i,j)*zVecElemNC(u,j);
    zVecElemNC(v,i) *= 2;
  }
  for( i=0; i<zVecSizeNC(v); i++ )
    for( j=from; j<to; j++ )
      zMatElemNC(p,i,j) -= zVecElemNC(v,i)*zVecElemNC(u,j);
}

/* projection vector of Householder transformation. */
zVec zHouseholderVec(zMat m, int col, int from, int to, zVec u)
{
  int i;
  double s;

  for( s=0, i=from; i<to; i++ )
    s += zSqr( zVecSetElemNC( u, i, zMatElemNC(m,i,col) ) );
  if( zIsTiny( s ) ) return u;
  zVecElemNC(u,from) -= ( s = -zSgn(zMatElemNC(m,from,col)) * sqrt(s) );
  zRawVecDivDRC( &zVecElemNC(u,from), sqrt(2*(s-zMatElemNC(m,from,col))*s), to-from );
  return u;
}

/* Hessenberg matrix using Householder transformation. */
zMat zHess(zMat m, zMat h, zMat p)
{
  int n, n2, i;
  zVec u, v, w;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  /* target hessian matrix */
  if( h ){
    if( !zMatSizeIsEqual( m, h ) ){
      ZRUNERROR( ZM_ERR_SIZMIS_MAT );
      return NULL;
    }
    zMatCopyNC( m, h );
  } else
    h = m;
  /* transformation matrix */
  if( p && !zMatSizeIsEqual( h, p ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }

  u = zVecAlloc( zMatRowSizeNC(h) );
  v = zVecAlloc( zMatRowSizeNC(h) );
  w = zVecAlloc( zMatRowSizeNC(h) );
  if( !u || !v || !w ) goto TERMINATE;
  if( p ) zMatIdentNC( p );
  n2 = zMatRowSizeNC(h) - 2;
  for( n=0; n<n2; n++ ){
 RETRY:
    zHouseholderVec( h, n, n+1, zVecSizeNC(u), u );
    zHouseholder( h, p, n+1, zVecSizeNC(u), u, v, w );
    for( i=n+2; i<zMatRowSizeNC(h); i++ ) /* need refinement? */
      if( !zIsTiny( zMatElemNC(h,i,n) ) ) goto RETRY;
  }
 TERMINATE:
  zVecFreeAO( 3, u, v, w );
  return h;
}

/* eigenvalue analysis by double QR method and inverse iteration */

/* Householder's method to transform a Hesseian matrix to a block uppertriangle matrix. */
static void _zEigDQRHouseholder(zMat a, int r, int c, double g1, double g2, double g3)
{
  double s, den, alpha;
  int i, c1, c2, n;

  s = sqrt( g1*g1 + g2*g2 + g3*g3 );
  den = g1 + ( g1 > 0 ? s : -s );
  g2 /= den;
  g3 /= den;
  alpha = 2.0 / ( 1.0 + g2*g2 + g3*g3 );

  c1 = c + 1;
  c2 = c + 2;
  for( i=zMax(c-1,0); i<=r; i++ ){
    s = (zMatElemNC(a,c,i)+g2*zMatElemNC(a,c1,i)+g3*zMatElemNC(a,c2,i))*alpha;
    zMatElemNC(a,c ,i) -= s;
    zMatElemNC(a,c1,i) -= s * g2;
    zMatElemNC(a,c2,i) -= s * g3;
  }
  n = zLimit( c+3, r, zMatRowSizeNC(a)-1 );
  for( i=0; i<=n; i++ ){
    s = (zMatElemNC(a,i,c)+g2*zMatElemNC(a,i,c1)+g3*zMatElemNC(a,i,c2))*alpha;
    zMatElemNC(a,i,c ) -= s;
    zMatElemNC(a,i,c1) -= s * g2;
    zMatElemNC(a,i,c2) -= s * g3;
  }
}

/* double QR method. */
bool zEigDQR(zMat m, zComplex z[], int iter)
{
  int  i, j, r, r1;
  double a0, a1, b, c;
  zMat a;

  if( !zMatIsSqr(m) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return false;
  }
  if( !( a = zMatClone(m) ) ){
    ZALLOCERROR();
    return false;
  }

  zHess( a, NULL, NULL );
  for( r=zMatRowSizeNC(a)-1; r>=0; ){
    if( r == 0 ){
      zComplexCreate( &z[0], zMatBuf(a)[0], 0.0 );
      break;
    }
    ZITERINIT( iter );
    for( i=0; i<iter; i++ ){
      r1 = r - 1;
      if( zIsTiny( zMatElemNC(a,r,r1) ) ){
        zComplexCreate( &z[r], zMatElemNC(a,r,r), 0.0 );
        r = r1;
        break;
      }
      b =-zMatElemNC(a,r1,r1)-zMatElemNC(a,r,r);
      c = zMatElemNC(a,r1,r1)*zMatElemNC(a,r,r)
         -zMatElemNC(a,r1,r) *zMatElemNC(a,r,r1);
      if( r1 == 0 || zIsTiny( zMatElemNC(a,r1,r1-1) ) ){
        zQESolve( 1.0, b, c, &z[r1] );
        r -= 2;
        break;
      }
      b += ( a0 = zMatBuf(a)[0] );
      a1 = zMatElemNC(a,1,0);
      _zEigDQRHouseholder( a, r, 0,
        zMatElemNC(a,0,1)*a1+a0*b+c, (zMatElemNC(a,1,1)+b)*a1, zMatElemNC(a,2,1)*a1 );
      for( j=1; j<=r1-1; j++ )
        _zEigDQRHouseholder( a, r, j,
          zMatElemNC(a,j,j-1), zMatElemNC(a,j+1,j-1), zMatElemNC(a,j+2,j-1) );
      _zEigDQRHouseholder( a, r, j,
        zMatElemNC(a,r1,r1-1), zMatElemNC(a,r ,r1-1), 0 );
    }
  }
  zMatFree( a );
  return true;
}

/* calculate the dominant eigenvalue. */
double zEigPower(zMat a, zVec evec, int iter)
{
  zVec ev, err;
  double eig = 1.0;
  int i;

  err = zVecAlloc( zVecSizeNC(evec) );
  ev = zVecAlloc( zVecSizeNC(evec) );
  if( !err || !ev ) goto TERMINATE;
  zVecRandUniform( evec, -1, 1 );

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecNormalize( evec, ev );
    zMulMatVecNC( a, ev, evec );
    eig = zVecInnerProd( ev, evec );
    zVecCatNC( evec, -eig, ev, err );
    if( zVecIsTiny( err ) ) goto TERMINATE;
  }
  ZITERWARN( iter );

 TERMINATE:
  zVecDivDRC( evec, eig );
  zVecFree( ev );
  zVecFree( err );
  return eig;
}

/* calculate the minimal eigenvalue. */
double zEigPowerInv(zMat a, zVec evec, int iter)
{
  zMat ai;
  double eig = 0;

  if( !zMatIsSqr(a) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return 0;
  }
  if( !( ai = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a) ) ) )
    return 0;
  if( zMatInv( a, ai ) )
    eig = zEigPower( ai, evec, iter );
  else
    ZRUNERROR( ZM_ERR_LE_SINGULAR );

  zMatFree( ai );
  return 1.0 / eig;
}

/* compute an eigenvector for a real-number eigenvalue. */
static int _zEigVecReal(zMat m, double eig, zCVec eigv, int iter)
{
  zMat ms, b;
  zVec eigv_r;
  double shift = zTOL * 1.0e1;
  int ret = 0;

  ms = zMatClone( m );
  b = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  eigv_r = zVecAlloc( zCVecSizeNC(eigv) );
  if( !ms || !b || !eigv_r ){
    ZALLOCERROR();
    ret = -1;
    goto TERMINATE;
  }

  zVecRandUniform( eigv_r, -1, 1 );
  zMatShift( ms, -eig - shift );
  while( !zMatInv( ms, b ) )
    zMatShift( ms, -( shift*=10 ) );
  zEigPower( b, eigv_r, iter );
  zVec2CVec( eigv_r, eigv );

 TERMINATE:
  zMatFree( ms );
  zMatFree( b );
  zVecFree( eigv_r );
  return ret;
}

/* compute eigenvectors for complex-number eigenvalues. */
static int _zEigVecComplex(zMat m, zComplex *eig, zCVec eigv1, zCVec eigv2, int iter)
{
  zMat ms, b;
  zVec eigv;
  double shift = zTOL * 1.0e1, s;
  int i, n, ret = 0;

  n = 2 * zMatRowSizeNC(m);
  ms = zMatAllocSqr( n );
  b = zMatAllocSqr( n );
  eigv = zVecAlloc( n );
  if( !ms || !b || !eigv ){
    ret = -1;
    goto TERMINATE;
  }

  zMatPutNC( ms, 0, 0, m );
  zMatPutNC( ms, zMatRowSizeNC(m), zMatColSizeNC(m), m );
  s = eig->re + shift;
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    zMatElemNC(ms,i,i) -= s;
    zMatSetElemNC( ms, i+zMatRowSizeNC(m), i,-eig->im );
    zMatSetElemNC( ms, i, i+zMatColSizeNC(m), eig->im );
    zMatElemNC(ms,i+zMatRowSizeNC(m),i+zMatColSizeNC(m)) -= s;
  }

  zVecRandUniform( eigv, -1, 1 );
  while( !zMatInv( ms, b ) )
    zMatShift( ms, -( shift *= 10 ) );
  zEigPower( b, eigv, iter );

  for( i=0; i<zCVecSizeNC(eigv1); i++ ){
    zComplexCreate( zCVecElem(eigv1,i),
      zVecElemNC(eigv,i), zVecElemNC(eigv,i+zCVecSizeNC(eigv1)) );
    zComplexCreate( zCVecElem(eigv2,i),
      zVecElemNC(eigv,i),-zVecElemNC(eigv,i+zCVecSizeNC(eigv1)) );
  }
  zCVecNormalizeDRC( eigv1 );
  zCVecNormalizeDRC( eigv2 );

 TERMINATE:
  zMatFree( ms );
  zMatFree( b );
  zVecFree( eigv );
  return ret;
}

/* compute eigenvalue and eigenvector sets of a square matrix. */
#define Z_EIGSYS_ITER_NUM 10
int zEigSystem(zMat m, zComplex eig[], zCVec eigv[], int iter)
{
  int i;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return -1;
  }
  /* eigen values by double QR method */
  zEigDQR( m, eig, iter );
  /* eigen vectors by inverse iteration method */
  for( i=0; i<zMatRowSizeNC(m); )
    if( eig[i].im == 0 ){
      if( _zEigVecReal( m, eig[i].re, eigv[i], Z_EIGSYS_ITER_NUM ) < 0 ){
        ZALLOCERROR();
        break;
      }
      i++;
    } else{
      if( _zEigVecComplex( m, &eig[i], eigv[i], eigv[i+1], Z_EIGSYS_ITER_NUM ) < 0 ){
        ZALLOCERROR();
        break;
      }
      i += 2;
    }
  /* return the number of eigenvalues and eigen vectors */
  return i;
}

/* diagonalization by bisection method */

/* compute the range in which all eigenvalues exist based on Gerschgorin's theorem. */
static void _zEigSymBisecRange(zMat m, double *emin, double *emax)
{
  int i;
  double ml, mr, e;

  *emin = *emax = 0;
  for( ml=0, i=0; i<zMatRowSizeNC(m); i++, ml=mr ){
    mr = i < zMatRowSizeNC(m) ? fabs( zMatElemNC(m,i,i+1) ) : 0;
    if( ( e = zMatElemNC(m,i,i) + ( ml + mr ) ) > *emax ) *emax = e;
    if( ( e = zMatElemNC(m,i,i) - ( ml + mr ) ) < *emin ) *emin = e;
  }
}

/* Sturm's function series. */
static int _zEigSysBisecSturm(zMat a, double e)
{
  int i, n = 0;
  double r, p1, p2;

  p1 = e - zMatBuf(a)[0];
  if( p1 <= 0 ) n++;
  for( i=1; i<zMatRowSizeNC(a); i++ ){
    r = zMatElemNC(a,i,i-1);
    p2 = e - zMatElemNC(a,i,i) - r*r / ( p1 == 0 ? zTOL : p1 );
    if( p2 < 0 ) n++;
    p1 = p2;
  }
  return n;
}

/* compute one eigenvalue by bisection shooting. */
static void _zEigSymBisecOne(zMat a, double emin, int nmin, double emax, int nmax, double *eig)
{
  double e, d, d_old;
  int n;

  for( d_old=HUGE_VAL; ; d_old=d ){
    e = emin + 0.5 * ( d = emax - emin );
    if( d == d_old ) break;
    n = _zEigSysBisecSturm( a, e );
    if( n == nmin ) emin = e;
    else emax = e;
  }
  *eig = e;
}

/* compute all eigenvalues recursively by bisection method. */
static void _zEigSymBisecRec(zMat a, double emin, int nmin, double emax, int nmax, double *eig)
{
  int n;
  double e;

  if( nmin - nmax == 1 )
    return _zEigSymBisecOne( a, emin, nmin, emax, nmax, eig );

  n = _zEigSysBisecSturm( a, ( e = 0.5 * ( emax + emin ) ) );
  if( n > nmax )
    _zEigSymBisecRec( a, e, n, emax, nmax, &eig[0] );
  if( n < nmin )
    _zEigSymBisecRec( a, emin, nmin, e, n, &eig[n-nmax] );
}

/* compute a unitary matrix consisting of eigenvectors for real-number eigenvalues. */
static bool _zEigSymBisecR(zMat m, zVec eig, zMat r, int iter)
{
  int i;
  zMat ms, b;
  zVec eigv;
  double shift = zTOL * 1.0e1;
  bool ret = true;

  ms = zMatAllocSqr( zMatRowSizeNC(m) );
  b = zMatAllocSqr( zMatRowSizeNC(m) );
  eigv = zVecAlloc( zMatRowSizeNC(m) );
  if( !ms || !b || !eigv ){
    ZALLOCERROR();
    ret = false;
    goto TERMINATE;
  }

  for( i=0; i<zVecSizeNC(eig); i++ ){
    zMatCopyNC( m, ms );
    zVecRandUniform( eigv, -1, 1 );
    zMatShift( ms, -zVecElemNC(eig,i)-shift );
    while( !zMatInv( ms, b ) )
      zMatShift( ms, -( shift*=10 ) );
    zEigPower( b, eigv, iter );
    zMatPutCol( r, i, eigv );
  }

 TERMINATE:
  zMatFree( ms );
  zMatFree( b );
  zVecFree( eigv );
  return ret;
}

/* compute eigenvalues and eigenvectors based on bisection method
 * (by J. W. Givens 1954).
 */
zVec zEigSymBisec(zMat m, zVec eig, zMat r)
{
  double emin, emax;
  int nmin, nmax;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( !r || !zMatSizeIsEqual( m, r ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  if( !eig || !zMatRowVecSizeIsEqual( m, eig ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  /* eigenvalues */
  zMatCopyNC( m, r ); /* temporary working space */
  zHess( r, NULL, NULL );
  _zEigSymBisecRange( r, &emin, &emax );
  nmin = _zEigSysBisecSturm( r, emin );
  nmax = _zEigSysBisecSturm( r, emax );
  _zEigSymBisecRec( r, emin, nmin, emax, nmax, zVecBuf(eig) );
  _zEigSymBisecR( m, eig, r, Z_EIGSYS_ITER_NUM );
  return eig;
}

/* diagonalization by Jacobi's method */

/* shift diagonal values of a matrix to accelerate diagonalization. */
static double _zEigSymJacobiShift(zMat m, double *shift)
{
  zComplex c_shift[2];
  double b, c;

  b = - zMatBuf(m)[0] - zMatElemNC(m,1,1);
  c = zMatBuf(m)[0] * zMatElemNC(m,1,1)
      - zMatElemNC(m,1,0) * zMatElemNC(m,0,1);
  zQESolve( 1, b, c, c_shift );
  c = fabs( c_shift[0].re - zMatBuf(m)[0] )
    < fabs( c_shift[1].re - zMatBuf(m)[0] ) ?
    c_shift[0].re : c_shift[1].re;
  zMatShift( m, -c );
  return shift ? ( *shift = c ) : c;
}

/* transform a matrix by Jacobi's rotation based on Wilkinson's formula. */
static void _zEigSymJacobiRot(zMat m, zMat r, int i, int j)
{
  int k;
  double as, ad, ti, c, s;
  double tmp1, tmp2;

  as = 0.5 * ( zMatElemNC(m,i,i) + zMatElemNC(m,j,j) );
  ad = 0.5 * ( zMatElemNC(m,i,i) - zMatElemNC(m,j,j) );
  ti = sqrt( ad*ad + zSqr( zMatElemNC(m,i,j) ) );
  if( ad < 0 ) ti = -ti;
  c = sqrt( 0.5 + 0.5*ad/ti );
  s = 0.5 * zMatElemNC(m,i,j) / ( ti * c );

  zMatSetElemNC( m, i, i, as + ti );
  zMatSetElemNC( m, j, j, as - ti );
  zMatSetElemNC( m, i, j, 0 );
  zMatSetElemNC( m, j, i, 0 );
  for( k=0; k<zMatRowSizeNC(m); k++ ){
    /* update of transformation matrix */
    tmp1 = zMatElemNC( r, k, i );
    tmp2 = zMatElemNC( r, k, j );
    zMatSetElemNC( r, k, i, c * tmp1 + s * tmp2 );
    zMatSetElemNC( r, k, j,-s * tmp1 + c * tmp2 );
    /* update of eigenmatrix */
    if( k == i || k == j ) continue;
    tmp1 = zMatElemNC( m, i, k );
    tmp2 = zMatElemNC( m, j, k );
    zMatSetElemNC( m, i, k, c * tmp1 + s * tmp2 );
    zMatSetElemNC( m, j, k,-s * tmp1 + c * tmp2 );
    zMatSetElemNC( m, k, i, zMatElemNC( m, i, k ) );
    zMatSetElemNC( m, k, j, zMatElemNC( m, j, k ) );
  }
}

/* diagonalize a symmetric matrix by Jacobi's method.
 * when m is not symmetric, it doesn't work as expected.
 */
zVec zEigSymJacobi(zMat m, zVec eig, zMat r)
{
  int i, j, n = 0;
  double shift;
  bool is_complete;
  zMat d;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  /* transformation matrix */
  if( !r || !zMatSizeIsEqual( m, r ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  if( !eig || !zMatRowVecSizeIsEqual( m, eig ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( !( d = zMatClone( m ) ) ){ /* diagonal matrix */
    ZALLOCERROR();
    return NULL;
  }

  zMatIdent( r );
  _zEigSymJacobiShift( d, &shift );
  do{
    is_complete = true;
    for( i=1; i<zMatRowSizeNC(d); i++ )
      for( j=0; j<i; j++ ){
        if( zIsTiny( zMatElemNC(d,i,j) ) ) continue;
	is_complete = false;
	/* iterative elimination of non-diagonal components */
        _zEigSymJacobiRot( d, r, i, j );
      }
    if( n++ > Z_MAX_ITER_NUM ){
      ZITERWARN( Z_MAX_ITER_NUM );
      goto TERMINATE;
    }
  } while( !is_complete );
 TERMINATE:
  for( i=0; i<zMatRowSizeNC(m); i++ )
    zVecSetElemNC( eig, i, zMatElemNC(d,i,i) + shift );
  zMatFree( d );
  return eig;
}

/* singular value decomposition */

/* sort singular values and corresponding bases. */
static int _zSVDSort(zVec sv, zMat u)
{
  int i, im;

  for( i=0; i<zVecSizeNC(sv); i++ ){
    if( zIsTiny( zDataMax( zVecBuf(sv)+i, zVecSizeNC(sv)-i, &im ) ) )
      break;
    im += i;
    zVecSwapNC( sv, i, im );
    zMatSwapColNC( u, i, im );
  }
  return i;
}

/* singular value decomposition. */
int zSVD(zMat m, zVec sv, zMat u, zMat v)
{
  int i, j, rank;
  zMat c, w = NULL;
  double s;

  if( !( c = zMatAllocSqr( zMatRowSizeNC(m) ) ) ){
    ZALLOCERROR();
    rank = -1;
    goto TERMINATE;
  }

  zMulMatMatTNC( m, m, c );
  zEigSymJacobi( c, sv, u );
  for( i=0; i<zVecSizeNC(sv); i++ )
    zVecSetElemNC( sv, i, ( zVecElemNC(sv,i) < zTOL ) ? 0 : sqrt(zVecElemNC(sv,i)) );
  rank = _zSVDSort( sv, u );
  if( !( w = zMatAlloc( zMatRowSizeNC(u), rank ) ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  zMatGetNC( u, 0, 0, w );
  zMatSetRowSize( v, rank );

  zMulMatTMatNC( w, m, v );
  for( i=0; i<rank; i++ ){
    s = 1.0 / zVecElemNC(sv,i);
    for( j=0; j<zMatColSizeNC(m); j++ )
      zMatElemNC(v,i,j) *= s;
  }
 TERMINATE:
  zMatFree( w );
  zMatFree( c );
  return rank;
}

/* allocate internal matrix and vector for workspace. */
static bool _zSVMat(zMat m, zMat *r, zVec *e)
{
  if( zMatRowSizeNC(m) > zMatColSizeNC(m) ){
    *r = zMatAllocSqr( zMatColSizeNC(m) );
    *e = zVecAlloc( zMatColSizeNC(m) );
    if( !*r || !*e ){
      ZALLOCERROR();
      return false;
    }
    zMulMatTMat( m, m, *r );
  } else{
    *r = zMatAllocSqr( zMatRowSizeNC(m) );
    *e = zVecAlloc( zMatRowSizeNC(m) );
    if( !*r || !*e ){
      ZALLOCERROR();
      return false;
    }
    zMulMatMatT( m, m, *r );
  }
  zVecRandUniform( *e, -1, 1 );
  return true;
}

/* maximum singular value of a matrix. */
double zSVMax(zMat m)
{
  zMat r;
  zVec e;
  double s = 0;

  if( _zSVMat( m, &r, &e ) )
    s = sqrt( zEigPower( r, e, 0 ) );
  zMatFree( r );
  zVecFree( e );
  return s;
}

/* minimum singular value of a matrix. */
double zSVMin(zMat m)
{
  zMat r;
  zVec e;
  double s = 0;

  if( _zSVMat( m, &r, &e ) )
    s = sqrt( zEigPowerInv( r, e, 0 ) );
  zMatFree( r );
  zVecFree( e );
  return s;
}

/* condition number of a matrix. */
double zMatCondNum(zMat m)
{
  zMat r;
  zVec e;
  double smax, smin;

  if( !_zSVMat( m, &r, &e ) ) return NAN;
  smax = sqrt( zEigPower( r, e, 0 ) );
  smin = sqrt( zEigPowerInv( r, e, 0 ) );
  zMatFree( r );
  zVecFree( e );
  return smin > zTOL ? smax / smin : HUGE_VAL;
}
