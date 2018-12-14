/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_eig - eigenvalue analysis.
 */

#include <zm/zm_eig.h>

/* Householder transformation to Hessian matrices */

/* zHouseholder
 * - Householder transformation (destructive).
 */
void zHouseholder(zMat m, zMat p, int from, int to, zVec u, zVec v, zVec w)
{
  register int i, j;
  double s, a;

  /* v = 2Mu, w = 2M^Tu */
  for( i=0; i<zVecSizeNC(v); i++ ){
    zVecSetElem( v, i, 0 );
    zVecSetElem( w, i, 0 );
    for( j=from; j<to; j++ ){
      zVecElem(v,i) += zMatElem(m,i,j)*zVecElem(u,j);
      zVecElem(w,i) += zMatElem(m,j,i)*zVecElem(u,j);
    }
    zVecElem(v,i) *= 2;
    zVecElem(w,i) *= 2;
  }
  /* s = v^T u */
  a = zRawVecInnerProd( &zVecElem(v,from), &zVecElem(u,from), to-from );
  /* v -= su, w -= su */
  for( i=from; i<to; i++ ){
    s = a * zVecElem(u,i);
    zVecElem(v,i) -= s;
    zVecElem(w,i) -= s;
  }
  /* M -= vu^T+uw^T */
  for( i=0; i<zVecSizeNC(v); i++ )
    for( j=from; j<to; j++ ){
      zMatElem(m,i,j) -= zVecElem(v,i)*zVecElem(u,j);
      zMatElem(m,j,i) -= zVecElem(u,j)*zVecElem(w,i);
    }
  /* transportation matrix */
  if( !p ) return;
  zVecClear( v );
  for( i=0; i<zVecSizeNC(v); i++ ){
    for( j=from; j<to; j++ )
      zVecElem(v,i) += zMatElem(p,i,j)*zVecElem(u,j);
    zVecElem(v,i) *= 2;
  }
  for( i=0; i<zVecSizeNC(v); i++ )
    for( j=from; j<to; j++ )
      zMatElem(p,i,j) -= zVecElem(v,i)*zVecElem(u,j);
}

/* (static)
 * zHouseholderVec
 * - projection vector of Householder transformation.
 */
zVec zHouseholderVec(zMat m, int col, int from, int to, zVec u)
{
  register int i;
  double s;

  for( s=0, i=from; i<to; i++ )
    s += zSqr( zVecSetElem( u, i, zMatElem(m,i,col) ) );
  if( zIsTiny( s ) ) return u;
  zVecElem(u,from) -= ( s = -zSgn(zMatElem(m,from,col)) * sqrt(s) );
  zRawVecDivDRC( &zVecElem(u,from), sqrt(2*(s-zMatElem(m,from,col))*s), to-from );
  return u;
}

/* zHess
 * - Hessenberg matrix using Householder transformation.
 */
zMat zHess(zMat m, zMat h, zMat p)
{
  register int n, n2, i;
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
  if( !u || !v || !w ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  if( p ) zMatIdentNC( p );
  n2 = zMatRowSizeNC(h) - 2;
  for( n=0; n<n2; n++ ){
 RETRY:
    zHouseholderVec( h, n, n+1, zVecSizeNC(u), u );
    zHouseholder( h, p, n+1, zVecSizeNC(u), u, v, w );
    for( i=n+2; i<zMatRowSizeNC(h); i++ ) /* need refinement? */
      if( !zIsTiny( zMatElem(h,i,n) ) ) goto RETRY;
  }
 TERMINATE:
  zVecFreeAO( 3, u, v, w );
  return h;
}

/* eigenvalue analysis by double QR method and inverse iteration */

static void _zEigDQRHouseholder(zMat a, int r, int c, double g1, double g2, double g3);

static int _zEigVecReal(zMat m, double eig, zCVec eigv, int iter);
static int _zEigVecComplex(zMat m, zComplex *eig, zCVec eigv1, zCVec eigv2, int iter);

/* (static)
 * _zEigDQRHouseholder
 * - Householder's method to transform a Hesseian matrix to
 *   a block uppertriangle matrix.
 */
void _zEigDQRHouseholder(zMat a, int r, int c, double g1, double g2, double g3)
{
  double s, den, alpha;
  register int i, c1, c2, n;

  s = sqrt( g1*g1 + g2*g2 + g3*g3 );
  den = g1 + ( g1 > 0 ? s : -s );
  g2 /= den;
  g3 /= den;
  alpha = 2.0 / ( 1.0 + g2*g2 + g3*g3 );

  c1 = c + 1;
  c2 = c + 2;
  for( i=zMax(c-1,0); i<=r; i++ ){
    s = (zMatElem(a,c,i)+g2*zMatElem(a,c1,i)+g3*zMatElem(a,c2,i))*alpha;
    zMatElem(a,c ,i) -= s;
    zMatElem(a,c1,i) -= s * g2;
    zMatElem(a,c2,i) -= s * g3;
  }
  n = zLimit( c+3, r, zMatRowSizeNC(a)-1 );
  for( i=0; i<=n; i++ ){
    s = (zMatElem(a,i,c)+g2*zMatElem(a,i,c1)+g3*zMatElem(a,i,c2))*alpha;
    zMatElem(a,i,c ) -= s;
    zMatElem(a,i,c1) -= s * g2;
    zMatElem(a,i,c2) -= s * g3;
  }
}

/* zEigDQR
 * - double QR method.
 */
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
      zComplexCreate( &z[0], zMatElem(a,0,0), 0.0 );
      break;
    }
    ZITERINIT( iter );
    for( i=0; i<iter; i++ ){
      r1 = r - 1;
      if( zIsTiny( zMatElem(a,r,r1) ) ){
        zComplexCreate( &z[r], zMatElem(a,r,r), 0.0 );
        r = r1;
        break;
      }
      b =-zMatElem(a,r1,r1)-zMatElem(a,r,r);
      c = zMatElem(a,r1,r1)*zMatElem(a,r,r)
         -zMatElem(a,r1,r) *zMatElem(a,r,r1);
      if( r1 == 0 || zIsTiny( zMatElem(a,r1,r1-1) ) ){
        zQESolve( 1.0, b, c, &z[r1] );
        r -= 2;
        break;
      }
      b += ( a0 = zMatElem(a,0,0) );
      a1 = zMatElem(a,1,0);
      _zEigDQRHouseholder( a, r, 0,
        zMatElem(a,0,1)*a1+a0*b+c, (zMatElem(a,1,1)+b)*a1, zMatElem(a,2,1)*a1 );
      for( j=1; j<=r1-1; j++ )
        _zEigDQRHouseholder( a, r, j,
          zMatElem(a,j,j-1), zMatElem(a,j+1,j-1), zMatElem(a,j+2,j-1) );
      _zEigDQRHouseholder( a, r, j,
        zMatElem(a,r1,r1-1), zMatElem(a,r ,r1-1), 0 );
    }
  }
  zMatFree( a );
  return true;
}

/* zEigPower
 * - calculate the dominant eigenvalue.
 */
double zEigPower(zMat a, zVec evec, int iter)
{
  zVec ev, err;
  double eig = 1.0;
  register int i;

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

/* zEigPowerInv
 * - calculate the minimal eigenvalue.
 */
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

/* (static)
 * _zEigVecReal
 * - compute an eigenvector for a real-number eigenvalue.
 */
int _zEigVecReal(zMat m, double eig, zCVec eigv, int iter)
{
  zMat ms, b;
  zVec eigv_r;
  double shift = zTOL * 1.0e1;
  int ret = 0;

  ms = zMatClone( m );
  b = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  eigv_r = zVecAlloc( _zCVecSize(eigv) );
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

/* (static)
 * _zEigVecComplex
 * - compute eigenvectors for complex-number eigenvalues.
 */
int _zEigVecComplex(zMat m, zComplex *eig, zCVec eigv1, zCVec eigv2, int iter)
{
  zMat ms, b;
  zVec eigv;
  double shift = zTOL * 1.0e1, s;
  int ret = 0;
  register int i, n;

  n = 2 * zMatRowSizeNC(m);
  ms = zMatAllocSqr( n );
  b = zMatAllocSqr( n );
  eigv = zVecAlloc( n );
  if( !ms || !b || !eigv ){
    ZALLOCERROR();
    ret = -1;
    goto TERMINATE;
  }

  zMatPutNC( ms, 0, 0, m );
  zMatPutNC( ms, zMatRowSizeNC(m), zMatColSizeNC(m), m );
  s = eig->re + shift;
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    zMatElem(ms,i,i) -= s;
    zMatSetElem( ms, i+zMatRowSizeNC(m), i,-eig->im );
    zMatSetElem( ms, i, i+zMatColSizeNC(m), eig->im );
    zMatElem(ms,i+zMatRowSizeNC(m),i+zMatColSizeNC(m)) -= s;
  }

  zVecRandUniform( eigv, -1, 1 );
  while( !zMatInv( ms, b ) )
    zMatShift( ms, -( shift *= 10 ) );
  zEigPower( b, eigv, iter );

  for( i=0; i<_zCVecSize(eigv1); i++ ){
    zComplexCreate( zCVecElem(eigv1,i),
      zVecElem(eigv,i), zVecElem(eigv,i+_zCVecSize(eigv1)) );
    zComplexCreate( zCVecElem(eigv2,i),
      zVecElem(eigv,i),-zVecElem(eigv,i+_zCVecSize(eigv1)) );
  }
  zCVecNormalizeDRC( eigv1 );
  zCVecNormalizeDRC( eigv2 );

 TERMINATE:
  zMatFree( ms );
  zMatFree( b );
  zVecFree( eigv );
  return ret;
}

/* zEigSystem
 * - compute eigenvalue and eigenvector sets of a square matrix.
 */
#define Z_EIGSYS_ITER_NUM 10
int zEigSystem(zMat m, zComplex eig[], zCVec eigv[], int iter)
{
  register int i;

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

static void _zEigSymBisecRange(zMat m, double *emin, double *emax);
static int _zEigSysBisecSturm(zMat a, double e);
static void _zEigSymBisecOne(zMat a, double emin, int nmin, double emax, int nmax, double *eig);
static void _zEigSymBisecRec(zMat a, double emin, int nmin, double emax, int nmax, double *eig);
static bool _zEigSymBisecR(zMat m, zVec eig, zMat r, int iter);

/* (static)
 * _zEigSymBisecRange
 * - compute the range in which all eigenvalues exist
 *   based on Gerschgorin's theorem.
 */
void _zEigSymBisecRange(zMat m, double *emin, double *emax)
{
  register int i;
  double ml, mr, e;

  *emin = *emax = 0;
  for( ml=0, i=0; i<zMatRowSizeNC(m); i++, ml=mr ){
    mr = i < zMatRowSizeNC(m) ? fabs( zMatElem(m,i,i+1) ) : 0;
    if( ( e = zMatElem(m,i,i) + ( ml + mr ) ) > *emax ) *emax = e;
    if( ( e = zMatElem(m,i,i) - ( ml + mr ) ) < *emin ) *emin = e;
  }
}

/* (static)
 * _zEigSysBisecSturm
 * - Sturm's function series.
 */
int _zEigSysBisecSturm(zMat a, double e)
{
  register int i, n = 0;
  double r, p1, p2;

  p1 = e - zMatElem(a,0,0);
  if( p1 <= 0 ) n++;
  for( i=1; i<zMatRowSizeNC(a); i++ ){
    r = zMatElem(a,i,i-1);
    p2 = e - zMatElem(a,i,i) - r*r / ( p1 == 0 ? zTOL : p1 );
    if( p2 < 0 ) n++;
    p1 = p2;
  }
  return n;
}

/* (static)
 * _zEigSymBisecOne
 * - compute one eigenvalue by bisection shooting.
 */
void _zEigSymBisecOne(zMat a, double emin, int nmin, double emax, int nmax, double *eig)
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

/* (static)
 * _zEigSymBisecRec
 * - compute all eigenvalues recursively by bisection method.
 */
void _zEigSymBisecRec(zMat a, double emin, int nmin, double emax, int nmax, double *eig)
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

/* (static)
 * _zEigSymBisecR
 * - compute a unitary matrix consisting of eigenvectors
 *   for real-number eigenvalues.
 */
bool _zEigSymBisecR(zMat m, zVec eig, zMat r, int iter)
{
  register int i;
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
    zMatShift( ms, -zVecElem(eig,i)-shift );
    while( !zMatInv( ms, b ) )
      zMatShift( ms, -( shift*=10 ) );
    zEigPower( b, eigv, iter );
    zMatSetCol( r, i, eigv );
  }

 TERMINATE:
  zMatFree( ms );
  zMatFree( b );
  zVecFree( eigv );
  return ret;
}

/* zEigSymBisec
 * - compute eigenvalues and eigenvectors based on bisection method
 *   (by J. W. Givens 1954).
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

static void _zEigSymJacobiRot(zMat m, zMat r, int i, int j);
static double _zEigSymJacobiShift(zMat m, double *shift);

/* (static)
 * _zEigSymJacobiShift
 * - shift diagonal values of a matrix to accelerate diagonalization.
 */
double _zEigSymJacobiShift(zMat m, double *shift)
{
  zComplex c_shift[2];
  double b, c;

  b = - zMatElem(m,0,0) - zMatElem(m,1,1);
  c = zMatElem(m,0,0) * zMatElem(m,1,1)
      - zMatElem(m,1,0) * zMatElem(m,0,1);
  zQESolve( 1, b, c, c_shift );
  c = fabs( c_shift[0].re - zMatElem(m,0,0) )
    < fabs( c_shift[1].re - zMatElem(m,0,0) ) ?
    c_shift[0].re : c_shift[1].re;
  zMatShift( m, -c );
  return shift ? ( *shift = c ) : c;
}

/* (static)
 * _zEigSymJacobiRot
 * - transform a matrix by Jacobi's rotation based on Wilkinson's formula.
 */
void _zEigSymJacobiRot(zMat m, zMat r, int i, int j)
{
  register int k;
  double as, ad, ti, c, s;
  double tmp1, tmp2;

  as = 0.5 * ( zMatElem(m,i,i) + zMatElem(m,j,j) );
  ad = 0.5 * ( zMatElem(m,i,i) - zMatElem(m,j,j) );
  ti = sqrt( ad*ad + zSqr( zMatElem(m,i,j) ) );
  if( ad < 0 ) ti = -ti;
  c = sqrt( 0.5 + 0.5*ad/ti );
  s = 0.5 * zMatElem(m,i,j) / ( ti * c );

  zMatSetElem( m, i, i, as + ti );
  zMatSetElem( m, j, j, as - ti );
  zMatSetElem( m, i, j, 0 );
  zMatSetElem( m, j, i, 0 );
  for( k=0; k<zMatRowSizeNC(m); k++ ){
    /* update of transformation matrix */
    tmp1 = zMatElem( r, k, i );
    tmp2 = zMatElem( r, k, j );
    zMatSetElem( r, k, i, c * tmp1 + s * tmp2 );
    zMatSetElem( r, k, j,-s * tmp1 + c * tmp2 );
    /* update of eigenmatrix */
    if( k == i || k == j ) continue;
    tmp1 = zMatElem( m, i, k );
    tmp2 = zMatElem( m, j, k );
    zMatSetElem( m, i, k, c * tmp1 + s * tmp2 );
    zMatSetElem( m, j, k,-s * tmp1 + c * tmp2 );
    zMatSetElem( m, k, i, zMatElem( m, i, k ) );
    zMatSetElem( m, k, j, zMatElem( m, j, k ) );
  }
}

/* zEigSymJacobi
 * - diagonalize a symmetric matrix by Jacobi's method.
 *   when 'm' is not symmetric, it doesn't work validly.
 */
zVec zEigSymJacobi(zMat m, zVec eig, zMat r)
{
  register int i, j;
  int n = 0;
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
        if( zIsTiny( zMatElem(d,i,j) ) ) continue;
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
    zVecSetElem( eig, i, zMatElem(d,i,i) + shift );
  zMatFree( d );
  return eig;
}

/* singular value decomposition */

static int _zSVDSort(zVec sv, zMat u);

/* (static)
 * _zSVDSort
 * - sort singular values and corresponding bases.
 */
int _zSVDSort(zVec sv, zMat u)
{
  register int i;
  int im;

  for( i=0; i<zVecSizeNC(sv); i++ ){
    if( zIsTiny( zDataMax( zVecBuf(sv)+i, zVecSizeNC(sv)-i, &im ) ) )
      break;
    im += i;
    zVecSwapNC( sv, i, im );
    zMatSwapColNC( u, i, im );
  }
  return i;
}

/* zSVD
 * - singular value decomposition.
 */
int zSVD(zMat m, zVec sv, zMat u, zMat v)
{
  register int i, j;
  zMat c, w = NULL;
  double s;
  int rank;

  if( !( c = zMatAllocSqr( zMatRowSizeNC(m) ) ) ){
    ZALLOCERROR();
    rank = -1;
    goto TERMINATE;
  }

  zMulMatMatTNC( m, m, c );
  zEigSymJacobi( c, sv, u );
  for( i=0; i<zVecSizeNC(sv); i++ )
    zVecSetElem( sv, i, ( zVecElem(sv,i) < zTOL ) ? 0 : sqrt(zVecElem(sv,i)) );
  rank = _zSVDSort( sv, u );
  if( !( w = zMatAlloc( zMatRowSizeNC(u), rank ) ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  zMatGetNC( u, 0, 0, w );
  zMatSetRowSize( v, rank );

  zMulMatTMatNC( w, m, v );
  for( i=0; i<rank; i++ ){
    s = 1.0 / zVecElem(sv,i);
    for( j=0; j<zMatColSizeNC(m); j++ )
      zMatElem(v,i,j) *= s;
  }
 TERMINATE:
  zMatFree( w );
  zMatFree( c );
  return rank;
}

static bool _zSVMat(zMat m, zMat *r, zVec *e);

/* (static)
 * _zSVMat
 * - allocate internal matrix and vector for workspace.
 */
bool _zSVMat(zMat m, zMat *r, zVec *e)
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

/* zSVMax
 * - maximum singular value of a matrix.
 */
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

/* zSVMin
 * - minimum singular value of a matrix.
 */
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

/* zMatCondNum
 * - condition number of a matrix.
 */
double zMatCondNum(zMat m)
{
  double smax, smin;

  smax = zSVMax( m );
  smin = zSVMin( m );
  if( smin/smin < zTOL ) return HUGE_VAL;
  return smax / smin;
}
