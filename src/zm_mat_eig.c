/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mat_eig - eigensystem of matrices.
 */

#include <zm/zm_mat_eig.h>

/* Householder transformation (destructive). */
static void _zMatToHessenbergHouseholder(zMat m, zMat p, int from, int to, zVec u, zVec v, zVec w)
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
  /* transformation matrix */
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
static zVec _zMatToHessenbergHouseholderVec(const zMat m, int col, int from, int to, zVec u)
{
  int i;
  double s;

  for( s=0, i=from; i<to; i++ )
    s += zSqr( zVecSetElemNC( u, i, zMatElemNC(m,i,col) ) );
  if( zIsTiny( s ) ) return u;
  zVecElemNC(u,from) -= ( s = -zSgn(zMatElemNC(m,from,col)) * sqrt(s) );
  zRawVecDivDRC( &zVecElemNC(u,from), sqrt( 2*(s-zMatElemNC(m,from,col))*s ), to-from );
  return u;
}

/* Hessenberg matrix using Householder transformation. */
zMat zMatToHessenberg(const zMat m, zMat h, zMat p)
{
  int n, n2, i;
  zVec u, v, w;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  /* target hessian matrix */
  if( h ){
    if( !zMatSizeEqual( m, h ) ){
      ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
      return NULL;
    }
    zMatCopyNC( m, h );
  } else
    h = m; /* if h is the null pointer, m is directly transformed. */
  /* transformation matrix */
  if( p && !zMatSizeEqual( h, p ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
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
    _zMatToHessenbergHouseholderVec( h, n, n+1, zVecSizeNC(u), u );
    _zMatToHessenbergHouseholder( h, p, n+1, zVecSizeNC(u), u, v, w );
    for( i=n+2; i<zMatRowSizeNC(h); i++ ) /* need refinement? */
      if( !zIsTiny( zMatElemNC(h,i,n) ) ) goto RETRY;
  }
 TERMINATE:
  zVecFreeAtOnce( 3, u, v, w );
  return h;
}

/* eigenvalue analysis by double QR method and inverse iteration */

/* Householder's method to transform a Hesseian matrix to a block uppertriangular matrix. */
static void _zMatEigDQRHouseholder(zMat a, int r, int c, double g1, double g2, double g3)
{
  double s, den, alpha;
  int i, c1, c2, n;

  s = ( s = g1*g1 + g2*g2 + g3*g3 ) < 0 ? 0 : sqrt( s );
  den = g1 + ( g1 > 0 ? s : -s );
  g2 /= den;
  g3 /= den;
  alpha = 2.0 / ( 1.0 + g2*g2 + g3*g3 );

  c1 = c + 1;
  c2 = c + 2;
  for( i=zMax(c-1,0); i<=r; i++ ){
    s = ( zMatElemNC(a,c,i) + g2*zMatElemNC(a,c1,i) + g3*zMatElemNC(a,c2,i) ) * alpha;
    zMatElemNC(a,c ,i) -= s;
    zMatElemNC(a,c1,i) -= s * g2;
    zMatElemNC(a,c2,i) -= s * g3;
  }
  n = zLimit( c+3, r, zMatRowSizeNC(a)-1 );
  for( i=0; i<=n; i++ ){
    s = ( zMatElemNC(a,i,c) + g2*zMatElemNC(a,i,c1) + g3*zMatElemNC(a,i,c2) ) * alpha;
    zMatElemNC(a,i,c ) -= s;
    zMatElemNC(a,i,c1) -= s * g2;
    zMatElemNC(a,i,c2) -= s * g3;
  }
}

/* double QR method. */
bool zMatEigDQR(const zMat m, zCVec eigval, int iter)
{
  int  i, j, r, r1;
  double a0, a1, b, c;
  zMat a;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return false;
  }
  if( zMatRowSizeNC(m) != zCVecSizeNC(eigval) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_CVEC );
    return false;
  }
  if( !( a = zMatClone( m ) ) ){
    ZALLOCERROR();
    return false;
  }

  zMatToHessenberg( a, NULL, NULL );
  for( r=zMatRowSizeNC(a)-1; r>=0; ){
    if( r == 0 ){
      zComplexCreate( zCVecElemNC(eigval,0), zMatBuf(a)[0], 0.0 );
      break;
    }
    ZITERINIT( iter );
    for( i=0; i<iter; i++ ){
      r1 = r - 1;
      if( zIsTiny( zMatElemNC(a,r,r1) ) ){
        zComplexCreate( zCVecElemNC(eigval,r), zMatElemNC(a,r,r), 0.0 );
        r = r1;
        break;
      }
      b = -zMatElemNC(a,r1,r1) - zMatElemNC(a,r,r);
      c =  zMatElemNC(a,r1,r1)*zMatElemNC(a,r,r) - zMatElemNC(a,r1,r)*zMatElemNC(a,r,r1);
      if( r1 == 0 || zIsTiny( zMatElemNC(a,r1,r1-1) ) ){
        zQESolve( 1.0, b, c, zCVecElemNC(eigval,r1) );
        r -= 2;
        break;
      }
      b += ( a0 = zMatBuf(a)[0] );
      a1 = zMatElemNC(a,1,0);
      _zMatEigDQRHouseholder( a, r, 0,
        zMatElemNC(a,0,1)*a1+a0*b+c, (zMatElemNC(a,1,1)+b)*a1, zMatElemNC(a,2,1)*a1 );
      for( j=1; j<=r1-1; j++ )
        _zMatEigDQRHouseholder( a, r, j,
          zMatElemNC(a,j,j-1), zMatElemNC(a,j+1,j-1), zMatElemNC(a,j+2,j-1) );
      _zMatEigDQRHouseholder( a, r, j,
        zMatElemNC(a,r1,r1-1), zMatElemNC(a,r ,r1-1), 0 );
    }
  }
  zMatFree( a );
  return true;
}

/* calculate the dominant eigenvalue. */
double zMatEigPower(const zMat m, zVec eigvec, int iter)
{
  zVec ev, err;
  double eigval = 1.0;
  int i;

  err = zVecAlloc( zVecSizeNC(eigvec) );
  ev = zVecAlloc( zVecSizeNC(eigvec) );
  if( !err || !ev ) goto TERMINATE;
  zVecRandUniform( eigvec, -1, 1 );

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecNormalize( eigvec, ev );
    zMulMatVecNC( m, ev, eigvec );
    eigval = zVecInnerProd( ev, eigvec );
    zVecCatNC( eigvec, -eigval, ev, err );
    if( zVecIsTiny( err ) ) goto TERMINATE;
  }
 TERMINATE:
  zVecDivDRC( eigvec, eigval );
  zVecFree( ev );
  zVecFree( err );
  return eigval;
}

/* calculate the minimal eigenvalue. */
double zMatEigPowerInv(const zMat m, zVec eigvec, int iter)
{
  zMat m_inv;
  double eigval = 0;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return 0;
  }
  if( !( m_inv = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) ) ) )
    return 0;
  if( zMatInv( m, m_inv ) )
    eigval = zMatEigPower( m_inv, eigvec, iter );
  else
    ZRUNERROR( ZM_ERR_MAT_SINGULAR );

  zMatFree( m_inv );
  return 1.0 / eigval;
}

/* compute an eigenvector for a real-number eigenvalue. */
static bool _zMatEigVecReal(const zMat m, double eigval, zCVec eigvec, int iter)
{
  zMat ms, b;
  zVec eigvec_real;
  double shift = zTOL * 1.0e1;
  bool ret = true;

  ms = zMatClone( m );
  b = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  eigvec_real = zVecAlloc( zCVecSizeNC(eigvec) );
  if( !ms || !b || !eigvec_real ){
    ZALLOCERROR();
    ret = false;
    goto TERMINATE;
  }

  zVecRandUniform( eigvec_real, -1, 1 );
  zMatShift( ms, -eigval - shift );
  while( !zMatInv( ms, b ) )
    zMatShift( ms, -( shift*=10 ) );
  zMatEigPower( b, eigvec_real, iter );
  zVecToCVec( eigvec_real, eigvec );

 TERMINATE:
  zMatFree( ms );
  zMatFree( b );
  zVecFree( eigvec_real );
  return ret;
}

/* compute eigenvectors for complex-number eigenvalues. */
static bool _zMatEigVecComplex(const zMat m, zComplex *eigval, zCVec eigvec1, zCVec eigvec2, int iter)
{
  zMat ms, b;
  zVec eigvec;
  double shift = zTOL * 1.0e1, s;
  int i, n;
  bool ret = true;

  n = 2 * zMatRowSizeNC(m);
  ms = zMatAllocSqr( n );
  b = zMatAllocSqr( n );
  eigvec = zVecAlloc( n );
  if( !ms || !b || !eigvec ){
    ret = false;
    goto TERMINATE;
  }

  zMatPutNC( ms, 0, 0, m );
  zMatPutNC( ms, zMatRowSizeNC(m), zMatColSizeNC(m), m );
  s = eigval->re + shift;
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    zMatElemNC(ms,i,i) -= s;
    zMatSetElemNC( ms, i+zMatRowSizeNC(m), i,-eigval->im );
    zMatSetElemNC( ms, i, i+zMatColSizeNC(m), eigval->im );
    zMatElemNC(ms,i+zMatRowSizeNC(m),i+zMatColSizeNC(m)) -= s;
  }
  zVecRandUniform( eigvec, -1, 1 );
  while( !zMatInv( ms, b ) )
    zMatShift( ms, -( shift *= 10 ) );
  zMatEigPower( b, eigvec, iter );

  for( i=0; i<zCVecSizeNC(eigvec1); i++ ){
    zComplexCreate( zCVecElem(eigvec1,i),
      zVecElemNC(eigvec,i), zVecElemNC(eigvec,i+zCVecSizeNC(eigvec1)) );
    zComplexCreate( zCVecElem(eigvec2,i),
      zVecElemNC(eigvec,i),-zVecElemNC(eigvec,i+zCVecSizeNC(eigvec1)) );
  }
  zCVecNormalizeDRC( eigvec1 );
  zCVecNormalizeDRC( eigvec2 );

 TERMINATE:
  zMatFree( ms );
  zMatFree( b );
  zVecFree( eigvec );
  return ret;
}

/* compute sets of eigenvalues and eigenvectors of a square matrix. */
#define Z_EIGSYS_ITER_NUM 10
bool zMatEig(const zMat m, zCVec eigval, zCMat eigbase, int iter)
{
  zCVec eigvec1, eigvec2;
  bool ret = false;
  int i;

  eigvec1 = zCVecAlloc( zCMatRowSizeNC(eigbase) );
  eigvec2 = zCVecAlloc( zCMatRowSizeNC(eigbase) );
  if( !eigvec1 || !eigvec2 ) goto TERMINATE;
  /* eigenvalues by double QR method */
  if( !zMatEigDQR( m, eigval, iter ) ) goto TERMINATE;
  /* eigenvectors by inverse iteration method */
  for( i=0; i<zMatRowSizeNC(m); ){
    if( zComplexIsReal( zCVecElemNC(eigval,i), zTOL ) ){
      if( !_zMatEigVecReal( m, zCVecElemNC(eigval,i)->re, eigvec1, Z_EIGSYS_ITER_NUM ) )
        break;
      zCMatPutColNC( eigbase, i, eigvec1 );
      i++;
    } else{
      if( !_zMatEigVecComplex( m, zCVecElemNC(eigval,i), eigvec1, eigvec2, Z_EIGSYS_ITER_NUM ) )
        break;
      zCMatPutCol( eigbase, i,   eigvec1 );
      zCMatPutCol( eigbase, i+1, eigvec2 );
      i += 2;
    }
  }
  ret = true;
  /* return the number of eigenvalues and eigenvectors */
 TERMINATE:
  zCVecFree( eigvec1 );
  zCVecFree( eigvec2 );
  return ret;
}

/* check if sizes of matrices and a vector for eigenvalues are consisten. */
#define _zMatSymEigCheckSize(m,eigval,eigbase) do{\
  if( !zMatIsSqr( m ) ){\
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );\
    return false;\
  }\
  if( !eigbase || !zMatSizeEqual( m, eigbase ) ){\
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );\
    return false;\
  }\
  if( !eigval || !zMatRowVecSizeEqual( m, eigval ) ){\
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );\
    return false;\
  }\
} while(0)

/* compute the range in which all eigenvalues exist based on Gerschgorin's theorem. */
static void _zMatSymEigBisecRange(const zMat m, double *emin, double *emax)
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
static int _zMatSymEigBisecSturm(const zMat a, double e)
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
static void _zMatSymEigBisecOne(const zMat a, double emin, int nmin, double emax, int nmax, double *eig)
{
  double e, d, d_old;
  int n;

  for( d_old=HUGE_VAL; ; d_old=d ){
    e = emin + 0.5 * ( d = emax - emin );
    if( d == d_old ) break;
    n = _zMatSymEigBisecSturm( a, e );
    if( n == nmin ) emin = e;
    else emax = e;
  }
  *eig = e;
}

/* compute all eigenvalues recursively by bisection method. */
static void _zMatSymEigBisecRec(const zMat a, double emin, int nmin, double emax, int nmax, double *eig)
{
  int n;
  double e;

  if( nmin - nmax == 1 )
    return _zMatSymEigBisecOne( a, emin, nmin, emax, nmax, eig );

  n = _zMatSymEigBisecSturm( a, ( e = 0.5 * ( emax + emin ) ) );
  if( n > nmax )
    _zMatSymEigBisecRec( a, e, n, emax, nmax, &eig[0] );
  if( n < nmin )
    _zMatSymEigBisecRec( a, emin, nmin, e, n, &eig[n-nmax] );
}

/* compute a unitary matrix consisting of eigenvectors for real-number eigenvalues. */
static bool _zMatSymEigBisecR(const zMat m, zVec eigval, zMat eigbase, int iter)
{
  int i;
  zMat ms, b;
  zVec eigvec;
  double shift = zTOL * 1.0e1;
  bool ret = true;

  ms = zMatAllocSqr( zMatRowSizeNC(m) );
  b = zMatAllocSqr( zMatRowSizeNC(m) );
  eigvec = zVecAlloc( zMatRowSizeNC(m) );
  if( !ms || !b || !eigvec ){
    ZALLOCERROR();
    ret = false;
    goto TERMINATE;
  }

  for( i=0; i<zVecSizeNC(eigval); i++ ){
    zMatCopyNC( m, ms );
    zVecRandUniform( eigvec, -1, 1 );
    zMatShift( ms, -zVecElemNC(eigval,i)-shift );
    while( !zMatInv( ms, b ) )
      zMatShift( ms, -( shift*=10 ) );
    zMatEigPower( b, eigvec, iter );
    zMatPutCol( eigbase, i, eigvec );
  }

 TERMINATE:
  zMatFree( ms );
  zMatFree( b );
  zVecFree( eigvec );
  return ret;
}

/* diagonalize a symmetric matrix by bisection method (J. W. Givens, 1954). */
bool zMatSymEigBisec(const zMat m, zVec eigval, zMat eigbase)
{
  double emin, emax;
  int nmin, nmax;

  _zMatSymEigCheckSize( m, eigval, eigbase );
  /* eigenvalues */
  zMatCopyNC( m, eigbase ); /* temporary working space */
  zMatToHessenberg( eigbase, NULL, NULL );
  _zMatSymEigBisecRange( eigbase, &emin, &emax );
  nmin = _zMatSymEigBisecSturm( eigbase, emin );
  nmax = _zMatSymEigBisecSturm( eigbase, emax );
  _zMatSymEigBisecRec( eigbase, emin, nmin, emax, nmax, zVecBuf(eigval) );
  _zMatSymEigBisecR( m, eigval, eigbase, Z_EIGSYS_ITER_NUM );
  return true;
}

/* shift diagonal values of a matrix to accelerate diagonalization. */
static double _zMatSymEigJacobiShift(zMat m, double *shift)
{
  zComplex c_shift[2];
  double b, c;

  b =-zMatBuf(m)[0] - zMatElemNC(m,1,1);
  c = zMatBuf(m)[0] * zMatElemNC(m,1,1) - zMatElemNC(m,1,0) * zMatElemNC(m,0,1);
  zQESolve( 1, b, c, c_shift );
  c = fabs( c_shift[0].re - zMatBuf(m)[0] ) < fabs( c_shift[1].re - zMatBuf(m)[0] ) ?
    c_shift[0].re : c_shift[1].re;
  zMatShift( m, -c );
  return shift ? ( *shift = c ) : c;
}

/* transform a matrix by Jacobi's rotation based on Wilkinson's formula. */
static void _zMatSymEigJacobiRot(zMat m, zMat eigbase, int i, int j)
{
  int k;
  double as, ad, ti, c, s;
  double tmp1, tmp2;

  as = 0.5 * ( zMatElemNC(m,i,i) + zMatElemNC(m,j,j) );
  ad = 0.5 * ( zMatElemNC(m,i,i) - zMatElemNC(m,j,j) );
  ti = sqrt( ad*ad + _zSqr( zMatElemNC(m,i,j) ) );
  if( ad < 0 ) ti = -ti;
  c = sqrt( 0.5 + 0.5*ad/ti );
  s = 0.5 * zMatElemNC(m,i,j) / ( ti * c );

  zMatSetElemNC( m, i, i, as + ti );
  zMatSetElemNC( m, j, j, as - ti );
  zMatSetElemNC( m, i, j, 0 );
  zMatSetElemNC( m, j, i, 0 );
  for( k=0; k<zMatRowSizeNC(m); k++ ){
    /* update of transformation matrix */
    tmp1 = zMatElemNC( eigbase, k, i );
    tmp2 = zMatElemNC( eigbase, k, j );
    zMatSetElemNC( eigbase, k, i, c * tmp1 + s * tmp2 );
    zMatSetElemNC( eigbase, k, j,-s * tmp1 + c * tmp2 );
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

/* diagonalize a symmetric matrix by Jacobi's method. */
bool zMatSymEigJacobi(const zMat m, zVec eigval, zMat eigbase)
{
  int i, j, n = 0;
  double shift;
  bool is_complete;
  zMat d;

  _zMatSymEigCheckSize( m, eigval, eigbase );
  if( !( d = zMatClone( m ) ) ){ /* diagonal matrix */
    ZALLOCERROR();
    return false;
  }
  zMatIdent( eigbase );
  _zMatSymEigJacobiShift( d, &shift );
  do{
    is_complete = true;
    for( i=1; i<zMatRowSizeNC(d); i++ )
      for( j=0; j<i; j++ ){
        if( zIsTiny( zMatElemNC(d,i,j) ) ) continue;
        is_complete = false;
        /* iterative elimination of non-diagonal components */
        _zMatSymEigJacobiRot( d, eigbase, i, j );
      }
    if( n++ > Z_MAX_ITER_NUM ){
      ZITERWARN( Z_MAX_ITER_NUM );
      goto TERMINATE;
    }
  } while( !is_complete );

 TERMINATE:
  for( j=0; j<zMatRowSizeNC(m); j++ )
    zVecSetElemNC( eigval, j, zMatElemNC(d,j,j) + shift );
  zMatFree( d );
  return true;
}

/* sort singular values and corresponding bases. */
static int _zMatSVDSort(const zMat m, zVec sv, zMat u)
{
  int i, im, rank;

  rank = zMin( zMatRowSizeNC(m), zMatColSizeNC(m) );
  for( i=0; i<rank; i++ ){
    if( zIsTiny( zDataMax( zVecBuf(sv)+i, zVecSizeNC(sv)-i, &im ) ) ) break;
    im += i;
    zVecSwapNC( sv, i, im );
    zVecSetElemNC( sv, i, sqrt( zVecElemNC(sv,i) ) );
    zMatSwapColNC( u, i, im );
  }
  if( ( rank = i ) < zVecSizeNC(sv) )
    for( i=rank; i<zVecSizeNC(sv); i++ ) /* replace theoretically zero singular values with zeroes */
      zVecSetElemNC( sv, i, 0 );
  return rank;
}

/* singular value decomposition. */
int zMatSVD(const zMat m, zMat u, zVec sv, zMat v)
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
  zMatSymEigJacobi( c, sv, u );
  rank = _zMatSVDSort( m, sv, u );
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
static bool _zMatSingularValueMat(const zMat m, zMat *r, zVec *e)
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
double zMatSingularValueMax(const zMat m)
{
  zMat r;
  zVec e;
  double s = 0;

  if( _zMatSingularValueMat( m, &r, &e ) )
    s = sqrt( zMatEigPower( r, e, 0 ) );
  zMatFree( r );
  zVecFree( e );
  return s;
}

/* minimum singular value of a matrix. */
double zMatSingularValueMin(const zMat m)
{
  zMat r;
  zVec e;
  double s = 0;

  if( _zMatSingularValueMat( m, &r, &e ) )
    s = sqrt( zMatEigPowerInv( r, e, 0 ) );
  zMatFree( r );
  zVecFree( e );
  return s;
}

/* condition number of a matrix. */
double zMatCondNum(const zMat m)
{
  zMat r;
  zVec e;
  double smax, smin;

  if( !_zMatSingularValueMat( m, &r, &e ) ) return NAN;
  smax = sqrt( zMatEigPower( r, e, 0 ) );
  smin = sqrt( zMatEigPowerInv( r, e, 0 ) );
  zMatFree( r );
  zVecFree( e );
  return smin > zTOL ? smax / smin : HUGE_VAL;
}
