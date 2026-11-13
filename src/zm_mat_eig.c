/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mat_eig - eigensystem of matrices.
 */

#include <zm/zm_mat_eig.h>

/* Householder reflection */

/* unit normal vector nu of the Householder reflection I - 2 nu nu^T that maps x to s e1, where |s|=|x|.
 * The sign of s is chosen to be opposite to that of the first component of x in order to avoid
 * cancellation. It returns the false value if x is the zero vector, for which no reflection is needed. */
static bool _zHouseholderVec(const double *x, int n, double *nu, double *s)
{
  double scale = 0, r = 0, d;
  int i;

  for( i=0; i<n; i++ )
    if( fabs( x[i] ) > scale ) scale = fabs( x[i] );
  if( scale == 0 ) return false;
  for( i=0; i<n; i++ ){ /* scaled to avoid overflow and underflow */
    nu[i] = x[i] / scale;
    r += nu[i] * nu[i];
  }
  r = sqrt( r );
  *s = x[0] >= 0 ? -r : r;
  nu[0] -= *s; /* |nu|^2 = 2 r ( r + |x0| ) */
  d = sqrt( 2 * r * ( r + fabs( x[0] / scale ) ) );
  for( i=0; i<n; i++ ) nu[i] /= d;
  *s *= scale;
  return true;
}

/* apply a Householder reflection to rows from r0 to r0+n-1 of a matrix from the left side,
 * where only columns from c0 to c1 are modified. w is a workspace with the column size of m. */
static void _zHouseholderApplyL(zMat m, const double *nu, int n, int r0, int c0, int c1, double *w)
{
  int i, j;
  double a, *row;

  for( j=c0; j<=c1; j++ ) w[j] = 0;
  for( i=0; i<n; i++ ){ /* w = nu^T m */
    row = zMatRowBufNC( m, r0+i );
    for( a=nu[i], j=c0; j<=c1; j++ ) w[j] += a * row[j];
  }
  for( i=0; i<n; i++ ){ /* m -= 2 nu w */
    row = zMatRowBufNC( m, r0+i );
    for( a=2*nu[i], j=c0; j<=c1; j++ ) row[j] -= a * w[j];
  }
}

/* apply a Householder reflection to columns from c0 to c0+n-1 of a matrix from the right side,
 * where only rows from r0 to r1 are modified. */
static void _zHouseholderApplyR(zMat m, const double *nu, int n, int c0, int r0, int r1)
{
  int i;
  double *p;

  for( i=r0; i<=r1; i++ ){
    p = &zMatElemNC(m,i,c0);
    zRawVecCatDRC( p, -2*zRawVecInnerProd( p, nu, n ), nu, n );
  }
}

/* Hessenberg matrix using Householder transformation. */
zMat zMatToHessenberg(const zMat m, zMat h, zMat p)
{
  int n, k, i, len;
  double s;
  zVec nu, w;

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
  if( p ) zMatIdentNC( p );
  n = zMatRowSizeNC(h);
  nu = zVecAlloc( n );
  w = zVecAlloc( n );
  if( !nu || !w ){
    ZALLOCERROR();
    h = NULL;
    goto TERMINATE;
  }
  for( k=0; k<n-2; k++ ){
    len = n - k - 1;
    for( i=0; i<len; i++ ) /* the k-th column below the subdiagonal component */
      zVecSetElemNC( w, i, zMatElemNC(h,k+1+i,k) );
    /* a reflection is skipped only if the column is exactly zero. */
    if( !_zHouseholderVec( zVecBuf(w), len, zVecBuf(nu), &s ) ) continue;
    _zHouseholderApplyL( h, zVecBuf(nu), len, k+1, k+1, n-1, zVecBuf(w) );
    _zHouseholderApplyR( h, zVecBuf(nu), len, k+1, 0, n-1 );
    /* the k-th column is theoretically mapped to s e1. */
    zMatSetElemNC( h, k+1, k, s );
    for( i=k+2; i<n; i++ ) zMatSetElemNC( h, i, k, 0 );
    if( p ) _zHouseholderApplyR( p, zVecBuf(nu), len, k+1, 0, n-1 );
  }
 TERMINATE:
  zVecFree( nu );
  zVecFree( w );
  return h;
}

/* eigenvalue analysis by double-shift QR method (J. G. F. Francis, 1961) */

/* find the top row of the lowest unreduced block of a Hessenberg matrix that ends at the r-th row.
 * A subdiagonal component is regarded as zero if it is negligible relative to its neighboring
 * diagonal components, and then zero is exactly assigned to it.
 * If the iteration stagnates, a subdiagonal component not greater than n eps times the norm of the
 * matrix is also regarded as zero, which is a perturbation comparable to the backward error of the
 * QR method itself. This happens for a cluster of close eigenvalues, for which the subdiagonal
 * components stay at the level of rounding errors. */
static int _zMatEigDQRSplit(zMat h, int r, double norm, bool stagnant)
{
  int l;
  double s, a;

  for( l=r; l>0; l-- ){
    if( ( s = fabs( zMatElemNC(h,l-1,l-1) ) + fabs( zMatElemNC(h,l,l) ) ) == 0 ) s = norm;
    if( ( a = fabs( zMatElemNC(h,l,l-1) ) ) <= DBL_EPSILON * s || ( stagnant && a <= zMatRowSizeNC(h) * DBL_EPSILON * norm ) ){
      zMatSetElemNC( h, l, l-1, 0 );
      break;
    }
  }
  return l;
}

/* one step of implicit double-shift QR iteration (bulge chasing) on an unreduced block from the l-th
 * to the r-th rows, where s and t are the sum and the product of the shifts, respectively.
 * If z is not the null pointer, the whole Hessenberg matrix is transformed and the transformation
 * is accumulated to z. Otherwise, only the active block is transformed. */
static void _zMatEigDQRStep(zMat h, zMat z, int l, int r, double s, double t, double *nu, double *w)
{
  int n, k, nr, cmax, rmin;
  double x[3], sv;

  n = zMatRowSizeNC(h);
  cmax = z ? n - 1 : r;
  rmin = z ? 0 : l;
  /* the first column of (H-mu1 I)(H-mu2 I) = H^2 - s H + t I */
  x[0] = zMatElemNC(h,l,l) * ( zMatElemNC(h,l,l) - s ) + zMatElemNC(h,l,l+1) * zMatElemNC(h,l+1,l) + t;
  x[1] = zMatElemNC(h,l+1,l) * ( zMatElemNC(h,l,l) + zMatElemNC(h,l+1,l+1) - s );
  x[2] = zMatElemNC(h,l+1,l) * zMatElemNC(h,l+2,l+1);
  for( k=l; k<r; k++ ){
    nr = _zMin( 3, r - k + 1 );
    if( k > l ){ /* the bulge to be chased */
      x[0] = zMatElemNC(h,k,k-1);
      x[1] = zMatElemNC(h,k+1,k-1);
      x[2] = nr == 3 ? zMatElemNC(h,k+2,k-1) : 0;
    }
    if( !_zHouseholderVec( x, nr, nu, &sv ) ) continue;
    _zHouseholderApplyL( h, nu, nr, k, k, cmax, w );
    if( k > l ){
      zMatSetElemNC( h, k, k-1, sv );
      zMatSetElemNC( h, k+1, k-1, 0 );
      if( nr == 3 ) zMatSetElemNC( h, k+2, k-1, 0 );
    }
    _zHouseholderApplyR( h, nu, nr, k, rmin, _zMin( k+3, r ) );
    if( z ) _zHouseholderApplyR( z, nu, nr, k, 0, n-1 );
  }
}

/* rotate the p-th and (p+1)-th rows and columns of a quasi-upper triangular matrix t by a Givens
 * rotation G = [ c -s ; s c ] as t <- G^T t G, and accumulate it to z as z <- z G. */
static void _zMatEigDQRRot(zMat t, zMat z, int p, double c, double s)
{
  int i, q = p + 1;
  double a, b;

  for( i=p; i<zMatColSizeNC(t); i++ ){
    a = zMatElemNC(t,p,i);
    b = zMatElemNC(t,q,i);
    zMatSetElemNC( t, p, i, c*a + s*b );
    zMatSetElemNC( t, q, i,-s*a + c*b );
  }
  for( i=0; i<=q; i++ ){
    a = zMatElemNC(t,i,p);
    b = zMatElemNC(t,i,q);
    zMatSetElemNC( t, i, p, c*a + s*b );
    zMatSetElemNC( t, i, q,-s*a + c*b );
  }
  for( i=0; i<zMatRowSizeNC(z); i++ ){
    a = zMatElemNC(z,i,p);
    b = zMatElemNC(z,i,q);
    zMatSetElemNC( z, i, p, c*a + s*b );
    zMatSetElemNC( z, i, q,-s*a + c*b );
  }
}

/* eigenvalues of a deflated 2x2 block at the p-th and (p+1)-th rows.
 * If z is not the null pointer and the eigenvalues are real, the block is triangularized by a Givens
 * rotation, which is accumulated to z. */
static void _zMatEigDQR2x2(zMat h, zMat z, int p, zComplex *e1, zComplex *e2)
{
  double a, b, c, d, hd, disc, q, x, l1, l2, v1, v2, r1, r2;

  a = zMatElemNC(h,p,p);   b = zMatElemNC(h,p,p+1);
  c = zMatElemNC(h,p+1,p); d = zMatElemNC(h,p+1,p+1);
  if( c == 0 ){
    _zComplexCreate( e1, a, 0 );
    _zComplexCreate( e2, d, 0 );
    return;
  }
  hd = 0.5 * ( a - d );
  if( ( disc = hd*hd + b*c ) < 0 ){ /* a pair of complex conjugate eigenvalues */
    _zComplexCreate( e1, d + hd, sqrt( -disc ) );
    _zComplexCreate( e2, d + hd,-sqrt( -disc ) );
    return;
  }
  /* real eigenvalues computed without cancellation */
  q = sqrt( disc );
  if( ( x = hd + ( hd >= 0 ? q : -q ) ) == 0 )
    l1 = l2 = d + hd;
  else{
    l1 = d + x;
    l2 = d - ( b / x ) * c;
  }
  if( z ){ /* triangularize the block by the eigenvector for l1 */
    r1 = _zSqr( b ) + _zSqr( l1 - a );
    r2 = _zSqr( l1 - d ) + _zSqr( c );
    if( r1 >= r2 ){
      v1 = b;      v2 = l1 - a;
    } else{
      v1 = l1 - d; v2 = c;
      r1 = r2;
    }
    if( r1 > 0 ){
      r1 = sqrt( r1 );
      _zMatEigDQRRot( h, z, p, v1/r1, v2/r1 );
      zMatSetElemNC( h, p+1, p, 0 );
      l1 = zMatElemNC(h,p,p);
      l2 = zMatElemNC(h,p+1,p+1);
    }
  }
  _zComplexCreate( e1, l1, 0 );
  _zComplexCreate( e2, l2, 0 );
}

/* complex arithmetics for back substitution */
static void _zComplexCMulSub(zComplex *c, double a, const zComplex *y) /* c -= a y */
{
  c->re -= a * y->re;
  c->im -= a * y->im;
}

/* replace a tiny pivot with a perturbed one. */
static void _zMatEigDQRPerturb(zComplex *c, double delta)
{
  if( _zComplexAbs( c ) < delta ) _zComplexCreate( c, delta, 0 );
}

/* back substitution to find an eigenvector y of a quasi-upper triangular matrix t for an eigenvalue
 * lambda, where the components of y from the (i+1)-th to the last-th are given. */
static void _zMatEigDQRBackSubst(zMat t, const zComplex *lambda, zComplex *y, int i, int last, double delta)
{
  int j, k;
  zComplex r1, r2, a11, a12, a21, a22, m;
  double ymax;

  for( ; i>=0; i-- ){
    if( i > 0 && zMatElemNC(t,i,i-1) != 0 ){ /* 2x2 block of complex conjugate eigenvalues */
      _zComplexZero( &r1 );
      _zComplexZero( &r2 );
      for( j=i+1; j<=last; j++ ){
        _zComplexCMulSub( &r1, zMatElemNC(t,i-1,j), &y[j] );
        _zComplexCMulSub( &r2, zMatElemNC(t,i,j), &y[j] );
      }
      _zComplexCreate( &a11, zMatElemNC(t,i-1,i-1) - lambda->re, -lambda->im );
      _zComplexCreate( &a12, zMatElemNC(t,i-1,i), 0 );
      _zComplexCreate( &a21, zMatElemNC(t,i,i-1), 0 );
      _zComplexCreate( &a22, zMatElemNC(t,i,i) - lambda->re, -lambda->im );
      /* Gaussian elimination with partial pivoting */
      if( _zComplexSqrAbs( &a11 ) < _zComplexSqrAbs( &a21 ) ){
        zSwap( zComplex, a11, a21 );
        zSwap( zComplex, a12, a22 );
        zSwap( zComplex, r1, r2 );
      }
      _zMatEigDQRPerturb( &a11, delta );
      zComplexCDiv( &a21, &a11, &m );
      _zComplexCMul( &m, &a12, &a21 ); /* a21 is reused as a temporary */
      a22.re -= a21.re; a22.im -= a21.im;
      _zComplexCMul( &m, &r1, &a21 );
      r2.re -= a21.re; r2.im -= a21.im;
      _zMatEigDQRPerturb( &a22, delta );
      zComplexCDiv( &r2, &a22, &y[i] );
      _zComplexCMul( &a12, &y[i], &a21 );
      r1.re -= a21.re; r1.im -= a21.im;
      zComplexCDiv( &r1, &a11, &y[i-1] );
      i--;
    } else{
      _zComplexZero( &r1 );
      for( j=i+1; j<=last; j++ )
        _zComplexCMulSub( &r1, zMatElemNC(t,i,j), &y[j] );
      _zComplexCreate( &a11, zMatElemNC(t,i,i) - lambda->re, -lambda->im );
      _zMatEigDQRPerturb( &a11, delta );
      zComplexCDiv( &r1, &a11, &y[i] );
    }
    /* rescale to avoid overflow */
    if( ( ymax = _zComplexAbs( &y[i] ) ) > 1.0e100 )
      for( k=i; k<=last; k++ ){
        y[k].re /= ymax;
        y[k].im /= ymax;
      }
  }
}

/* transform an eigenvector y of a quasi-upper triangular matrix to that of the original matrix by z,
 * normalize it, and put it into the col-th column of eigbase. */
static void _zMatEigDQRPutVec(zMat z, zComplex *y, int last, zCMat eigbase, int col, bool conj)
{
  int i, j, imax = 0;
  double norm = 0, amax = -1, a;
  zComplex *x, f;

  for( i=0; i<zMatRowSizeNC(z); i++ ){
    x = zCMatElemNC(eigbase,i,col);
    _zComplexZero( x );
    for( j=0; j<=last; j++ ){
      x->re += zMatElemNC(z,i,j) * y[j].re;
      x->im += zMatElemNC(z,i,j) * y[j].im;
    }
    if( conj ) x->im = -x->im;
    norm += ( a = _zComplexSqrAbs( x ) );
    if( a > amax ){
      amax = a;
      imax = i;
    }
  }
  /* normalize the vector so that the absolutely maximum component is real and positive. */
  amax = sqrt( amax );
  norm = sqrt( norm );
  x = zCMatElemNC(eigbase,imax,col);
  _zComplexCreate( &f, x->re / ( amax * norm ), -x->im / ( amax * norm ) );
  for( i=0; i<zMatRowSizeNC(z); i++ ){
    x = zCMatElemNC(eigbase,i,col);
    a = x->re * f.re - x->im * f.im;
    x->im = x->re * f.im + x->im * f.re;
    x->re = a;
  }
  zCMatElemNC(eigbase,imax,col)->im = 0;
}

/* eigenvectors from a real Schur form t and Schur vectors z. */
static bool _zMatEigDQRVec(zMat t, zMat z, zCVec eigval, zCMat eigbase, double norm)
{
  int k, j, n;
  zComplex *y;
  double a, b, c, d, delta;
  zComplex *lambda;

  n = zMatRowSizeNC(t);
  if( !( y = zAlloc( zComplex, n ) ) ){
    ZALLOCERROR();
    return false;
  }
  delta = _zMax( DBL_EPSILON * norm, DBL_MIN );
  for( k=0; k<n; ){
    lambda = zCVecElemNC(eigval,k);
    for( j=0; j<n; j++ ) _zComplexZero( &y[j] );
    if( k < n-1 && zMatElemNC(t,k+1,k) != 0 ){ /* a pair of complex conjugate eigenvalues */
      a = zMatElemNC(t,k,k);   b = zMatElemNC(t,k,k+1);
      c = zMatElemNC(t,k+1,k); d = zMatElemNC(t,k+1,k+1);
      /* an eigenvector of the 2x2 block, chosen from two candidates to avoid cancellation */
      if( _zSqr(b) + _zSqr(lambda->re-a) + _zSqr(lambda->im) >= _zSqr(lambda->re-d) + _zSqr(lambda->im) + _zSqr(c) ){
        _zComplexCreate( &y[k], b, 0 );
        _zComplexCreate( &y[k+1], lambda->re - a, lambda->im );
      } else{
        _zComplexCreate( &y[k], lambda->re - d, lambda->im );
        _zComplexCreate( &y[k+1], c, 0 );
      }
      _zMatEigDQRBackSubst( t, lambda, y, k-1, k+1, delta );
      _zMatEigDQRPutVec( z, y, k+1, eigbase, k, false );
      _zMatEigDQRPutVec( z, y, k+1, eigbase, k+1, true );
      k += 2;
    } else{
      _zComplexCreate( &y[k], 1, 0 );
      _zMatEigDQRBackSubst( t, lambda, y, k-1, k, delta );
      _zMatEigDQRPutVec( z, y, k, eigbase, k, false );
      k++;
    }
  }
  zFree( y );
  return true;
}

/* eigenvalues and eigenvectors by double-shift QR method. */
#define ZM_EIGDQR_EXCEPTIONAL_SHIFT_1 10
#define ZM_EIGDQR_EXCEPTIONAL_SHIFT_2 20
bool zMatEigDQR(const zMat m, zCVec eigval, zCMat eigbase, int iter)
{
  int r, l, its;
  double norm, s, t, w, scale;
  zMat h, z = NULL;
  zVec nu = NULL, ws = NULL;
  bool ret = false;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return false;
  }
  if( zMatRowSizeNC(m) != zCVecSizeNC(eigval) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_CVEC );
    return false;
  }
  if( eigbase && ( zCMatRowSizeNC(eigbase) != zMatRowSizeNC(m) || zCMatColSizeNC(eigbase) != zMatColSizeNC(m) ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return false;
  }
  h = zMatClone( m );
  nu = zVecAlloc( 3 );
  ws = zVecAlloc( zMatRowSizeNC(m) );
  if( eigbase ) z = zMatAllocSqr( zMatRowSizeNC(m) );
  if( !h || !nu || !ws || ( eigbase && !z ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  /* scale the matrix by a power of two, which is exact, to avoid overflow and underflow */
  for( scale=0, r=0; r<zMatRowSizeNC(h); r++ )
    for( l=0; l<zMatColSizeNC(h); l++ )
      if( fabs( zMatElemNC(h,r,l) ) > scale ) scale = fabs( zMatElemNC(h,r,l) );
  if( scale > 0 ){
    frexp( scale, &l );
    scale = ldexp( 1.0, l );
    for( r=0; r<zMatRowSizeNC(h); r++ )
      for( l=0; l<zMatColSizeNC(h); l++ )
        zMatElemNC(h,r,l) /= scale;
  } else
    scale = 1;
  if( !zMatToHessenberg( h, NULL, z ) ) goto TERMINATE;
  norm = zMatNorm( h );
  ZITERINIT( iter );
  for( r=zMatRowSizeNC(h)-1, its=0; r>=0; ){
    if( ( l = _zMatEigDQRSplit( h, r, norm, its > ZM_EIGDQR_EXCEPTIONAL_SHIFT_2 ) ) == r ){ /* a real eigenvalue */
      zComplexCreate( zCVecElemNC(eigval,r), zMatElemNC(h,r,r), 0 );
      r--;
      its = 0;
      continue;
    }
    if( l == r - 1 ){ /* a 2x2 block */
      _zMatEigDQR2x2( h, z, l, zCVecElemNC(eigval,l), zCVecElemNC(eigval,r) );
      r -= 2;
      its = 0;
      continue;
    }
    if( its >= iter ){
      ZITERWARN( iter );
      goto TERMINATE;
    }
    if( ++its == ZM_EIGDQR_EXCEPTIONAL_SHIFT_1 || its == ZM_EIGDQR_EXCEPTIONAL_SHIFT_2 ){ /* exceptional shift to break a stagnation */
      w = fabs( zMatElemNC(h,r,r-1) ) + fabs( zMatElemNC(h,r-1,r-2) );
      s = 1.5 * w;
      t = w * w;
    } else{ /* eigenvalues of the trailing 2x2 block */
      s = zMatElemNC(h,r-1,r-1) + zMatElemNC(h,r,r);
      t = zMatElemNC(h,r-1,r-1) * zMatElemNC(h,r,r) - zMatElemNC(h,r-1,r) * zMatElemNC(h,r,r-1);
    }
    _zMatEigDQRStep( h, z, l, r, s, t, zVecBuf(nu), zVecBuf(ws) );
  }
  ret = eigbase ? _zMatEigDQRVec( h, z, eigval, eigbase, norm ) : true;
  for( r=0; r<zCVecSizeNC(eigval); r++ ){
    zCVecElemNC(eigval,r)->re *= scale;
    zCVecElemNC(eigval,r)->im *= scale;
  }

 TERMINATE:
  zMatFree( h );
  zMatFree( z );
  zVecFree( nu );
  zVecFree( ws );
  return ret;
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
    zVecNormalize( eigvec, 0, ev );
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

#define ZM_EIGSYS_ITER_NUM 10

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
    mr = i < zMatRowSizeNC(m)-1 ? fabs( zMatElemNC(m,i,i+1) ) : 0;
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
  _zMatSymEigBisecR( m, eigval, eigbase, ZM_EIGSYS_ITER_NUM );
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
    if( zIsTol( zDataMax( zVecBuf(sv)+i, zVecSizeNC(sv)-i, &im ), ZM_MAT_SVD_TOL ) ) break;
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
