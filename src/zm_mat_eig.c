/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mat_eig - eigensystem of matrices.
 */

#include <zm/zm_mat_eig.h>

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
    h = NULL;
    goto TERMINATE;
  }
  for( k=0; k<n-2; k++ ){
    len = n - k - 1;
    for( i=0; i<len; i++ ) /* the k-th column below the subdiagonal component */
      zVecSetElemNC( w, i, zMatElemNC(h,k+1+i,k) );
    /* a reflection is skipped only if the column is exactly zero. */
    if( !_zMatHouseholderUnitNormal( zVecBuf(w), len, zVecBuf(nu), &s ) ) continue;
    _zMatHouseholderL( h, zVecBuf(nu), len, k+1, k+1, n-1, zVecBuf(w) );
    _zMatHouseholderR( h, zVecBuf(nu), len, k+1, 0, n-1 );
    /* the k-th column is theoretically mapped to s e1. */
    zMatSetElemNC( h, k+1, k, s );
    for( i=k+2; i<n; i++ ) zMatSetElemNC( h, i, k, 0 );
    if( p ) _zMatHouseholderR( p, zVecBuf(nu), len, k+1, 0, n-1 );
  }
 TERMINATE:
  zVecFree( nu );
  zVecFree( w );
  return h;
}

/* eigenvalue analysis by double-shift QR method by J. G. F. Francis. */

/* find the top row of the lowest unreduced block of a Hessenberg matrix. */
static int _zMatEigDQRSplit(zMat h, int r, double norm, bool stagnant)
{
  int l;
  double s, a;

  for( l=r; l>0; l-- ){
    if( ( s = fabs( zMatElemNC(h,l-1,l-1) ) + fabs( zMatElemNC(h,l,l) ) ) == 0 ) s = norm;
    if( ( a = fabs( zMatElemNC(h,l,l-1) ) ) <= DBL_EPSILON * s ||
        ( stagnant && a <= zMatRowSizeNC(h) * DBL_EPSILON * norm ) ){
      zMatSetElemNC( h, l, l-1, 0 );
      break;
    }
  }
  return l;
}

/* one step of implicit double-shift QR iteration (bulge chasing) on an unreduced block.
 * The shifts mu1 and mu2 are the eigenvalues of a 2x2 matrix [a b; c d], which are given by a, d and
 * bc = b c. The first column of (H-mu1 I)(H-mu2 I) = H^2 - (a+d) H + (ad-bc) I is computed with the
 * differences of the diagonal components from a and d in order to avoid cancellation when the shifts are
 * close to the diagonal components (R. S. Martin, G. Peters and J. H. Wilkinson, 1970). */
static void _zMatEigDQRStep(zMat h, zMat z, int l, int r, double a, double d, double bc, double *nu, double *w)
{
  int n, k, nr, cmax, rmin;
  double x[3], sv;

  n = zMatRowSizeNC(h);
  cmax = z ? n - 1 : r;
  rmin = z ? 0 : l;
  x[0] = ( zMatElemNC(h,l,l) - a ) * ( zMatElemNC(h,l,l) - d ) - bc + zMatElemNC(h,l,l+1) * zMatElemNC(h,l+1,l);
  x[1] = zMatElemNC(h,l+1,l) * ( ( zMatElemNC(h,l,l) - a ) + ( zMatElemNC(h,l+1,l+1) - d ) );
  x[2] = zMatElemNC(h,l+1,l) * zMatElemNC(h,l+2,l+1);
  for( k=l; k<r; k++ ){
    nr = _zMin( 3, r - k + 1 );
    if( k > l ){ /* the bulge to be chased */
      x[0] = zMatElemNC(h,k,k-1);
      x[1] = zMatElemNC(h,k+1,k-1);
      x[2] = nr == 3 ? zMatElemNC(h,k+2,k-1) : 0;
    }
    if( !_zMatHouseholderUnitNormal( x, nr, nu, &sv ) ) continue;
    _zMatHouseholderL( h, nu, nr, k, k, cmax, w );
    if( k > l ){
      zMatSetElemNC( h, k, k-1, sv );
      zMatSetElemNC( h, k+1, k-1, 0 );
      if( nr == 3 ) zMatSetElemNC( h, k+2, k-1, 0 );
    }
    _zMatHouseholderR( h, nu, nr, k, rmin, _zMin( k+3, r ) );
    if( z ) _zMatHouseholderR( z, nu, nr, k, 0, n-1 );
  }
}

/* rotate a quasi-upper triangular matrix by a Givens rotation. */
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

/* eigenvalues of a deflated 2x2 block. */
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
    d += hd;
    q = sqrt( -disc );
    _zComplexCreate( e1, d, q );
    _zComplexCreate( e2, d,-q );
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
      v1 = b;
      v2 = l1 - a;
    } else{
      v1 = l1 - d;
      v2 = c;
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

/* backward substitution to find an eigenvector of a quasi-upper triangular matrix for an eigenvalue. */
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
      zComplexCreate( &a11, zMatElemNC(t,i,i) - lambda->re, -lambda->im );
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

/* transform an eigenvector of a quasi-upper triangular matrix to that of the original matrix,
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
        zComplexCreate( &y[k+1], lambda->re - a, lambda->im );
      } else{
        zComplexCreate( &y[k], lambda->re - d, lambda->im );
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
  double norm, a, d, bc, w, scale;
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
  if( !h || !nu || !ws || ( eigbase && !z ) ) goto TERMINATE;

  scale = _zMatScalePow2( h );
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
    if( ++its == ZM_EIGDQR_EXCEPTIONAL_SHIFT_1 || its == ZM_EIGDQR_EXCEPTIONAL_SHIFT_2 ){
      /* exceptional shift to break a stagnation, which is relative to the last diagonal component */
      w = fabs( zMatElemNC(h,r,r-1) ) + fabs( zMatElemNC(h,r-1,r-2) );
      a = d = zMatElemNC(h,r,r) + 0.75 * w;
      bc = -0.4375 * w * w;
    } else{ /* eigenvalues of the trailing 2x2 block */
      a = zMatElemNC(h,r-1,r-1);
      d = zMatElemNC(h,r,r);
      bc = zMatElemNC(h,r-1,r) * zMatElemNC(h,r,r-1);
    }
    _zMatEigDQRStep( h, z, l, r, a, d, bc, zVecBuf(nu), zVecBuf(ws) );
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

/* seed of the random number generator for the initial vector of the power method. */
#define ZM_EIG_POWER_SEED 1
/* minimum ratio of the norm of an orthogonalized initial vector to the original one, and the maximum
 * number of trials to generate it. */
#define ZM_EIG_POWER_ORTH_MIN 1.0e-3
#define ZM_EIG_POWER_ORTH_TRIAL 100

/* power method: the dominant eigenvalue mu of a linear map v -> w and the corresponding unit eigenvector.
 * The initial vector is random to avoid being orthogonal to the dominant eigenvector, where a local random
 * number generator with a fixed seed is used, so that the result is reproducible and the random number
 * sequence of the caller is not disturbed.
 * The iteration is regarded as converged if the residual |w - mu v| for a unit vector v is not greater
 * than zTOL times |mu|. It returns 1 if it converges, 0 if it does not converge within iter steps, or -1
 * if it fails to allocate working memory. In the latter two cases, the last estimates are stored. */
static int _zMatEigPowerIter(void (* mul)(void*,const zVec,zVec), void (* orth)(void*,zVec), void *util, zVec v, int iter, double *mu)
{
  zRandMT mt; /* a local random number generator */
  zVec w, r;
  double norm, norm0;
  int i, k, ret = -1;

  *mu = 0;
  w = zVecAlloc( zVecSizeNC(v) );
  r = zVecAlloc( zVecSizeNC(v) );
  if( !w || !r ) goto TERMINATE;
  zRandInitSeedMT( &mt, ZM_EIG_POWER_SEED );
  for( k=0; ; k++ ){ /* a random initial unit vector, which is orthogonalized by orth if given */
    for( i=0; i<zVecSizeNC(v); i++ )
      zVecSetElemNC( v, i, zRandMTF( &mt, -1, 1 ) );
    norm0 = zVecNorm( v );
    if( orth ) orth( util, v );
    /* retry if the vector is almost in the subspace excluded by orth */
    if( ( norm = zVecNorm( v ) ) > ZM_EIG_POWER_ORTH_MIN * norm0 && norm > 0 ) break;
    if( k >= ZM_EIG_POWER_ORTH_TRIAL ){
      ret = -1; /* failure, which never happens unless a numerical breakdown occurs */
      goto TERMINATE;
    }
  }
  zVecDivDRC( v, norm );
  for( ret=0, i=0; i<iter; i++ ){
    mul( util, v, w );
    *mu = zVecInnerProdNC( v, w );
    zVecCatNC( w, -*mu, v, r );
    if( zVecNorm( r ) <= zTOL * fabs( *mu ) ){ /* including the case where w = O */
      /* the latest iterate, which is more accurate than v, is adopted unless it is the zero vector. */
      if( ( norm = zVecNorm( w ) ) > 0 ) zVecDivNC( w, norm, v );
      ret = 1;
      break;
    }
    zVecDivNC( w, zVecNorm( w ), v );
  }
 TERMINATE:
  zVecFree( w );
  zVecFree( r );
  return ret;
}

/* multiply a vector by a matrix for the power method. */
static void _zMatEigPowerMulMat(void *util, const zVec v, zVec w)
{
  zMulMatVecNC( (zMat)util, v, w );
}

/* check if a matrix is square and the size of a vector is consistent with it. */
#define _zMatEigPowerCheckSize(m,eigvec) do{\
  if( !zMatIsSqr( m ) ){\
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );\
    return 0;\
  }\
  if( zMatColSizeNC(m) != zVecSizeNC(eigvec) ){\
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );\
    return 0;\
  }\
} while(0)

/* calculate the dominant eigenvalue. */
double zMatEigPower(const zMat m, zVec eigvec, int iter)
{
  double eigval;

  _zMatEigPowerCheckSize( m, eigvec );
  ZITERINIT( iter );
  if( _zMatEigPowerIter( _zMatEigPowerMulMat, NULL, m, eigvec, iter, &eigval ) == 0 ) ZITERWARN( iter );
  return eigval;
}

/* LU decomposition of a matrix for the inverse power method. */
typedef struct{
  zMat l, u;
  zIndex index;
  zVec y;
} _zMatEigPowerInvLU;

/* LU-decompose a square matrix for the inverse power method. It returns 1 if it succeeds, 0 if the
 * matrix is singular, or -1 if it fails to allocate working memory. */
static int _zMatEigPowerInvLUAlloc(_zMatEigPowerInvLU *lu, const zMat m)
{
  int rank;

  lu->y = NULL;
  if( ( rank = zMatDecompLUAlloc( m, &lu->l, &lu->u, &lu->index ) ) < 0 ) return -1;
  if( rank < zMatRowSizeNC(m) ) return 0;
  return ( lu->y = zVecAlloc( zMatRowSizeNC(m) ) ) ? 1 : -1;
}

/* free the LU decomposition for the inverse power method. */
static void _zMatEigPowerInvLUFree(_zMatEigPowerInvLU *lu)
{
  zMatFree( lu->l );
  zMatFree( lu->u );
  zIndexFree( lu->index );
  zVecFree( lu->y );
}

/* multiply a vector by the inverse of a matrix, namely, solve the equation l u w = v. */
static void _zMatEigPowerMulInv(void *util, const zVec v, zVec w)
{
  _zMatEigPowerInvLU *lu;

  lu = (_zMatEigPowerInvLU *)util;
  zLESolveL( lu->l, v, lu->y, lu->index );
  zLESolveU( lu->u, lu->y, w );
}

/* calculate the minimal eigenvalue. */
double zMatEigPowerInv(const zMat m, zVec eigvec, int iter)
{
  _zMatEigPowerInvLU lu;
  double eigval = 0;
  int ret;

  _zMatEigPowerCheckSize( m, eigvec );
  if( ( ret = _zMatEigPowerInvLUAlloc( &lu, m ) ) <= 0 ){
    if( ret == 0 ) ZRUNERROR( ZM_ERR_MAT_SINGULAR );
    zVecZero( eigvec );
    goto TERMINATE;
  }
  ZITERINIT( iter );
  if( ( ret = _zMatEigPowerIter( _zMatEigPowerMulInv, NULL, &lu, eigvec, iter, &eigval ) ) == 0 )
    ZITERWARN( iter );
  if( ret >= 0 )
    eigval = 1.0 / eigval; /* the dominant eigenvalue of the inverse is never zero. */
 TERMINATE:
  _zMatEigPowerInvLUFree( &lu );
  return eigval;
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

/* the range of eigenvalues of a symmetric tridiagonal matrix by Gerschgorin's theorem, which is slightly
 * widened so that the Sturm's counts at the ends include all the eigenvalues. It returns the true value
 * if the range is a point, namely, the matrix is a scalar matrix, the eigenvalues of which are emin. */
static bool _zMatSymEigBisecRange(const zMat m, double *emin, double *emax)
{
  int i;
  double ml, mr, e;

  *emin = *emax = zMatElemNC(m,0,0);
  for( ml=0, i=0; i<zMatRowSizeNC(m); i++, ml=mr ){
    mr = i < zMatRowSizeNC(m)-1 ? fabs( zMatElemNC(m,i,i+1) ) : 0;
    if( ( e = zMatElemNC(m,i,i) + ( ml + mr ) ) > *emax ) *emax = e;
    if( ( e = zMatElemNC(m,i,i) - ( ml + mr ) ) < *emin ) *emin = e;
  }
  if( *emin == *emax ) return true;
  e = 2 * DBL_EPSILON * ( fabs( *emin ) + fabs( *emax ) ) + DBL_MIN;
  *emin -= e;
  *emax += e;
  return false;
}

/* Sturm's function series, which counts the eigenvalues not less than e. A zero pivot is replaced
 * by pivmin, which is tiny relative to the norm of the matrix. */
static int _zMatSymEigBisecSturm(const zMat a, double e, double pivmin)
{
  int i, n = 0;
  double r, p1, p2;

  p1 = e - zMatBuf(a)[0];
  if( p1 <= 0 ) n++;
  for( i=1; i<zMatRowSizeNC(a); i++ ){
    r = zMatElemNC(a,i,i-1);
    p2 = e - zMatElemNC(a,i,i) - r*r / ( p1 == 0 ? pivmin : p1 );
    if( p2 < 0 ) n++;
    p1 = p2;
  }
  return n;
}

/* compute one eigenvalue by bisection shooting. */
static void _zMatSymEigBisecOne(const zMat a, double emin, int nmin, double emax, double pivmin, double *eig)
{
  double e, d, d_old;

  for( d_old=HUGE_VAL; ; d_old=d ){
    e = emin + 0.5 * ( d = emax - emin );
    if( d == d_old ) break;
    if( _zMatSymEigBisecSturm( a, e, pivmin ) == nmin ) emin = e;
    else emax = e;
  }
  *eig = e;
}

/* compute all eigenvalues recursively by bisection method. The eigenvalues in [emin,emax) are stored
 * in eig in descending order. */
static void _zMatSymEigBisecRec(const zMat a, double emin, int nmin, double emax, int nmax, double pivmin, double *eig)
{
  int n;
  double e;

  if( nmin - nmax == 1 )
    return _zMatSymEigBisecOne( a, emin, nmin, emax, pivmin, eig );
  if( !( ( e = 0.5 * ( emax + emin ) ) > emin && e < emax ) ){
    /* the interval cannot be bisected any more, where multiple eigenvalues are regarded as the same. */
    for( n=0; n<nmin-nmax; n++ ) eig[n] = e;
    return;
  }
  n = _zMatSymEigBisecSturm( a, e, pivmin );
  n = _zLimit( n, nmax, nmin ); /* to keep the consistency of the counts against rounding errors */
  if( n > nmax )
    _zMatSymEigBisecRec( a, e, n, emax, nmax, pivmin, &eig[0] );
  if( n < nmin )
    _zMatSymEigBisecRec( a, emin, nmin, e, n, pivmin, &eig[n-nmax] );
}

/* relative threshold to regard eigenvalues as a cluster, the eigenvectors of which are orthogonalized. */
#define ZM_EIG_BISEC_CLUSTER_TOL ( 1.0e-3 )

/* inverse iteration for the eigenvectors in a cluster of eigenvalues. */
typedef struct{
  _zMatEigPowerInvLU lu;
  zMat w; /* eigenvectors stored in rows */
  int c0, c1; /* the eigenvectors in the c0-th to (c1-1)-th rows are in the cluster */
} _zMatSymEigBisecInvIter;

/* orthogonalize a vector against the eigenvectors in the same cluster. */
static void _zMatSymEigBisecOrth(void *util, zVec v)
{
  _zMatSymEigBisecInvIter *ii;
  zVecStruct u;
  int j;

  ii = (_zMatSymEigBisecInvIter *)util;
  for( j=ii->c0; j<ii->c1; j++ ){
    zVecAssignArray( &u, zMatColSizeNC(ii->w), zMatRowBufNC(ii->w,j) );
    zVecOrthogonalize( v, &u, v );
  }
}

/* multiply a vector by the inverse of a matrix, and orthogonalize it against the eigenvectors in
 * the same cluster. */
static void _zMatSymEigBisecMulInv(void *util, const zVec v, zVec w)
{
  _zMatEigPowerMulInv( &((_zMatSymEigBisecInvIter *)util)->lu, v, w );
  _zMatSymEigBisecOrth( util, w );
}

/* compute a unitary matrix consisting of eigenvectors for real-number eigenvalues by the inverse
 * iteration with shifts slightly different from the eigenvalues, where the eigenvectors for a cluster
 * of eigenvalues are orthogonalized with each other (G. Peters and J. H. Wilkinson, 1971). */
static bool _zMatSymEigBisecR(const zMat m, double scale, zVec eigval, zMat eigbase, int iter)
{
  _zMatSymEigBisecInvIter ii;
  zMat ms;
  zVecStruct eigvec;
  double norm, shift, mu;
  int i, ret;
  bool result = false;

  ms = zMatAllocSqr( zMatRowSizeNC(m) );
  ii.w = zMatAllocSqr( zMatRowSizeNC(m) );
  if( !ms || !ii.w ) goto TERMINATE;
  zMatDiv( m, scale, ms ); /* the same scaling with the eigenvalues */
  if( !( ( norm = zMatNorm( ms ) ) > 0 ) ) norm = 1; /* zero matrix */
  for( ii.c0=0, i=0; i<zVecSizeNC(eigval); i++ ){
    if( i > 0 && !( zVecElemNC(eigval,i-1) - zVecElemNC(eigval,i) <= ZM_EIG_BISEC_CLUSTER_TOL * norm ) )
      ii.c0 = i; /* a new cluster */
    ii.c1 = i;
    for( shift=1.0e1*zTOL*norm; ; shift*=10 ){ /* the shift relative to the norm of the matrix */
      zMatDiv( m, scale, ms );
      zMatShift( ms, -zVecElemNC(eigval,i)-shift );
      if( ( ret = _zMatEigPowerInvLUAlloc( &ii.lu, ms ) ) > 0 ) break;
      _zMatEigPowerInvLUFree( &ii.lu );
      if( ret < 0 ) goto TERMINATE;
    }
    zVecAssignArray( &eigvec, zMatColSizeNC(ii.w), zMatRowBufNC(ii.w,i) );
    ret = _zMatEigPowerIter( _zMatSymEigBisecMulInv, _zMatSymEigBisecOrth, &ii, &eigvec, iter, &mu );
    _zMatEigPowerInvLUFree( &ii.lu );
    if( ret < 0 ) goto TERMINATE;
  }
  zMatTNC( ii.w, eigbase );
  result = true;

 TERMINATE:
  zMatFree( ms );
  zMatFree( ii.w );
  return result;
}

/* diagonalize a symmetric matrix by bisection method (J. W. Givens, 1954). */
bool zMatSymEigBisec(const zMat m, zVec eigval, zMat eigbase)
{
  double emin, emax, scale, pivmin;
  int i, nmin, nmax;

  _zMatSymEigCheckSize( m, eigval, eigbase );
  /* eigenvalues */
  zMatCopyNC( m, eigbase ); /* temporary working space */
  scale = _zMatScalePow2( eigbase ); /* to avoid overflow and underflow */
  zMatToHessenberg( eigbase, NULL, NULL );
  pivmin = zMax( DBL_EPSILON * zMatNorm( eigbase ), DBL_MIN );
  if( _zMatSymEigBisecRange( eigbase, &emin, &emax ) ){ /* a scalar matrix */
    zVecSetAll( eigval, emin );
  } else{
    nmin = _zMatSymEigBisecSturm( eigbase, emin, pivmin );
    nmax = _zMatSymEigBisecSturm( eigbase, emax, pivmin );
    _zMatSymEigBisecRec( eigbase, emin, nmin, emax, nmax, pivmin, zVecBuf(eigval) );
  }
  /* eigenvectors of the scaled matrix */
  if( !_zMatSymEigBisecR( m, scale, eigval, eigbase, ZM_EIGSYS_ITER_NUM ) ) return false;
  for( i=0; i<zVecSizeNC(eigval); i++ ) zVecElemNC(eigval,i) *= scale;
  return true;
}

/* cosine c and sine s of the rotation angle phi of Jacobi's method, which is the smaller one
 * (|phi| <= pi/4) of the solutions of tan^2 phi + 2 theta tan phi - 1 = 0. It returns tan phi.
 * The half-angle formulae from cos 2phi = |theta| / r and sin 2phi = 1 / r, where r = sqrt(1+theta^2),
 * give a more orthogonal rotation than c = 1 / sqrt(1+t^2) from t = tan phi. For a large theta,
 * theta^2 would overflow, and r is approximated by |theta|. */
static double _zMatRotJacobiCS(double theta, double *c, double *s)
{
  double r;

  r = fabs( theta ) < 1.0e150 ? sqrt( 1 + theta*theta ) : fabs( theta );
  *c = sqrt( 0.5 + 0.5 * fabs( theta ) / r );
  *s = 0.5 / ( r * *c );
  if( theta < 0 ) *s = -*s;
  return *s / *c;
}

/* Jacobi's rotation to eliminate the (i,j) component of a symmetric matrix with the smaller rotation
 * angle (|angle| <= pi/4), where the diagonal components are updated by Rutishauser's formula. Ref.:
 *  H. Rutishauser, The Jacobi method for real symmetric matrices, Handbook for Automatic Computation
 *  (Wilkinson-Reinsch eds.), Vol. II, 1971.
 */
static void _zMatSymEigJacobiRot(zMat m, zMat eigbase, int i, int j)
{
  int k;
  double t, c, s, tmp1, tmp2;

  t = _zMatRotJacobiCS( 0.5 * ( zMatElemNC(m,j,j) - zMatElemNC(m,i,i) ) / zMatElemNC(m,i,j), &c, &s );
  zMatElemNC(m,i,i) -= t * zMatElemNC(m,i,j);
  zMatElemNC(m,j,j) += t * zMatElemNC(m,i,j);
  zMatSetElemNC( m, i, j, 0 );
  zMatSetElemNC( m, j, i, 0 );
  for( k=0; k<zMatRowSizeNC(m); k++ ){
    /* update of transformation matrix */
    tmp1 = zMatElemNC( eigbase, k, i );
    tmp2 = zMatElemNC( eigbase, k, j );
    zMatSetElemNC( eigbase, k, i, c * tmp1 - s * tmp2 );
    zMatSetElemNC( eigbase, k, j, s * tmp1 + c * tmp2 );
    /* update of the matrix being diagonalized */
    if( k == i || k == j ) continue;
    tmp1 = zMatElemNC( m, k, i );
    tmp2 = zMatElemNC( m, k, j );
    zMatSetElemNC( m, k, i, c * tmp1 - s * tmp2 );
    zMatSetElemNC( m, k, j, s * tmp1 + c * tmp2 );
    zMatSetElemNC( m, i, k, zMatElemNC( m, k, i ) );
    zMatSetElemNC( m, j, k, zMatElemNC( m, k, j ) );
  }
}

/* diagonalize a symmetric matrix by the cyclic Jacobi's method. */
#define _zMatSymEigJacobiSqrtDiag(m,i) sqrt( fabs( zMatElemNC(m,i,i) ) )
bool zMatSymEigJacobi(const zMat m, zVec eigval, zMat eigbase)
{
  int i, j, n;
  double scale, tol_i;
  bool is_complete, ret = false;
  zMat d;

  _zMatSymEigCheckSize( m, eigval, eigbase );
  if( !( d = zMatClone( m ) ) ){ /* to be diagonalized */
    ZALLOCERROR();
    return false;
  }
  scale = _zMatScalePow2( d );
  zMatIdent( eigbase );
  for( n=0; ; n++ ){
    is_complete = true;
    for( i=1; i<zMatRowSizeNC(d); i++ ){
      tol_i = DBL_EPSILON * _zMatSymEigJacobiSqrtDiag(d,i);
      for( j=0; j<i; j++ ){
        /* the threshold relative to the diagonal components (J. Demmel and K. Veselic, 1992) does not
         * depend on the scale of the matrix, and keeps the relative accuracy of small eigenvalues of
         * a positive definite matrix. */
        if( fabs( zMatElemNC(d,i,j) ) <= tol_i * _zMatSymEigJacobiSqrtDiag(d,j) ) continue;
        is_complete = false;
        _zMatSymEigJacobiRot( d, eigbase, i, j );
        tol_i = DBL_EPSILON * _zMatSymEigJacobiSqrtDiag(d,i); /* m_ii is changed by the rotation */
      }
    }
    if( is_complete ) break;
    if( n >= Z_MAX_ITER_NUM ){
      ZITERWARN( Z_MAX_ITER_NUM );
      goto TERMINATE;
    }
  }
  ret = true;

 TERMINATE:
  for( j=0; j<zMatRowSizeNC(m); j++ )
    zVecSetElemNC( eigval, j, zMatElemNC(d,j,j) * scale );
  zMatFree( d );
  return ret;
}

/* one-sided Jacobi's method (M. R. Hestenes, 1958) to orthogonalize the rows of the k x k leading
 * block of a matrix a by rotations from the left side, which are accumulated to the first k columns
 * of u. It returns the false value if it does not converge. */
static bool _zMatSVDJacobi(zMat a, int k, zMat u)
{
  int i, j, l, n;
  double *ai, *aj, alpha, beta, gamma, tol, tol_i, norm_min, c, s, tmp1, tmp2;
  bool is_complete;

  /* the inner product of two orthogonal rows computed in floating-point arithmetic is up to k eps
   * times the product of their norms. */
  tol = k * DBL_EPSILON;
  /* a row not longer than k eps times the Frobenius norm of the block, which is invariant under the
   * rotations, consists of rounding errors, and is not orthogonalized, since the rotations with other
   * rows would only exchange rounding errors. */
  for( norm_min=0, i=0; i<k; i++ )
    norm_min += zRawVecInnerProd( zMatRowBufNC(a,i), zMatRowBufNC(a,i), k );
  norm_min = tol * sqrt( norm_min );
  for( n=0; ; n++ ){
    is_complete = true;
    for( i=1; i<k; i++ ){
      ai = zMatRowBufNC( a, i );
      alpha = zRawVecInnerProd( ai, ai, k );
      for( j=0; j<i; j++ ){
        if( ( tol_i = sqrt( alpha ) ) <= norm_min ) break;
        tol_i *= tol;
        aj = zMatRowBufNC( a, j );
        if( ( beta = sqrt( zRawVecInnerProd( aj, aj, k ) ) ) <= norm_min ) continue;
        gamma = zRawVecInnerProd( ai, aj, k );
        if( fabs( gamma ) <= tol_i * beta ) continue;
        beta *= beta;
        is_complete = false;
        /* the rotation that makes the i-th and j-th rows orthogonal */
        _zMatRotJacobiCS( 0.5 * ( beta - alpha ) / gamma, &c, &s );
        for( l=0; l<k; l++ ){
          tmp1 = ai[l];
          tmp2 = aj[l];
          ai[l] = c * tmp1 - s * tmp2;
          aj[l] = s * tmp1 + c * tmp2;
        }
        for( l=0; l<zMatRowSizeNC(u); l++ ){
          tmp1 = zMatElemNC( u, l, i );
          tmp2 = zMatElemNC( u, l, j );
          zMatSetElemNC( u, l, i, c * tmp1 - s * tmp2 );
          zMatSetElemNC( u, l, j, s * tmp1 + c * tmp2 );
        }
        alpha = zRawVecInnerProd( ai, ai, k ); /* the i-th row is changed by the rotation */
      }
    }
    if( is_complete ) return true;
    if( n >= Z_MAX_ITER_NUM ){
      ZITERWARN( Z_MAX_ITER_NUM );
      return false;
    }
  }
}

/* sort the rows of the k x k leading block of a matrix a in descending order of their norms with the
 * corresponding columns of u, and store the norms into sv. */
static void _zMatSVDSort(zMat a, int k, zMat u, zVec sv)
{
  int i, im;

  for( i=0; i<k; i++ )
    zVecSetElemNC( sv, i, sqrt( zRawVecInnerProd( zMatRowBufNC(a,i), zMatRowBufNC(a,i), k ) ) );
  for( i=0; i<k-1; i++ ){
    zDataMax( zVecBuf(sv)+i, k-i, &im );
    if( ( im += i ) == i ) continue;
    zVecSwapNC( sv, i, im );
    zMatSwapRowNC( a, i, im );
    zMatSwapColNC( u, i, im );
  }
}

/* singular value decomposition. */
int zMatSVD(const zMat m, zMat u, zVec sv, zMat v)
{
  int r, c, k, i, j, rank = -1;
  double scale, d, *nu;
  zMat a, hv = NULL;
  zVec w = NULL;

  r = zMatRowSizeNC(m);
  c = zMatColSizeNC(m);
  k = zMin( r, c );
  if( !zMatIsSqr(u) || zMatRowSizeNC(u) != r || zMatColSizeNC(v) != c || zMatRowCapacity(v) < k ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  if( zVecSizeNC(sv) != r ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return -1;
  }
  a = zMatClone( m );
  hv = zMatAlloc( k, r + c - k ); /* Householder vectors */
  w = zVecAlloc( zMax( r, c ) );
  if( !a || !hv || !w ) goto TERMINATE;
  scale = _zMatScalePow2( a );
  zMatIdent( u );
  zMatZero( hv );
  /* reduce m to a k x k triangular matrix by Householder reflections, namely, m = [L 0] H for r <= c
   * (LQ decomposition), and m = H [R; 0] for r > c (QR decomposition), where H is orthogonal. A null
   * vector nu, which is left for a zero column or row, represents the identity transformation. */
  for( i=0; i<k; i++ ){
    nu = zMatRowBufNC( hv, i ) + i;
    if( r <= c ){ /* the i-th row is mapped to s e_i from the right side */
      if( !_zMatHouseholderUnitNormal( zMatRowBufNC(a,i)+i, c-i, nu, &d ) ) continue;
      if( i < r-1 ) _zMatHouseholderR( a, nu, c-i, i, i+1, r-1 );
      zMatSetElemNC( a, i, i, d );
      for( j=i+1; j<c; j++ ) zMatSetElemNC( a, i, j, 0 );
    } else{ /* the i-th column is mapped to s e_i from the left side */
      for( j=i; j<r; j++ ) zVecSetElemNC( w, j, zMatElemNC(a,j,i) );
      if( !_zMatHouseholderUnitNormal( zVecBuf(w)+i, r-i, nu, &d ) ) continue;
      if( i < c-1 ) _zMatHouseholderL( a, nu, r-i, i, i+1, c-1, zVecBuf(w) );
      zMatSetElemNC( a, i, i, d );
      for( j=i+1; j<r; j++ ) zMatSetElemNC( a, j, i, 0 );
    }
  }
  /* the rows of the triangular matrix are orthogonalized as J^T T = W */
  if( !_zMatSVDJacobi( a, k, u ) ) goto TERMINATE;
  _zMatSVDSort( a, k, u, sv );
  /* the rank and the singular values; a singular value relatively smaller than zTOL is regarded as zero */
  for( rank=0; rank<k; rank++ )
    if( !( zVecElemNC(sv,rank) > zTOL * zVecElemNC(sv,0) ) ) break;
  /* the right singular vectors, which are the normalized rows of W */
  zMatSetRowSizeNC( v, rank );
  zMatZero( v );
  for( i=0; i<rank; i++ )
    for( j=0; j<k; j++ )
      zMatSetElemNC( v, i, j, zMatElemNC(a,i,j) / zVecElemNC(sv,i) );
  for( i=0; i<rank; i++ ) zVecElemNC(sv,i) *= scale;
  for( ; i<r; i++ ) zVecSetElemNC( sv, i, 0 );
  if( r <= c ){ /* m = J W H, where the right singular vectors are the rows of W H */
    if( rank > 0 )
      for( i=k-1; i>=0; i-- )
        _zMatHouseholderR( v, zMatRowBufNC(hv,i)+i, c-i, i, 0, rank-1 );
  } else{ /* m = H [J 0; 0 I] [W; 0], where the left singular vectors are the columns of H [J 0; 0 I] */
    for( i=k-1; i>=0; i-- )
      _zMatHouseholderL( u, zMatRowBufNC(hv,i)+i, r-i, i, 0, r-1, zVecBuf(w) );
  }

 TERMINATE:
  zMatFree( a );
  zMatFree( hv );
  zVecFree( w );
  return rank;
}

/* singular values of a matrix by SVD. It returns the rank of the matrix, or -1 if it fails. The singular
 * values are stored in *sv in descending order, which has to be freed by the caller if it succeeds. */
static int _zMatSingularValue(const zMat m, zVec *sv)
{
  zMat u, v;
  int rank = -1;

  u = zMatAllocSqr( zMatRowSizeNC(m) );
  v = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *sv = zVecAlloc( zMatRowSizeNC(m) );
  if( u && v && *sv ) rank = zMatSVD( m, u, *sv, v );
  zMatFree( u );
  zMatFree( v );
  if( rank < 0 ){
    zVecFree( *sv );
    *sv = NULL;
  }
  return rank;
}

/* maximum singular value of a matrix. */
double zMatSingularValueMax(const zMat m)
{
  zVec sv;
  double s;

  if( _zMatSingularValue( m, &sv ) < 0 ) return NAN;
  s = zMatMinSizeNC(m) > 0 ? zVecElemNC(sv,0) : 0;
  zVecFree( sv );
  return s;
}

/* minimum singular value of a matrix. */
double zMatSingularValueMin(const zMat m)
{
  zVec sv;
  int k;
  double s;

  if( _zMatSingularValue( m, &sv ) < 0 ) return NAN;
  s = ( k = zMatMinSizeNC(m) ) > 0 ? zVecElemNC(sv,k-1) : 0;
  zVecFree( sv );
  return s;
}

/* condition number of a matrix. */
double zMatCondNum(const zMat m)
{
  zVec sv;
  int rank;
  double c;

  if( ( rank = _zMatSingularValue( m, &sv ) ) < 0 ) return NAN;
  c = rank > 0 && rank == zMatMinSizeNC(m) ? zVecElemNC(sv,0) / zVecElemNC(sv,rank-1) : HUGE_VAL;
  zVecFree( sv );
  return c;
}
