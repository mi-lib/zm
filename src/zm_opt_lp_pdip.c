/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lp_pdip - optimization tools:
 * linear programming by primal-dual interior-point method.
 */

#include <zm/zm_opt.h>

/* ratio of a step length to the maximum one to keep the variables positive */
#define ZM_LP_PDIP_STEP_RATIO 0.99
/* number of iterations without any decrease of the error, after which the iteration is given up */
#define ZM_LP_PDIP_STALL_ITER 15

/* workspace for primal-dual interior-point method with Mehrotra's predictor-corrector, which solves
 * the problem normalized as c^T x / |c| -> min. subject to A x = b / |b| and x >= 0. */
typedef struct{
  zMat a;
  zVec b, c;            /* normalized constant vector and cost coefficient vector */
  zVec x, y, z;         /* primal variable, Lagrange multiplier and dual slack variable */
  zVec rp, rd, rc;      /* residuals of the primal, dual and complementarity conditions */
  zVec dx, dy, dz;      /* predictor */
  zVec dx2, dy2, dz2;   /* corrector */
  zVec w, t, s, e;      /* working vectors */
  zMat at, aw, l, q;    /* transpose of a, weighted matrix S A (X Z^-1)^1/2 and its LQ decomposition */
  zIndex idx;
  int rank;
} _zLP_PDIP;

/* free the workspace for PD-IP. */
static void _zLP_PDIPFree(_zLP_PDIP *dat)
{
  zVecFreeAtOnce( 18, dat->b, dat->c, dat->x, dat->y, dat->z, dat->rp, dat->rd, dat->rc, dat->dx, dat->dy, dat->dz, dat->dx2, dat->dy2, dat->dz2, dat->w, dat->t, dat->s, dat->e );
  zMatFreeAtOnce( 4, dat->at, dat->aw, dat->l, dat->q );
  zIndexFree( dat->idx );
}

/* allocate the workspace for PD-IP. */
static _zLP_PDIP *_zLP_PDIPAlloc(_zLP_PDIP *dat, const zMat a)
{
  int m, n;

  m = zMatRowSizeNC(a);
  n = zMatColSizeNC(a);
  dat->a = a;
  dat->b = zVecAlloc( m );
  dat->c = zVecAlloc( n );
  dat->x = zVecAlloc( n );
  dat->y = zVecAlloc( m );
  dat->z = zVecAlloc( n );
  dat->rp = zVecAlloc( m );
  dat->rd = zVecAlloc( n );
  dat->rc = zVecAlloc( n );
  dat->dx = zVecAlloc( n );
  dat->dy = zVecAlloc( m );
  dat->dz = zVecAlloc( n );
  dat->dx2 = zVecAlloc( n );
  dat->dy2 = zVecAlloc( m );
  dat->dz2 = zVecAlloc( n );
  dat->w = zVecAlloc( n );
  dat->t = zVecAlloc( m );
  dat->s = zVecAlloc( m );
  dat->e = zVecAlloc( m );
  dat->at = zMatAlloc( n, m );
  dat->aw = zMatAlloc( m, n );
  dat->l = zMatAlloc( m, zMin( m, n ) );
  dat->q = zMatAlloc( zMin( m, n ), n );
  dat->idx = zIndexCreate( m );
  if( !dat->b || !dat->c || !dat->x || !dat->y || !dat->z || !dat->rp || !dat->rd || !dat->rc ||
      !dat->dx || !dat->dy || !dat->dz || !dat->dx2 || !dat->dy2 || !dat->dz2 ||
      !dat->w || !dat->t || !dat->s || !dat->e || !dat->at || !dat->aw || !dat->l || !dat->q || !dat->idx ){
    ZALLOCERROR();
    _zLP_PDIPFree( dat );
    return NULL;
  }
  return dat;
}

/* initial interior point for PD-IP by Mehrotra's heuristic, where x and y are the minimum-norm
 * least-square solutions of A x = b and A^T y = c, respectively, and x and z = c - A^T y are shifted to
 * be positive and to be balanced. */
static bool _zLP_PDIPInit(_zLP_PDIP *dat)
{
  double dx, dz, xz;

  zMatT( dat->a, dat->at );
  if( !zLESolveMP( dat->a, dat->b, NULL, NULL, dat->x ) ||
      !zLESolveMP( dat->at, dat->c, NULL, NULL, dat->y ) ){
    ZRUNERROR( ZM_ERR_OPT_NOINITIALPOINT );
    return false;
  }
  zMulMatVecNC( dat->at, dat->y, dat->z );
  zVecSubNC( dat->c, dat->z, dat->z );
  zVecShiftDRC( dat->x, zMax( -1.5 * zVecElemMin( dat->x, NULL ), 0 ) );
  zVecShiftDRC( dat->z, zMax( -1.5 * zVecElemMin( dat->z, NULL ), 0 ) );
  xz = zVecInnerProd( dat->x, dat->z );
  dx = zVecElemSum( dat->x );
  dz = zVecElemSum( dat->z );
  if( dx > 0 && dz > 0 ){
    zVecShiftDRC( dat->x, 0.5 * xz / dz );
    zVecShiftDRC( dat->z, 0.5 * xz / dx );
  }
  /* in case that x or z vanishes */
  if( !( zVecElemMin( dat->x, NULL ) > 0 ) ) zVecShiftDRC( dat->x, 1 );
  if( !( zVecElemMin( dat->z, NULL ) > 0 ) ) zVecShiftDRC( dat->z, 1 );
  return true;
}

/* residuals of the primal and dual conditions rp = b - A x and rd = c - A^T y - z. */
static void _zLP_PDIPResidual(_zLP_PDIP *dat)
{
  zMulMatVecNC( dat->a, dat->x, dat->rp );
  zVecSubNC( dat->b, dat->rp, dat->rp );
  zMulMatVecNC( dat->at, dat->y, dat->rd );
  zVecAddNCDRC( dat->rd, dat->z );
  zVecSubNC( dat->c, dat->rd, dat->rd );
}

/* decomposition of the normal matrix A X Z^-1 A^T of the Newton's equation. It is symmetrically scaled
 * by S so that the diagonal components are one, and is decomposed as S A X Z^-1 A^T S = L L^T, where L
 * is found by the LQ decomposition of S A (X Z^-1)^1/2 without forming the normal matrix, so that the
 * condition number is not squared. The rank-deficiency due to dependent constraints is detected in
 * the decomposition. */
static bool _zLP_PDIPDecomp(_zLP_PDIP *dat)
{
  int i, j;
  double d;

  for( j=0; j<zVecSizeNC(dat->w); j++ )
    zVecSetElemNC( dat->w, j, sqrt( zVecElemNC(dat->x,j) / zVecElemNC(dat->z,j) ) );
  for( i=0; i<zMatRowSizeNC(dat->a); i++ ){
    for( j=0; j<zMatColSizeNC(dat->a); j++ )
      zMatElemNC(dat->aw,i,j) = zMatElemNC(dat->a,i,j) * zVecElemNC(dat->w,j);
    zVecSetElemNC( dat->s, i, ( d = zRawVecNorm( zMatRowBufNC(dat->aw,i), zMatColSizeNC(dat->aw) ) ) > 0 ? 1.0 / d : 1.0 );
    zRawVecMulDRC( zMatRowBufNC(dat->aw,i), zVecElemNC(dat->s,i), zMatColSizeNC(dat->aw) );
  }
  /* restore the sizes truncated by the last decomposition */
  zMatSetColSizeNC( dat->l, zMin( zMatRowSizeNC(dat->a), zMatColSizeNC(dat->a) ) );
  zMatSetRowSizeNC( dat->q, zMatColSizeNC(dat->l) );
  return ( dat->rank = zMatDecompLQ( dat->aw, dat->l, dat->q, dat->idx ) ) >= 0;
}

/* solve the Newton's equation A dx = rp, A^T dy + dz = rd, Z dx + X dz = rc for dx, dy and dz.
 * dy is a solution of the normal equation A X Z^-1 A^T dy = rp + A Z^-1 ( X rd - rc ). Since the
 * accuracy of A dx = rp is limited by the rounding errors of the normal equation, the solution is
 * refined once with respect to the residual rp - A dx. */
static void _zLP_PDIPSolve(_zLP_PDIP *dat, zVec dx, zVec dy, zVec dz)
{
  int i;

  for( i=0; i<zVecSizeNC(dat->x); i++ )
    zVecSetElemNC( dat->w, i, ( zVecElemNC(dat->x,i) * zVecElemNC(dat->rd,i) - zVecElemNC(dat->rc,i) ) / zVecElemNC(dat->z,i) );
  zMulMatVecNC( dat->a, dat->w, dat->t );
  zVecAddNCDRC( dat->t, dat->rp );
  zVecAmpNCDRC( dat->t, dat->s );
  zLESolveCholeskyGen( dat->l, dat->t, dy, dat->idx );
  zVecAmpNCDRC( dy, dat->s );
  zMulMatVecNC( dat->at, dy, dz );
  zVecSubNC( dat->rd, dz, dz );
  for( i=0; i<zVecSizeNC(dat->x); i++ )
    zVecSetElemNC( dx, i, ( zVecElemNC(dat->rc,i) - zVecElemNC(dat->x,i) * zVecElemNC(dz,i) ) / zVecElemNC(dat->z,i) );
  /* refinement: dy += e, dz -= A^T e, dx += X Z^-1 A^T e, where A X Z^-1 A^T e = rp - A dx */
  zMulMatVecNC( dat->a, dx, dat->t );
  zVecSubNC( dat->rp, dat->t, dat->t );
  zVecAmpNCDRC( dat->t, dat->s );
  zLESolveCholeskyGen( dat->l, dat->t, dat->e, dat->idx );
  zVecAmpNCDRC( dat->e, dat->s );
  zVecAddNCDRC( dy, dat->e );
  zMulMatVecNC( dat->at, dat->e, dat->w );
  zVecSubNCDRC( dz, dat->w );
  for( i=0; i<zVecSizeNC(dat->x); i++ )
    zVecElemNC(dx,i) += zVecElemNC(dat->x,i) * zVecElemNC(dat->w,i) / zVecElemNC(dat->z,i);
}

/* maximum step length a such that v + a dv >= 0. HUGE_VAL is returned if no component decreases. */
static double _zLP_PDIPStepLen(const zVec v, const zVec dv)
{
  int i;
  double a = HUGE_VAL, ai;

  for( i=0; i<zVecSizeNC(v); i++ )
    if( zVecElemNC(dv,i) < 0 && ( ai = -zVecElemNC(v,i) / zVecElemNC(dv,i) ) < a ) a = ai;
  return a;
}

/* error of the normalized problem, which is the maximum of the norms of the residuals rp and rd and the
 * duality gap relative to the scale of its rounding error. */
static double _zLP_PDIPError(_zLP_PDIP *dat)
{
  double gap, scale;

  gap = fabs( zVecInnerProd( dat->c, dat->x ) - zVecInnerProd( dat->b, dat->y ) );
  if( ( scale = zVecNorm( dat->c ) * zVecNorm( dat->x ) + zVecNorm( dat->b ) * zVecNorm( dat->y ) ) > 0 )
    gap /= scale;
  return zMax3( zVecNorm( dat->rp ), zVecNorm( dat->rd ), gap );
}

/* one iteration of Mehrotra's predictor-corrector method. */
static bool _zLP_PDIPUpdate(_zLP_PDIP *dat)
{
  int i, n;
  double mu, mu_aff, sigma, ap, ad;

  n = zVecSizeNC(dat->x);
  mu = zVecInnerProd( dat->x, dat->z ) / n;
  if( !_zLP_PDIPDecomp( dat ) ) return false;
  /* predictor (affine-scaling direction) */
  zVecAmpNC( dat->x, dat->z, dat->rc );
  zVecRevNCDRC( dat->rc );
  _zLP_PDIPSolve( dat, dat->dx, dat->dy, dat->dz );
  ap = zMin( _zLP_PDIPStepLen( dat->x, dat->dx ), 1 );
  ad = zMin( _zLP_PDIPStepLen( dat->z, dat->dz ), 1 );
  for( mu_aff=0, i=0; i<n; i++ )
    mu_aff += ( zVecElemNC(dat->x,i) + ap * zVecElemNC(dat->dx,i) ) * ( zVecElemNC(dat->z,i) + ad * zVecElemNC(dat->dz,i) );
  sigma = zCube( mu_aff / n / mu );
  /* corrector (combined direction) */
  for( i=0; i<n; i++ )
    zVecElemNC(dat->rc,i) -= zVecElemNC(dat->dx,i) * zVecElemNC(dat->dz,i) - sigma * mu;
  _zLP_PDIPSolve( dat, dat->dx2, dat->dy2, dat->dz2 );
  /* step lengths, which keep x and z strictly positive */
  ap = zMin( ZM_LP_PDIP_STEP_RATIO * _zLP_PDIPStepLen( dat->x, dat->dx2 ), 1 );
  ad = zMin( ZM_LP_PDIP_STEP_RATIO * _zLP_PDIPStepLen( dat->z, dat->dz2 ), 1 );
  zVecCatNCDRC( dat->x, ap, dat->dx2 );
  zVecCatNCDRC( dat->y, ad, dat->dy2 );
  zVecCatNCDRC( dat->z, ad, dat->dz2 );
  return zIsFinite( zVecInnerProd( dat->x, dat->z ) ) && zIsFinite( zVecNorm( dat->y ) );
}

/* linear programming solver based on primal-dual interior-point method with Mehrotra's predictor-corrector. */
bool zLPSolvePDIP_PC(const zMat a, const zVec b, const zVec c, zVec x, double *cost)
{
  _zLP_PDIP dat;
  double nb, nc, err, err_min = HUGE_VAL;
  int i, i_min = 0, iter = 0;
  bool ret = false;

  if( !zMatRowVecSizeEqual( a, b ) || !zMatColVecSizeEqual( a, c ) || !zVecSizeEqual( c, x ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return false;
  }
  if( !_zLP_PDIPAlloc( &dat, a ) ) return false;
  /* normalization */
  if( ( nb = zVecNorm( b ) ) == 0 ) nb = 1;
  if( ( nc = zVecNorm( c ) ) == 0 ) nc = 1;
  zVecDiv( b, nb, dat.b );
  zVecDiv( c, nc, dat.c );
  if( !_zLP_PDIPInit( &dat ) ) goto TERMINATE;
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    _zLP_PDIPResidual( &dat );
    if( ( err = _zLP_PDIPError( &dat ) ) <= zTOL ){
      ret = true;
      break;
    }
    if( err < err_min ){
      err_min = err;
      i_min = i;
    } else
    if( i - i_min > ZM_LP_PDIP_STALL_ITER ){ /* the error stops decreasing */
      ZRUNWARN( ZM_ERR_OPT_UNSOLVABLE );
      break;
    }
    if( !_zLP_PDIPUpdate( &dat ) ){
      ZRUNWARN( ZM_ERR_OPT_UNSOLVABLE );
      break;
    }
  }
  if( i == iter ) ZITERWARN( iter );
 TERMINATE:
  zVecMul( dat.x, nb, x );
  if( cost ) *cost = zVecInnerProd( c, x );
  _zLP_PDIPFree( &dat );
  return ret;
}
