/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_dm - optimization tools: descent method.
 */

#include <zm/zm_opt.h>

/* step of numerical differentiation for the i-th component of the variable var, which is eps s_i for
   s_i = max{ |x_i|, ratio |x|inf } (s_i = 1 if x = 0), where xmax is |x|inf.
   the actual step is returned to cancel the round-off error in representing the variable. */
static double _zOptDMDiffStep(const zVec var, int i, double xmax, double eps, double ratio)
{
  double x, xh;

  x = zVecElemNC(var,i);
  xh = x + eps * ( xmax == 0 ? 1.0 : _zMax( fabs( x ), ratio * xmax ) );
  return xh - x;
}

/* gradient vector function */

static zVec _zOptDMGrad(zOptDM *opt, const zVec var, zVec g, void *util){
  return opt->grad( var, g, util );
}
/* numerical gradient vector by the central difference of the function. */
static zVec _zOptDMGradNG(zOptDM *opt, const zVec var, zVec g, void *util)
{
  int i;
  double adv, prv, org, h, xmax;

  xmax = zVecElemAbsMax( var, NULL );
  for( i=0; i<zVecSizeNC(var); i++ ){
    h = _zOptDMDiffStep( var, i, xmax, opt->eps_min, ZM_OPT_DM_DIFF_SCALE_RATIO );
    org = zVecElemNC(var,i);
    zVecSetElemNC( var, i, org + h );
    adv = opt->eval( var, util );
    zVecSetElemNC( var, i, org - h );
    prv = opt->eval( var, util );
    zVecSetElemNC( g, i, 0.5 * ( adv - prv ) / h );
    zVecSetElemNC( var, i, org );
  }
  return g;
}

/* Hessian matix function */

static zMat _zOptDMHess(zOptDM *opt, const zVec var, zMat h, void *util){
  return opt->hess( var, h, util );
}
/* numerical Hessian matrix by the central difference of the given gradient vector,
   where _p and _q are used as workspaces. */
static zMat _zOptDMHessNG(zOptDM *opt, const zVec var, zMat h, void *util)
{
  int i, j;
  double org, d, xmax;

  xmax = zVecElemAbsMax( var, NULL );
  for( i=0; i<zVecSizeNC(var); i++ ){
    d = _zOptDMDiffStep( var, i, xmax, opt->eps_min, ZM_OPT_DM_DIFF_SCALE_RATIO );
    org = zVecElemNC(var,i);
    zVecSetElemNC( var, i, org + d );
    opt->grad( var, opt->_p, util );
    zVecSetElemNC( var, i, org - d );
    opt->grad( var, opt->_q, util );
    zVecSetElemNC( var, i, org );
    zRawVecSub( zVecBufNC(opt->_p), zVecBufNC(opt->_q), zMatRowBuf(h,i), zMatRowSizeNC(h) );
    zRawVecDivDRC( zMatRowBuf(h,i), 2*d, zMatRowSizeNC(h) );
  }
  for( i=0; i<zMatRowSizeNC(h); i++ ) /* symmetrize */
    for( j=i+1; j<zMatColSizeNC(h); j++ )
      zMatElemNC(h,i,j) = zMatElemNC(h,j,i) = 0.5 * ( zMatElemNC(h,i,j) + zMatElemNC(h,j,i) );
  return h;
}
/* numerical Hessian matrix by the second-order central difference of the function. */
static zMat _zOptDMHessNG2(zOptDM *opt, const zVec var, zMat h, void *util)
{
  int i, j;
  double orgi, orgj, di, dj, f0, fpp, fpm, fmp, fmm, xmax, eps, ratio;

  xmax = zVecElemAbsMax( var, NULL );
  eps = pow( opt->eps_min, 0.75 ); /* the optimal ratio for the second-order difference */
  ratio = sqrt( ZM_OPT_DM_DIFF_SCALE_RATIO ); /* see the comment of ZM_OPT_DM_DIFF_SCALE_RATIO */
  f0 = opt->eval( var, util );
  for( i=0; i<zVecSizeNC(var); i++ ){
    di = _zOptDMDiffStep( var, i, xmax, eps, ratio );
    orgi = zVecElemNC(var,i);
    zVecSetElemNC( var, i, orgi + di );
    fpp = opt->eval( var, util );
    zVecSetElemNC( var, i, orgi - di );
    fmm = opt->eval( var, util );
    zMatSetElemNC( h, i, i, ( fpp - 2*f0 + fmm ) / ( di*di ) );
    for( j=i+1; j<zVecSizeNC(var); j++ ){
      dj = _zOptDMDiffStep( var, j, xmax, eps, ratio );
      orgj = zVecElemNC(var,j);
      zVecSetElemNC( var, i, orgi + di ); zVecSetElemNC( var, j, orgj + dj ); fpp = opt->eval( var, util );
      zVecSetElemNC( var, j, orgj - dj ); fpm = opt->eval( var, util );
      zVecSetElemNC( var, i, orgi - di ); fmm = opt->eval( var, util );
      zVecSetElemNC( var, j, orgj + dj ); fmp = opt->eval( var, util );
      zVecSetElemNC( var, j, orgj );
      zMatSetElemNC( h, i, j, ( fpp - fpm - fmp + fmm ) / ( 4*di*dj ) );
      zMatSetElemNC( h, j, i, zMatElemNC(h,i,j) );
    }
    zVecSetElemNC( var, i, orgi );
  }
  return h;
}

/* line search */

/* the upper bound of step length, where every component of a step is limited by _scale. */
static double _zOptDMStepMax(zOptDM *opt)
{
  return opt->_scale / zVecElemAbsMax( opt->_d, NULL );
}

/* evaluate the function at a trial step length. */
static double _zOptDMStepEval(zOptDM *opt, const zVec var, void *util, double a)
{
  zVecMulNC( opt->_d, -a, opt->_p );
  zVecAddNC( var, opt->_p, opt->_x );
  return opt->eval( opt->_x, util );
}

/* check if a step is negligibly small in the sense of the termination condition of zOptDMSolve. */
static bool _zOptDMStepIsNegligible(zOptDM *opt, const zVec var, double tol, double a)
{
  return a * zVecElemAbsMax( opt->_d, NULL ) <= tol * ( 1 + zVecElemAbsMax( var, NULL ) );
}

/* stay at the current variable, since no further decrease of the function is detected. */
static zVec _zOptDMStepStay(zOptDM *opt, const zVec var, double e0, double *e1)
{
  zVecZero( opt->_p );
  zVecZero( opt->_q );
  zVecCopyNC( var, opt->_x );
  *e1 = e0;
  return var;
}

/* step with a fixed length without line search, where _scale is used as the step length. */
static zVec _zOptDMStepFixed(zOptDM *opt, const zVec var, double tol, void *util, double e0, double *e1)
{
  zVecMulNC( opt->_d, -opt->_scale, opt->_p ); /* descent step */
  zVecAddNC( var, opt->_p, opt->_x );
  *e1 = opt->eval( opt->_x, util );
  if( zIsInf( *e1 ) || zIsNan( *e1 ) ){
    ZRUNERROR( ZM_ERR_OPT_CANNOTFINDSTEP );
    return NULL;
  }
  opt->_grad( opt, opt->_x, opt->_r, util );
  zVecSubNC( opt->_r, opt->_g, opt->_q );
  zVecCopyNC( opt->_r, opt->_g );
  return zVecCopyNC( opt->_x, var );
}

/* line search by More-Thuente's method
 * (J. J. More and D. J. Thuente, Line Search Algorithms with Guaranteed Sufficient Decrease,
 *  ACM Transactions on Mathematical Software, 20(3), pp. 286-307, 1994)
 * A step length a along the line phi(a) = f(var - a d) is searched so as to satisfy the strong
 * Wolfe's conditions
 *  phi(a) <= phi(0) + _df a phi'(0) (sufficient decrease), and
 *  |phi'(a)| <= _cf |phi'(0)| (curvature condition).
 * The auxiliary function psi(a) = phi(a) - phi(0) - _df a phi'(0) is used instead of phi until
 * a trial step length satisfies psi(a) <= 0 and phi'(a) >= 0.
 */
#define ZM_OPT_DM_MT_EVAL_MAX 20

/* evaluate the function and its derivative along the line at a trial step length. */
static double _zOptDMStepMTEval(zOptDM *opt, const zVec var, void *util, double a, double *g)
{
  double f;

  f = _zOptDMStepEval( opt, var, util, a );
  if( zIsInf( f ) || zIsNan( f ) ) return f;
  opt->_grad( opt, opt->_x, opt->_r, util );
  *g = -zVecInnerProdNC( opt->_r, opt->_d );
  return f;
}

/* minimizer of the cubic that interpolates f1, g1 at a1 and f2, g2 at a2.
 * the cubic is f1 + g1 t + c2 t^2 + c1/3 t^3 for t = a - a1, and its minimizer is
 * t = (-c2+sqrt(c2^2-c1 g1))/c1, which is rationalized to be valid also for c1=0.
 * the coefficient c1 is stored in c if c is not the null pointer. */
static double _zOptDMStepMTCubic(double a1, double f1, double g1, double a2, double f2, double g2, double *c)
{
  double da, v1, v2, c1, c2;

  da = a2 - a1;
  v1 = f2 - f1 - g1*da;
  v2 = g2 - g1;
  c1 = 3 * ( v2*da - 2*v1 ) / ( da*da*da );
  c2 =     ( 3*v1 - v2*da ) / ( da*da );
  if( c ) *c = c1;
  return a1 - g1 / ( c2 + sqrt( c2*c2 - c1*g1 ) );
}

/* minimizer of the quadratic that interpolates f1, g1 at a1 and f2 at a2. */
static double _zOptDMStepMTQuadratic(double a1, double f1, double g1, double a2, double f2)
{
  double da;

  da = a2 - a1;
  return a1 - 0.5 * g1*da*da / ( f2 - f1 - g1*da );
}

/* minimizer of the quadratic that interpolates g1 at a1 and g2 at a2 (secant method). */
static double _zOptDMStepMTSecant(double a1, double g1, double a2, double g2)
{
  return a1 - g1 * ( a2 - a1 ) / ( g2 - g1 );
}

/* choose the next trial step length from the current endpoints al, au of the interval of
 * uncertainty and the current trial step length at (Section 4 of More and Thuente). */
static double _zOptDMStepMTTrial(double al, double fl, double gl, double au, double fu, double gu, double at, double ft, double gt, bool bracketed, double amax)
{
  double ac, aq, as, c;

  if( ft > fl ){ /* case 1 */
    ac = _zOptDMStepMTCubic( al, fl, gl, at, ft, gt, NULL );
    aq = _zOptDMStepMTQuadratic( al, fl, gl, at, ft );
    return fabs( ac - al ) < fabs( aq - al ) ? ac : 0.5 * ( aq + ac );
  }
  if( gt * gl < 0 ){ /* case 2 */
    ac = _zOptDMStepMTCubic( al, fl, gl, at, ft, gt, NULL );
    as = _zOptDMStepMTSecant( al, gl, at, gt );
    return fabs( ac - at ) >= fabs( as - at ) ? ac : as;
  }
  if( fabs( gt ) <= fabs( gl ) ){ /* case 3 */
    ac = _zOptDMStepMTCubic( al, fl, gl, at, ft, gt, &c );
    as = _zOptDMStepMTSecant( al, gl, at, gt );
    /* the cubic minimizer is adopted only if the cubic tends to infinity in the direction
       of the step, the minimizer is beyond at, and it is closer to at than the secant step. */
    if( !( c * ( at - al ) > 0 && ( ac - at ) * ( at - al ) > 0 && fabs( ac - at ) < fabs( as - at ) ) )
      ac = as;
    if( bracketed )
      return at > al ? _zMin( at + 0.66 * ( au - at ), ac ) : _zMax( at + 0.66 * ( au - at ), ac );
  } else /* case 4 */
    ac = bracketed ? _zOptDMStepMTCubic( at, ft, gt, au, fu, gu, NULL ) : amax;
  if( bracketed ) return ac;
  /* extrapolation within [at+1.1(at-al),at+4(at-al)] */
  if( zIsNan( ac ) ) ac = amax;
  return _zMin( zLimit( ac, at + 1.1 * ( at - al ), at + 4 * ( at - al ) ), amax );
}

static zVec _zOptDMStepMT(zOptDM *opt, const zVec var, double tol, void *util, double e0, double *e1)
{
  double g0, gtest, gtol, amax, m = 0, w, w1, wn, afin = HUGE_VAL;
  double al, fl, gl, au = 0, fu = 0, gu = 0, at, ft, gt = 0, ae, an;
  bool bracketed = false, stage1 = true;
  int i;

  if( zIsNan( e0 ) ){
    ZRUNWARN( ZM_WARN_OPT_BADINI );
    return NULL;
  }
  g0 = -zVecInnerProdNC( opt->_g, opt->_d ); /* phi'(0) < 0 */
  gtest = opt->_df * g0;
  gtol = -opt->_cf * g0;
  amax = _zOptDMStepMax( opt );
  al = 0; fl = e0; gl = g0;
  at = _zMin( 1.0, amax );
  w = w1 = HUGE_VAL;
  for( i=0; i<ZM_OPT_DM_MT_EVAL_MAX; i++ ){
    ft = _zOptDMStepMTEval( opt, var, util, ( ae = at ), &gt );
    if( zIsInf( ft ) || zIsNan( ft ) || zIsInf( gt ) || zIsNan( gt ) ){
      /* shrink the upper bound of step length if the function is not evaluated at the trial point */
      at = amax = al + 0.5 * ( at - al );
      continue;
    }
    if( at < afin ) afin = at; /* the shortest step length where the function is evaluated */
    if( ft <= e0 + at*gtest && fabs( gt ) <= gtol ) goto TERMINATE; /* strong Wolfe's conditions */
    if( at == amax && ft <= e0 + at*gtest && gt < gtest ) goto TERMINATE; /* the upper bound of step length */
    if( stage1 && ft <= e0 + at*gtest && gt >= 0 ){
      stage1 = false; /* phi is used instead of psi hereafter */
      m = 0;
    } else
    if( stage1 ) m = gtest; /* psi(a) = phi(a) - m a, psi'(a) = phi'(a) - m except a constant */
    an = _zOptDMStepMTTrial( al, fl-m*al, gl-m, au, fu-m*au, gu-m, at, ft-m*at, gt-m, bracketed, amax );
    /* update the interval of uncertainty (Section 2 of More and Thuente) */
    if( ft - m*at > fl - m*al ){ /* case U1 */
      au = at; fu = ft; gu = gt;
      bracketed = true;
    } else{
      if( ( gt - m ) * ( al - at ) < 0 ){ /* case U3 */
        au = al; fu = fl; gu = gl;
        bracketed = true;
      } /* otherwise, case U2 */
      al = at; fl = ft; gl = gt;
    }
    if( bracketed ){
      if( ( wn = fabs( au - al ) ) <= zTOL * _zMax( fabs( al ), fabs( au ) ) ) break; /* no further progress */
      /* bisection if the interval does not shrink enough by two trials or the trial step length is not inside */
      if( wn > 0.66 * w1 || zIsNan( an ) || ( an - al ) * ( an - au ) >= 0 )
        an = 0.5 * ( al + au );
      w1 = w;
      w = wn;
    }
    at = an;
  }
  /* the best step length so far is accepted if it satisfies the sufficient decrease condition. */
  if( al <= 0 || fl > e0 + al*gtest ){
    /* no further decrease is detected even by a negligibly small step. */
    if( _zOptDMStepIsNegligible( opt, var, tol, afin ) ) return _zOptDMStepStay( opt, var, e0, e1 );
    ZRUNERROR( ZM_ERR_OPT_CANNOTFINDSTEP );
    return NULL;
  }
  if( ae != al ) ft = _zOptDMStepMTEval( opt, var, util, al, &gt );

 TERMINATE:
  *e1 = ft;
  zVecSubNC( opt->_r, opt->_g, opt->_q );
  zVecCopyNC( opt->_r, opt->_g );
  return zVecCopyNC( opt->_x, var );
}

/* line search by Dennis-Schnabel's backtracking */

/* limit a trial step length of Dennis-Schnabel's backtracking within [0.1a,0.5a] for the previous one a. */
static double _zOptDMStepDSLimit(double at, double a)
{
  if( zIsNan( at ) || at > 0.5*a ) return 0.5*a;
  return at < 0.1*a ? 0.1*a : at;
}

/* test the sufficient decrease condition. a trial step length where the function is not evaluated is rejected. */
static bool _zOptDMStepDSTest(double e0, double e1, double a, double df)
{
  return !zIsInf( e1 ) && !zIsNan( e1 ) && e0 - e1 >= a*df;
}

static zVec _zOptDMStepDS(zOptDM *opt, const zVec var, double tol, void *util, double e0, double *e1)
{
  double a1, a2, w1, w2, s, df, c1, c2;

  if( zIsNan( e0 ) ){
    ZRUNWARN( ZM_WARN_OPT_BADINI );
    return NULL;
  }
  s = zVecInnerProdNC( opt->_g, opt->_d ); /* ascendant direction */
  df = opt->_df * s;
  /* initial test from full-length step */
  a1 = _zMin( 1.0, _zOptDMStepMax( opt ) );
  if( _zOptDMStepDSTest( e0, ( *e1 = _zOptDMStepEval( opt, var, util, a1 ) ), a1, df ) ) goto TERMINATE;
  /* w = ( f(a) - f(0) - a f'(0) ) / a^2, where f'(0) = -s.
     the cubic approximation f(a) = f(0) - s a + c2 a^2 + c1/3 a^3 satisfies w = c2 + c1/3 a. */
  w1 = ( *e1 - e0 + s*a1 ) / ( a1*a1 );
  /* quadratic approximation test */
  a2 = a1;
  w2 = w1;
  a1 = _zOptDMStepDSLimit( 0.5 * s / w1, a2 );
  /* cubic approximation test */
  while( !_zOptDMStepIsNegligible( opt, var, tol, a1 ) ){
    if( _zOptDMStepDSTest( e0, ( *e1 = _zOptDMStepEval( opt, var, util, a1 ) ), a1, df ) ) goto TERMINATE;
    w1 = ( *e1 - e0 + s*a1 ) / ( a1*a1 );
    c1 = 3 * ( w1 - w2 ) / ( a1 - a2 );
    c2 = ( w2*a1 - w1*a2 ) / ( a1 - a2 );
    a2 = a1;
    w2 = w1;
    /* minimizer of the cubic approximation, (-c2+sqrt(c2^2+c1 s))/c1, rationalized to be valid also for c1=0 */
    a1 = _zOptDMStepDSLimit( s / ( c2 + sqrt( c2*c2 + c1*s ) ), a2 );
  }
  /* no further decrease is detected even by a negligibly small step. */
  if( !zIsInf( *e1 ) && !zIsNan( *e1 ) ) return _zOptDMStepStay( opt, var, e0, e1 );
  ZRUNERROR( ZM_ERR_OPT_CANNOTFINDSTEP );
  return NULL;

 TERMINATE:
  opt->_grad( opt, opt->_x, opt->_r, util );
  zVecSubNC( opt->_r, opt->_g, opt->_q );
  zVecCopyNC( opt->_r, opt->_g );
  return zVecCopyNC( opt->_x, var );
}

/* line search by Brent's strict minimization within an interval bracketing a minimum.
 * a minimum is bracketed by expanding or shrinking the step length by the golden ratio from
 * the initial trial step length, and then is searched by Brent's method within the interval.
 */
typedef struct{
  zOptDM *opt;
  const zVec *var;
  void *util;
  double scale; /* the step length is normalized by this value in Brent's method */
} _zOptDMStepData;

/* the function value at a step length, where a non-evaluable value is replaced by the infinity. */
static double _zOptDMStepBrentVal(zOptDM *opt, const zVec var, void *util, double a)
{
  double f;

  f = _zOptDMStepEval( opt, var, util, a );
  return zIsInf( f ) || zIsNan( f ) ? HUGE_VAL : f;
}

static double _zOptDMStepBrentEval(double u, void *util)
{
  _zOptDMStepData *dat;
  dat = util;
  return _zOptDMStepBrentVal( dat->opt, *dat->var, dat->util, u * dat->scale );
}

static zVec _zOptDMStepBrent(zOptDM *opt, const zVec var, double tol, void *util, double e0, double *e1)
{
  double amax, a0, a1, a2, f0, f1, f2, a, fa;
  _zOptDMStepData opt_data;

  if( zIsNan( e0 ) ){
    ZRUNWARN( ZM_WARN_OPT_BADINI );
    return NULL;
  }
  amax = _zOptDMStepMax( opt );
  a0 = 0; f0 = e0;
  f1 = _zOptDMStepBrentVal( opt, var, util, ( a1 = _zMin( 1.0, amax ) ) );
  if( f1 >= f0 ){
    /* shrink the step length until the function decreases */
    do{
      a2 = a1; f2 = f1;
      if( _zOptDMStepIsNegligible( opt, var, tol, ( a1 /= ZM_GOLDENRATIO ) ) ){
        /* no further decrease is detected even by a negligibly small step. */
        if( f2 < HUGE_VAL ) return _zOptDMStepStay( opt, var, e0, e1 );
        ZRUNERROR( ZM_ERR_OPT_CANNOTFINDSTEP );
        return NULL;
      }
      f1 = _zOptDMStepBrentVal( opt, var, util, a1 );
    } while( f1 >= f0 );
  } else{
    /* expand the step length until the function increases */
    for( ; ; a0 = a1, f0 = f1, a1 = a2, f1 = f2 ){
      if( a1 >= amax ){ /* the step length is not expanded beyond the upper bound */
        a2 = a1; f2 = f1;
        f1 = _zOptDMStepBrentVal( opt, var, util, ( a1 = a0 + ( a2 - a0 ) / ZM_GOLDENRATIO ) );
        if( f1 < f2 ) break;
        fa = _zOptDMStepEval( opt, var, util, a2 ); /* the function decreases up to the upper bound */
        goto TERMINATE;
      }
      f2 = _zOptDMStepBrentVal( opt, var, util, ( a2 = _zMin( a1 + ZM_GOLDENRATIO * ( a1 - a0 ), amax ) ) );
      if( f2 >= f1 ) break;
    }
  }
  /* a minimum is bracketed in [a0,a2], where f(a1) < f(a0) and f(a1) <= f(a2) */
  opt_data.opt = opt;
  opt_data.var = &var;
  opt_data.util = util;
  /* the step length is normalized by the upper end of the interval, since the tolerance of
     zOptLineBrent includes an absolute value. */
  opt_data.scale = a2;
  a = zOptLineBrent( _zOptDMStepBrentEval, a0 / a2, 1.0, &opt_data, 100 ) * a2;
  if( !( ( fa = _zOptDMStepBrentVal( opt, var, util, a ) ) <= f1 ) )
    fa = _zOptDMStepEval( opt, var, util, a1 ); /* the best step length while bracketing */

 TERMINATE:
  *e1 = fa;
  opt->_grad( opt, opt->_x, opt->_r, util );
  zVecSubNC( opt->_r, opt->_g, opt->_q );
  zVecCopyNC( opt->_r, opt->_g );
  return zVecCopyNC( opt->_x, var );
}

/* update pseudo-Hessian matix */

static zMat _zOptDMUpdateNone(zOptDM *opt, int count){
  return opt->_h;
}
/* the coefficient of conjugate gradient method by Fletcher-Reeves' formula b = |g|^2 / |g_prev|^2,
   where g_prev = g - q. it restarts from the steepest descent direction every dim iterations. */
static zMat _zOptDMUpdateCG(zOptDM *opt, int count)
{
  double gp;

  zVecSubNC( opt->_g, opt->_q, opt->_r ); /* gradient vector at the previous iteration */
  opt->_b = ( ( count + 1 ) % zVecSizeNC(opt->_d) == 0 || ( gp = zVecSqrNorm( opt->_r ) ) == 0 ) ?
    0 : zVecSqrNorm( opt->_g ) / gp;
  return opt->_h;
}
/* the inverse of the Hessian matrix H is updated by the step p and the difference of gradients q.
   the update is skipped unless p^T q > 0, since H would lose the positive definiteness. */
static zMat _zOptDMUpdateDFP(zOptDM *opt, int count) /* Davidon-Fletcher-Powell */
{
  double k, l;

  /* H += p p^T / p^T q - H q q^T H / q^T H q */
  if( ( k = zVecInnerProdNC( opt->_p, opt->_q ) ) <= 0 ) return opt->_h;
  zMulMatVecNC( opt->_h, opt->_q, opt->_r );
  if( ( l = zVecInnerProdNC( opt->_q, opt->_r ) ) <= 0 ) return opt->_h;
  zMatCatDyadNC( opt->_h, 1.0/k, opt->_p, opt->_p );
  return zMatCatDyadNC( opt->_h,-1.0/l, opt->_r, opt->_r );
}
static zMat _zOptDMUpdateBFGS(zOptDM *opt, int count) /* Broyden-Fletcher-Goldfalb-Shanno */
{
  int i, j;
  double k, l, pi, pj, ri, rj;

  /* H += ( ( 1 + q^T H q / p^T q ) p p^T - p q^T H - H q p^T ) / p^T q */
  if( ( k = zVecInnerProdNC( opt->_p, opt->_q ) ) <= 0 ) return opt->_h;
  k = 1.0 / k;
  zMulMatVecNC( opt->_h, opt->_q, opt->_r );
  l = zVecInnerProdNC( opt->_q, opt->_r )*k + 1.0;
  for( i=0; i<zVecSizeNC(opt->_p); i++ ){
    pi = zVecElemNC( opt->_p, i );
    ri = zVecElemNC( opt->_r, i );
    for( j=0; j<zVecSizeNC(opt->_p); j++ ){
      pj = zVecElemNC( opt->_p, j );
      rj = zVecElemNC( opt->_r, j );
      zMatElemNC(opt->_h,i,j) += ( l*pi*pj - ri*pj - rj*pi ) * k;
    }
  }
  return opt->_h;
}

/* descent vector */

/* Steepest descent method */
static zVec _zOptDMVecSD(zOptDM *opt, const zVec var, zVec d, void *util)
{
  return zVecCopyNC( opt->_g, d );
}

/* damping coefficient of Levenberg-Marquardt method: the norm of the gradient vector */
static double _zOptDMDampGrad(zOptDM *opt, void *dummy)
{
  return zVecNorm( opt->_g );
}

/* damping coefficient of Levenberg-Marquardt method: the smaller of the norm of the gradient vector
   and the Newton decrement g^T H^-1 g, where _x and _p are used as workspaces. */
static double _zOptDMDampDec(zOptDM *opt, void *dummy)
{
  double mg, md;

  mg = zVecNorm( opt->_g );
  zMatCopyNC( opt->_h, opt->_hw );
  zVecCopyNC( opt->_g, opt->_x );
  if( !zLESolveGaussDST( opt->_hw, opt->_x, opt->_p, opt->_idx, opt->_r ) ) return mg;
  md = zVecInnerProdNC( opt->_g, opt->_p );
  return md > 0 && md < mg ? md : mg;
}

/* Levenberg-Marquardt method */
static zVec _zOptDMVecLM(zOptDM *opt, const zVec var, zVec d, void *util)
{
  int i;
  double m;

  opt->_hess( opt, var, opt->_h, util );
  m = opt->_damp( opt, util );
  for( i=0; i<zVecSizeNC(opt->_g); i++ )
    zMatElemNC(opt->_h,i,i) += m;
  /* the right-hand-side vector is copied to _q since zLESolveGaussDST destroys it. */
  zVecCopyNC( opt->_g, opt->_q );
  return zLESolveGaussDST( opt->_h, opt->_q, opt->_d, opt->_idx, opt->_r );
}

/* variable metric method */
static zVec _zOptDMVecVM(zOptDM *opt, const zVec var, zVec d, void *util)
{
  return zMulMatVecNC( opt->_h, opt->_g, opt->_d );
}

/* conjugate gradient method */
static zVec _zOptDMVecCG(zOptDM *opt, const zVec var, zVec d, void *util)
{
  zVecMulDRC( opt->_d, opt->_b );
  return zVecAddNCDRC( opt->_d, opt->_g );
}

/* constructor and destructor */

/* create an unconstrained nonlinear optimization solver based on descent method. */
zOptDM *zOptDMCreate(zOptDM *opt, int dim, double (* eval)(const zVec,void*), zVec (* grad)(const zVec,zVec,void*), zMat (* hess)(const zVec,zMat,void*))
{
  if( !( opt->eval = eval ) ){
    ZRUNERROR( ZM_ERR_OPT_NOEVALUATOR );
    return NULL;
  }
  opt->_grad = ( opt->grad = grad ) ? _zOptDMGrad : _zOptDMGradNG;
  opt->_hess = ( opt->hess = hess ) ? _zOptDMHess : ( grad ? _zOptDMHessNG : _zOptDMHessNG2 );
  opt->eps_min = ZM_OPT_DM_EPS;
  opt->tol_grad = ZM_OPT_DM_GRAD_MIN;
  opt->_damp = _zOptDMDampGrad;
  opt->_hw = NULL;
  opt->_x = zVecAlloc( dim );
  opt->_d = zVecAlloc( dim );
  opt->_g = zVecAlloc( dim );
  opt->_p = zVecAlloc( dim );
  opt->_q = zVecAlloc( dim );
  opt->_r = zVecAlloc( dim );
  opt->_h = zMatAllocSqr( dim );
  opt->_idx = zIndexCreate( dim );
  if( !opt->_x || !opt->_d || !opt->_g || !opt->_p || !opt->_q || !opt->_r || !opt->_h || !opt->_idx ){
    zOptDMDestroy( opt );
    return NULL;
  }
  opt->_b = 0.0;
  opt->_df = 1.0e-4;
  return zOptDMAssignSD( opt, NULL, 0 );
}

/* assign a step method. scale is the upper bound of components of a step for line search methods,
   and the step length for the fixed step. the default value is used if it is zero. */
static void _zOptDMAssignStep(zOptDM *opt, const char *method, double scale)
{
  if( !method ) goto ZOPT_DM_ASSIGN_STEP_MT;
  if( strcmp( method, "fixed" ) == 0 ){ /* fixed step without line search */
    opt->_step = _zOptDMStepFixed;
    opt->_scale = scale == 0 ? ZM_OPT_DM_DEFAULT_STEP : scale;
    return;
  }
  if( strcmp( method, "DS" ) == 0 ){ /* Dennis-Schnabel's method */
    opt->_step = _zOptDMStepDS;
  } else
  if( strcmp( method, "Brent" ) == 0 ) /* Brent's method */
    opt->_step = _zOptDMStepBrent;
  else{
    ZOPT_DM_ASSIGN_STEP_MT:
    opt->_step = _zOptDMStepMT; /* More-Thuente's method as default */
  }
  opt->_scale = scale == 0 ? ZM_OPT_DM_DEFAULT_SCALE : scale;
}

zOptDM *zOptDMAssignSD(zOptDM *opt, const char *stepmethod, double scale)
{
  opt->_vec = _zOptDMVecSD;
  opt->_update_h = _zOptDMUpdateNone;
  _zOptDMAssignStep( opt, stepmethod, scale );
  opt->_cf = 0.9;
  return opt;
}

zOptDM *zOptDMAssignLM(zOptDM *opt, const char *dampmethod, const char *stepmethod, double scale)
{
  if( dampmethod && strcmp( dampmethod, "dec" ) == 0 ){ /* Newton decrement */
    if( !opt->_hw && !( opt->_hw = zMatAllocSqr( zVecSizeNC(opt->_g) ) ) ){
      ZALLOCERROR();
      return NULL;
    }
    opt->_damp = _zOptDMDampDec;
  } else
    opt->_damp = _zOptDMDampGrad; /* norm of the gradient vector as default */
  opt->_vec = _zOptDMVecLM;
  opt->_update_h = _zOptDMUpdateNone;
  _zOptDMAssignStep( opt, stepmethod, scale );
  opt->_cf = 0.9;
  return opt;
}

zOptDM *zOptDMAssignVM(zOptDM *opt, const char *updatemethod, const char *stepmethod, double scale)
{
  opt->_vec = _zOptDMVecVM;
  _zOptDMAssignStep( opt, stepmethod, scale );
  opt->_update_h = ( updatemethod && strcmp( updatemethod, "DFP" ) == 0 ) ?
    _zOptDMUpdateDFP : _zOptDMUpdateBFGS;
  opt->_cf = 0.9;
  return opt;
}

zOptDM *zOptDMAssignCG(zOptDM *opt, const char *stepmethod, double scale)
{
  opt->_vec = _zOptDMVecCG;
  _zOptDMAssignStep( opt, stepmethod ? stepmethod : "Brent", scale ); /* Brent's method as default */
  opt->_update_h = _zOptDMUpdateCG;
  opt->_cf = 0.4;
  return opt;
}

void zOptDMDestroy(zOptDM *opt)
{
  zVecFreeAtOnce( 6, opt->_x, opt->_d, opt->_g, opt->_p, opt->_q, opt->_r );
  zMatFree( opt->_h );
  zMatFree( opt->_hw );
  zIndexFree( opt->_idx );
  opt->eval = NULL;
  opt->grad = NULL;
  opt->hess = NULL;
  opt->_grad = NULL;
  opt->_hess = NULL;
  opt->_vec = NULL;
  opt->_step = NULL;
  opt->_update_h = NULL;
  opt->_damp = NULL;
}

/* the iteration terminates if
 * - the gradient vector becomes small relatively to the initial one, namely, |g|inf <= tol_grad |g0|inf,
 * - the step becomes small relatively to the variable, namely, |p|inf <= tol ( 1 + |x|inf ), or
 * - the number of iterations exceeds the given upper limit iter.
 * it also fails if the line search does not find any step. */
int zOptDMSolve(zOptDM *opt, const zVec var, void *util, double tol, int iter, double *eval)
{
  int i;
  double _e_dummy, e, gmin;

  if( !eval ) eval = &_e_dummy;
  zMatIdent( opt->_h );
  opt->_b = 0; /* for conjugate gradient method */
  ZITERINIT( iter );
  *eval = opt->eval( var, util );
  opt->_grad( opt, var, opt->_g, util );
  if( ( gmin = zVecElemAbsMax( opt->_g, NULL ) ) == 0 )
    return 0; /* no need for iteration. */
  gmin *= opt->tol_grad;
  for( i=0; i<iter; i++ ){
    opt->_vec( opt, var, opt->_d, util ); /* descent vector */
    if( zVecInnerProdNC( opt->_g, opt->_d ) <= 0 ){
      /* restart from the steepest descent direction if not descent. */
      zMatIdent( opt->_h );
      zVecCopyNC( opt->_g, opt->_d );
    }
    if( !opt->_step( opt, var, tol, util, *eval, &e ) ) return -1; /* line search failed. */
    *eval = e;
    if( zVecElemAbsMax( opt->_g, NULL ) <= gmin ||
        zVecElemAbsMax( opt->_p, NULL ) <= tol * ( 1 + zVecElemAbsMax( var, NULL ) ) )
      return i + 1; /* succeed. */
    opt->_update_h( opt, i ); /* update quasi-hessian matrix */
  }
  ZITERWARN( iter );
  return -1;
}
