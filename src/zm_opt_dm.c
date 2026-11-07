/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_dm - optimization tools: descent method.
 */

#include <zm/zm_opt.h>

/* gradient vector function */

static zVec _zOptDMGrad(zOptDM *opt, const zVec var, zVec g, void *util){
  return opt->grad( var, g, util );
}
static zVec _zOptDMGradNG(zOptDM *opt, const zVec var, zVec g, void *util)
{
  int i;
  double adv, prv, org;

  for( i=0; i<zVecSizeNC(var); i++ ){
    org = zVecElemNC(var,i);
    zVecSetElemNC( var, i, org + ZM_OPT_EPS );
    adv = opt->eval( var, util );
    zVecSetElemNC( var, i, org - ZM_OPT_EPS );
    prv = opt->eval( var, util );
    zVecSetElemNC( g, i, 0.5 * ( adv - prv ) / ZM_OPT_EPS );
    zVecSetElemNC( var, i, org );
  }
  return g;
}

/* Hessian matix function */

static zMat _zOptDMHess(zOptDM *opt, const zVec var, zMat h, void *util){
  return opt->hess( var, h, util );
}
static zMat _zOptDMHessNG(zOptDM *opt, const zVec var, zMat h, void *util)
{
  int i;
  double org;
  zVec adg, prg;

  adg = zVecAlloc( zVecSizeNC(var) );
  prg = zVecAlloc( zVecSizeNC(var) );
  if( !adg || !prg ){
    ZALLOCERROR();
    h = NULL;
    goto TERMINATE;
  }
  for( i=0; i<zVecSizeNC(var); i++ ){
    org = zVecElemNC(var,i);
    zVecSetElemNC( var, i, org + ZM_OPT_EPS );
    opt->_grad( opt, var, adg, util );
    zVecSetElemNC( var, i, org - ZM_OPT_EPS );
    opt->_grad( opt, var, prg, util );
    zRawVecSub( zVecBuf(adg), zVecBuf(prg), zMatRowBuf(h,i), zMatRowSizeNC(h) );
    zVecSetElemNC( var, i, org );
  }
  zMatMulDRC( h, 0.5/ZM_OPT_EPS );
 TERMINATE:
  zVecFree( adg );
  zVecFree( prg );
  return h;
}

/* line search */

static bool _zOptDMStepTest(zOptDM *opt, const zVec var, void *util, double *a, double e0, double df, double *e)
{
  double max;

  while( *a > zTOL ){
    zVecMulNC( opt->_d, -*a, opt->_p );
    zVecAddNC( var, opt->_p, opt->_x );
    *e = opt->eval( opt->_x, util );
    if( zIsInf(*e) || zIsNan(*e) )
      *a *= 0.1;
    else
    if( ( max = zVecElemAbsMax( opt->_p, NULL ) ) > opt->_scale )
      *a *= opt->_scale / max;
    else
      break;
  }
  return ( e0 - *e >= *a*df ) ? true : false;
}

/* directly add a step vector to a variable */
static zVec _zOptDMStepPlain(zOptDM *opt, const zVec var, void *util, double e0, double *e1)
{
  zVecRevNC( opt->_d, opt->_p ); /* descent step */
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

  zVecMulNC( opt->_d, -a, opt->_p );
  zVecAddNC( var, opt->_p, opt->_x );
  f = opt->eval( opt->_x, util );
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

static zVec _zOptDMStepMT(zOptDM *opt, const zVec var, void *util, double e0, double *e1)
{
  double g0, gtest, gtol, amax, m = 0, w, w1, wn;
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
  amax = opt->_scale / zVecElemAbsMax( opt->_d, NULL ); /* every component of a step is limited by _scale */
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
#define ZM_OPT_DM_STEP_MIN 0.001

/* limit a trial step length of Dennis-Schnabel's backtracking within [0.1a,0.5a] for the previous one a. */
static double _zOptDMStepDSLimit(double at, double a)
{
  if( zIsNan( at ) || at > 0.5*a ) return 0.5*a;
  return at < 0.1*a ? 0.1*a : at;
}

static zVec _zOptDMStepDS(zOptDM *opt, const zVec var, void *util, double e0, double *e1)
{
  double a1 = 1.0, a2, w1, w2, s, df, c1, c2;

  if( zIsNan( e0 ) ){
    ZRUNWARN( ZM_WARN_OPT_BADINI );
    return NULL;
  }
  s = zVecInnerProdNC( opt->_g, opt->_d ); /* ascendant direction */
  df = opt->_df * s;
  /* initial test from full-length step */
  if( _zOptDMStepTest( opt, var, util, &a1, e0, df, e1 ) ) goto TERMINATE;
  /* w = ( f(a) - f(0) - a f'(0) ) / a^2, where f'(0) = -s.
     the cubic approximation f(a) = f(0) - s a + c2 a^2 + c1/3 a^3 satisfies w = c2 + c1/3 a. */
  w1 = ( *e1 - e0 + s*a1 ) / ( a1*a1 );
  /* quadratic approximation test */
  a2 = a1;
  w2 = w1;
  a1 = _zOptDMStepDSLimit( 0.5 * s / w1, a2 );
  /* cubic approximation test */
  while( a1 > zTOL ){
    if( _zOptDMStepTest( opt, var, util, &a1, e0, df, e1 ) ) goto TERMINATE;
    w1 = ( *e1 - e0 + s*a1 ) / ( a1*a1 );
    c1 = 3 * ( w1 - w2 ) / ( a1 - a2 );
    c2 = ( w2*a1 - w1*a2 ) / ( a1 - a2 );
    a2 = a1;
    w2 = w1;
    /* minimizer of the cubic approximation, (-c2+sqrt(c2^2+c1 s))/c1, rationalized to be valid also for c1=0 */
    a1 = _zOptDMStepDSLimit( s / ( c2 + sqrt( c2*c2 + c1*s ) ), a2 );
  }
  ZRUNERROR( ZM_ERR_OPT_CANNOTFINDSTEP );
  return NULL;

 TERMINATE:
  opt->_grad( opt, opt->_x, opt->_r, util );
  zVecSubNC( opt->_r, opt->_g, opt->_q );
  zVecCopyNC( opt->_r, opt->_g );
  return zVecCopyNC( opt->_x, var );
}

/* line search by Nocedal's backtracking based on Wolfe's condition */
static zVec _zOptDMStepNoc(zOptDM *opt, const zVec var, void *util, double e0, double *e1)
{
  double a, s, df, cf;

  if( zIsNan( e0 ) ){
    ZRUNWARN( ZM_WARN_OPT_BADINI );
    return NULL;
  }
  s = zVecInnerProdNC( opt->_g, opt->_d );
  df = opt->_df * s;
  cf = opt->_cf * s;
  for( a=1.0; a>=ZM_OPT_DM_STEP_MIN; ){
    if( !_zOptDMStepTest( opt, var, util, &a, e0, df, e1 ) ){
      a *= 0.5;
      continue;
    }
    opt->_grad( opt, opt->_x, opt->_r, util );
    if( a > ZM_OPT_DM_STEP_MIN && zVecInnerProdNC( opt->_r, opt->_d ) > cf ){
      a *= 2.1;
      continue;
    }
    zVecSubNC( opt->_r, opt->_g, opt->_q );
    zVecCopyNC( opt->_r, opt->_g );
    return zVecCopyNC( opt->_x, var );
  }
  ZRUNERROR( ZM_ERR_OPT_CANNOTFINDSTEP );
  return NULL;
}

/* line search by Brent's strict minimization (particularly for CG) */
typedef struct{
  zOptDM *opt;
  const zVec *var;
  void *util;
} _zOptDMStepData;

#define _zm_opt_data(u) ( (_zOptDMStepData *)u )
static double _zOptDMStepBrentEval(double x, void *util)
{
  zVecCatNC( *_zm_opt_data(util)->var, -x, _zm_opt_data(util)->opt->_d, _zm_opt_data(util)->opt->_x );
  return _zm_opt_data(util)->opt->eval( _zm_opt_data(util)->opt->_x, _zm_opt_data(util)->util );
}
#undef _zm_opt_data

static zVec _zOptDMStepBrent(zOptDM *opt, const zVec var, void *util, double e0, double *e1)
{
  double a;
  _zOptDMStepData opt_data;

  if( zIsNan( e0 ) ){
    ZRUNWARN( ZM_WARN_OPT_BADINI );
    return NULL;
  }
  opt_data.opt = opt;
  opt_data.var = &var;
  opt_data.util = util;
  a = zOptLineBrent( _zOptDMStepBrentEval, 0, 1, &opt_data, 100 );
  zVecMulNC( opt->_d, -a, opt->_p );
  zVecAddNC( var, opt->_p, opt->_x );
  *e1 = opt->eval( opt->_x, util );
  opt->_grad( opt, opt->_x, opt->_r, util );
  zVecSubNC( opt->_r, opt->_g, opt->_q );
  opt->_b = zVecSqrNorm(opt->_r) / zVecSqrNorm(opt->_g); /* Fletcher-Reeves */
  zVecCopyNC( opt->_r, opt->_g );
  return zVecCopyNC( opt->_x, var );
}

/* update pseudo-Hessian matix */

static zMat _zOptDMUpdateNone(zOptDM *opt, int count){
  if( ( count - 1 ) % zVecSizeNC(opt->_d) == 0 )
    opt->_b = 0; /* for CG restart */
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

/* Levenberg-Marquardt method */
static zVec _zOptDMVecLM(zOptDM *opt, const zVec var, zVec d, void *util)
{
  int i;
  double m;

  opt->_hess( opt, var, opt->_h, util );
  m = zVecSqrNorm( opt->_g );
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

#define Z_OPT_DM_SCALE ( 1.0e3 )
zOptDM *zOptDMCreate(zOptDM *opt, int dim, double scale, double (* eval)(const zVec,void*), zVec (* grad)(const zVec,zVec,void*), zMat (* hess)(const zVec,zMat,void*))
{
  if( !( opt->eval = eval ) ){
    ZRUNERROR( ZM_ERR_OPT_NOEVALUATOR );
    return NULL;
  }
  opt->_grad = ( opt->grad = grad ) ? _zOptDMGrad : _zOptDMGradNG;
  opt->_hess = ( opt->hess = hess ) ? _zOptDMHess : _zOptDMHessNG;
  opt->_scale = ( scale == 0 ) ? Z_OPT_DM_SCALE : scale;
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
  return zOptDMAssignSD( opt, NULL );
}

static void _zOptDMAssignStep(zOptDM *opt, const char *method)
{
  if( !method )
    opt->_step = _zOptDMStepPlain; /* line search unassigned */
  else
  if( strcmp( method, "DS" ) == 0 ) /* Dennis-Schnabel */
    opt->_step = _zOptDMStepDS;
  else
  if( strcmp( method, "Noc" ) == 0 ) /* Nocedal */
    opt->_step = _zOptDMStepNoc;
  else
  if( strcmp( method, "Brent" ) == 0 ) /* Brent */
    opt->_step = _zOptDMStepBrent;
  else
    opt->_step = _zOptDMStepMT; /* More-Thuente */
}

zOptDM *zOptDMAssignSD(zOptDM *opt, const char *stepmethod)
{
  opt->_vec = _zOptDMVecSD;
  opt->_update_h = _zOptDMUpdateNone;
  _zOptDMAssignStep( opt, stepmethod );
  opt->_cf = 0.9;
  return opt;
}

zOptDM *zOptDMAssignLM(zOptDM *opt, const char *stepmethod)
{
  opt->_vec = _zOptDMVecLM;
  opt->_update_h = _zOptDMUpdateNone;
  _zOptDMAssignStep( opt, stepmethod );
  opt->_cf = 0.9;
  return opt;
}

zOptDM *zOptDMAssignVM(zOptDM *opt, const char *stepmethod, const char *updatemethod)
{
  opt->_vec = _zOptDMVecVM;
  _zOptDMAssignStep( opt, stepmethod );
  opt->_update_h = ( updatemethod && strcmp( updatemethod, "DFP" ) == 0 ) ?
    _zOptDMUpdateDFP : _zOptDMUpdateBFGS;
  opt->_cf = 0.9;
  return opt;
}

zOptDM *zOptDMAssignCG(zOptDM *opt)
{
  opt->_vec = _zOptDMVecCG;
  opt->_step = _zOptDMStepBrent;
  opt->_update_h = _zOptDMUpdateNone;
  opt->_cf = 0.4;
  return opt;
}

void zOptDMDestroy(zOptDM *opt)
{
  zVecFreeAtOnce( 6, opt->_x, opt->_d, opt->_g, opt->_p, opt->_q, opt->_r );
  zMatFree( opt->_h );
  zIndexFree( opt->_idx );
  opt->eval = NULL;
  opt->grad = NULL;
  opt->hess = NULL;
  opt->_grad = NULL;
  opt->_hess = NULL;
  opt->_vec = NULL;
  opt->_step = NULL;
  opt->_update_h = NULL;
}

/* the iteration terminates if
 * - the gradient vector becomes small relatively to the initial one, namely, |g|inf <= ZM_OPT_DM_GRAD_MIN |g0|inf,
 * - the step becomes small relatively to the variable, namely, |p|inf <= tol ( 1 + |x|inf ), or
 * - the number of iterations exceeds the given upper limit iter.
 * it also fails if the line search does not find any step. */
int zOptDMSolve(zOptDM *opt, const zVec var, void *util, double tol, int iter, double *eval)
{
  int i;
  double _e_dummy, e, gmin;

  if( !eval ) eval = &_e_dummy;
  zMatIdent( opt->_h );
  ZITERINIT( iter );
  *eval = opt->eval( var, util );
  opt->_grad( opt, var, opt->_g, util );
  if( ( gmin = zVecElemAbsMax( opt->_g, NULL ) ) == 0 )
    return 0; /* no need for iteration. */
  gmin *= ZM_OPT_DM_GRAD_MIN;
  for( i=0; i<iter; i++ ){
    opt->_vec( opt, var, opt->_d, util ); /* descent vector */
    if( zVecInnerProdNC( opt->_g, opt->_d ) <= 0 ){
      /* restart from the steepest descent direction if not descent. */
      zMatIdent( opt->_h );
      zVecCopyNC( opt->_g, opt->_d );
    }
    if( !opt->_step( opt, var, util, *eval, &e ) ) return -1; /* line search failed. */
    *eval = e;
    if( zVecElemAbsMax( opt->_g, NULL ) <= gmin ||
        zVecElemAbsMax( opt->_p, NULL ) <= tol * ( 1 + zVecElemAbsMax( var, NULL ) ) )
      return i + 1; /* succeed. */
    opt->_update_h( opt, i ); /* update quasi-hessian matrix */
  }
  ZITERWARN( iter );
  return -1;
}
