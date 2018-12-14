/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_dm - optimization tools: descent method.
 */

#include <zm/zm_opt.h>

static zVec _zOptDMGrad(zOptDM *opt, zVec var, zVec g, void *util);
static zVec _zOptDMGradNG(zOptDM *opt, zVec var, zVec g, void *util);

static zMat _zOptDMHess(zOptDM *opt, zVec var, zMat h, void *util);
static zMat _zOptDMHessNG(zOptDM *opt, zVec var, zMat h, void *util);

static bool _zOptDMStepTest(zOptDM *opt, zVec var, void *util, double *a, double e0, double df, double *e);
static zVec _zOptDMStepPlain(zOptDM *opt, zVec var, void *util, double e0, double *e1);
static zVec _zOptDMStepMT(zOptDM *opt, zVec var, void *util, double e0, double *e1);
static zVec _zOptDMStepDS(zOptDM *opt, zVec var, void *util, double e0, double *e1);
static zVec _zOptDMStepNoc(zOptDM *opt, zVec var, void *util, double e0, double *e1);
static zVec _zOptDMStepBrent(zOptDM *opt, zVec var, void *util, double e0, double *e1);

static zMat _zOptDMUpdateNone(zOptDM *opt, int count);
static zMat _zOptDMUpdateDFP(zOptDM *opt, int count);
static zMat _zOptDMUpdateBFGS(zOptDM *opt, int count);

static zVec _zOptDMVecSD(zOptDM *opt, zVec var, zVec d, void *util);
static zVec _zOptDMVecLM(zOptDM *opt, zVec var, zVec d, void *util);
static zVec _zOptDMVecVM(zOptDM *opt, zVec var, zVec d, void *util);
static zVec _zOptDMVecCG(zOptDM *opt, zVec var, zVec d, void *util);

static void _zOptDMAssignStep(zOptDM *opt, const char *method);

/* *** gradient vector function *** */

zVec _zOptDMGrad(zOptDM *opt, zVec var, zVec g, void *util){
  return opt->grad( var, g, util );
}
zVec _zOptDMGradNG(zOptDM *opt, zVec var, zVec g, void *util)
{
  register int i;
  double adv, prv, org;

  for( i=0; i<zVecSizeNC(var); i++ ){
    org = zVecElem(var,i);
    zVecSetElem( var, i, org+Z_OPT_EPS );
    adv = opt->eval( var, util );
    zVecSetElem( var, i, org-Z_OPT_EPS );
    prv = opt->eval( var, util );
    zVecSetElem( g, i, 0.5 * ( adv - prv ) / Z_OPT_EPS );
    zVecSetElem( var, i, org );
  }
  return g;
}

/* *** Hessian matix function *** */

zMat _zOptDMHess(zOptDM *opt, zVec var, zMat h, void *util){
  return opt->hess( var, h, util );
}
zMat _zOptDMHessNG(zOptDM *opt, zVec var, zMat h, void *util)
{
  register int i;
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
    org = zVecElem(var,i);
    zVecSetElem( var, i, org+Z_OPT_EPS );
    opt->_grad( opt, var, adg, util );
    zVecSetElem( var, i, org-Z_OPT_EPS );
    opt->_grad( opt, var, prg, util );
    zRawVecSub( zVecBuf(adg), zVecBuf(prg), zMatRowBuf(h,i), zMatRowSizeNC(h) );
    zVecSetElem( var, i, org );
  }
  zMatMulDRC( h, 0.5/Z_OPT_EPS );
 TERMINATE:
  zVecFree( adg );
  zVecFree( prg );
  return h;
}

/* *** line search *** */

bool _zOptDMStepTest(zOptDM *opt, zVec var, void *util, double *a, double e0, double df, double *e)
{
  double max;

  while( *a > zTOL ){
    zVecMulNC( opt->_d, -*a, opt->_p );
    zVecAddNC( var, opt->_p, opt->_x );
    *e = opt->eval( opt->_x, util );
    if( zIsInf(*e) || zIsNan(*e) )
      *a *= 0.1;
    else
    if( ( max = zVecAbsMax( opt->_p, NULL ) ) > opt->_scale )
      *a *= opt->_scale / max;
    else
      break;
  }
  return ( e0 - *e >= *a*df ) ? true : false;
}

/* directly add a step vector to a variable */
zVec _zOptDMStepPlain(zOptDM *opt, zVec var, void *util, double e0, double *e1)
{
  zVecAddNCDRC( var, opt->_d );
  *e1 = opt->eval( var, util );
  return var;
}

/* line search by More-Thuente's method */
zVec _zOptDMStepMT(zOptDM *opt, zVec var, void *util, double e0, double *e1)
{
  double al = 0.0, au = 1.0, at, ac, aq, as, da, s, df;
  double fl, gl, gt, c1, c2, v1, v2;

  if( zIsNan(e0) ){
    ZRUNWARN( ZM_WARN_OPT_BADINI );
    return NULL;
  }
  s = zVecInnerProdNC( opt->_g, opt->_d ); /* ascendant direction */
  df = opt->_df * s;
  /* initial test from full-length step */
  if( _zOptDMStepTest( opt, var, util, &au, e0, df, e1 ) ) goto TERMINATE;
  fl = e0;
  gl =-s;
  at = au;
  opt->_grad( opt, opt->_x, opt->_r, util );
  gt =-zVecInnerProdNC( opt->_g, opt->_d );
  da = 1.0;
  while( da > zTOL ){
    /* STEP1: choose well-estimated trial value (at) in [al,au] */
    v1 =*e1 - fl - gl*da;
    v2 = gt - gl;
    c1 = 3 * (-2*v1 + v2*da ) / ( da*da*da );
    c2 =     ( 3*v1 - v2*da ) / ( da*da );
    ac = al + ( -c2 + sqrt(c2*c2-c1*gl) ) / c1; /* cubic minimizer */
    aq = al - 0.5 * gl*da*da / v1;  /* quadratic minimizer (type1) */
    as = al -       gl*da    / v2;  /* quadratic minimizer (type2) */
    if( *e1 > fl )
      at = ( fabs( ac - al ) < fabs( aq - al ) ) ? ac : 0.5*(aq+ac);
    else{
      if( gt * gl < 0 )
        at = ( fabs( ac - at ) >= fabs( as - at ) ) ? ac : as;
      else{
        if( fabs( gt ) <= fabs( gl ) ){
          at = ( fabs( ac - at ) < fabs( as - at ) ) ? ac : as;
          v1 = at + 0.66 * (au-at); /* v1 temporarily used */
          at = ( at > al ) ? zMin( v1, at ) : zMax( v1, at );
        } else
          at = ac;
      }
    }
    /* STEP2: check if (at) satisfies Wolfe's condition */
    if( _zOptDMStepTest( opt, var, util, &at, e0, df, e1 ) ) goto TERMINATE;
    /* STEP3: update interval [al,au] */
    da = at - al;
    opt->_grad( opt, opt->_x, opt->_r, util );
    gt =-zVecInnerProdNC( opt->_g, opt->_d );
    if( *e1 > fl ){
      au = at;
    } else{
      al = at;
      fl =*e1;
      gl = gt;
      if( gt * da > 0 ) au = al;
    }
  }

 TERMINATE:
  opt->_grad( opt, opt->_x, opt->_r, util );
  zVecSubNC( opt->_r, opt->_g, opt->_q );
  zVecCopyNC( opt->_r, opt->_g );
  return zVecCopyNC( opt->_x, var );
}

/* line search by Dennis-Schnabel's backtracking */
#define Z_OPT_DM_STEP_MIN 0.001
zVec _zOptDMStepDS(zOptDM *opt, zVec var, void *util, double e0, double *e1)
{
  double a1, a2=1.0, da, v1, v2, s, df, c1, c2;

  if( zIsNan(e0) ){
    ZRUNWARN( ZM_WARN_OPT_BADINI );
    return NULL;
  }
  s = zVecInnerProdNC( opt->_g, opt->_d ); /* ascendant direction */
  df = opt->_df * s;
  /* initial test from full-length step */
  if( _zOptDMStepTest( opt, var, util, &a2, e0, df, e1 ) ) goto TERMINATE;
  /* quadratic approximation test */
  a1 = zMax( 0.5 * s / ( v2 = ( *e1 - e0 + a2*s ) ), Z_OPT_DM_STEP_MIN );
  if( _zOptDMStepTest( opt, var, util, &a1, e0, df, e1 ) ) goto TERMINATE;
  /* cubic approximation test */
  for( ; ; v2=v1 ){
    if( zIsTiny( ( da = a1 - a2 ) ) ) break;
    v1 = ( *e1 - e0 + s*a2 ) / ( a2*a2 );
    c1 = 3 * ( v1 - v2 ) / da;
    c2 = ( v2*a1 - v1*a2 ) / da;
    a2 = a1;
    a1 = zMax( (-c2+sqrt(c2*c2+c1*s)) / c1, Z_OPT_DM_STEP_MIN );
    if( _zOptDMStepTest( opt, var, util, &a1, e0, df, e1 ) ) break;
  }

 TERMINATE:
  opt->_grad( opt, opt->_x, opt->_r, util );
  zVecSubNC( opt->_r, opt->_g, opt->_q );
  zVecCopyNC( opt->_r, opt->_g );
  return zVecCopyNC( opt->_x, var );
}

/* line search by Nocedal's backtracking based on Wolfe's condition */
zVec _zOptDMStepNoc(zOptDM *opt, zVec var, void *util, double e0, double *e1)
{
  double a, s, df, cf;

  if( zIsNan(e0) ){
    ZRUNWARN( ZM_WARN_OPT_BADINI );
    return NULL;
  }
  s = zVecInnerProdNC( opt->_g, opt->_d );
  df = opt->_df * s;
  cf = opt->_cf * s;
  for( a=1.0; a>=Z_OPT_DM_STEP_MIN; ){
    if( !_zOptDMStepTest( opt, var, util, &a, e0, df, e1 ) ){
      a *= 0.5;
      continue;
    }
    opt->_grad( opt, opt->_x, opt->_r, util );
    if( a > Z_OPT_DM_STEP_MIN && zVecInnerProdNC( opt->_r, opt->_d ) > cf ){
      a *= 2.1;
      continue;
    }
    zVecSubNC( opt->_r, opt->_g, opt->_q );
    zVecCopyNC( opt->_r, opt->_g );
    return zVecCopyNC( opt->_x, var );
  }
  ZRUNWARN( ZM_WARN_OPT_BADINI );
  return NULL;
}

/* line search by Brent's strict minimization (particularly for CG) */
typedef struct{
  zOptDM *opt;
  zVec *var;
  void *util;
} _zOptDMStepData;

static double _zOptDMStepBrentEval(double x, void *util);
double _zOptDMStepBrentEval(double x, void *util)
{
  _zOptDMStepData *opt_data;

  opt_data = util;
  zVecCatNC( *opt_data->var, -x, opt_data->opt->_d, opt_data->opt->_x );
  return opt_data->opt->eval( opt_data->opt->_x, opt_data->util );
}
zVec _zOptDMStepBrent(zOptDM *opt, zVec var, void *util, double e0, double *e1)
{
  double a;
  _zOptDMStepData opt_data;

  if( zIsNan(e0) ){
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

/* *** update pseudo-Hessian matix *** */

zMat _zOptDMUpdateNone(zOptDM *opt, int count){
  if( ( count - 1 ) % zVecSizeNC(opt->_d) == 0 )
    opt->_b = 0; /* for CG restart */
  return opt->_h;
}
zMat _zOptDMUpdateDFP(zOptDM *opt, int count) /* Davidon-Fletcher-Powell */
{
  double k, l;

  l = 1.0 / ( k = zVecInnerProdNC( opt->_p, opt->_q ) );
  zMulMatVecNC( opt->_h, opt->_q, opt->_r );
  zMatCatDyadNC( opt->_h, (1+zVecInnerProdNC(opt->_r,opt->_q)/k)*l, opt->_p, opt->_p );
  zMatCatDyadNC( opt->_h,-l, opt->_p, opt->_r );
  return zMatCatDyadNC( opt->_h,-l, opt->_r, opt->_p );
}
zMat _zOptDMUpdateBFGS(zOptDM *opt, int count) /* Broyden-Fletcher-Goldfalb-Shanno */
{
  register int i, j;
  double k, l, pi, pj, ri, rj;

  zMulMatVecNC( opt->_h, opt->_q, opt->_r );
  k = 1.0 / zVecInnerProdNC( opt->_p, opt->_q );
  l = zVecInnerProdNC( opt->_q, opt->_r )*k + 1.0;
  for( i=0; i<zVecSizeNC(opt->_p); i++ ){
    pi = zVecElem( opt->_p, i );
    ri = zVecElem( opt->_r, i );
    for( j=0; j<zVecSizeNC(opt->_p); j++ ){
      pj = zVecElem( opt->_p, j );
      rj = zVecElem( opt->_r, j );
      zMatElem(opt->_h,i,j) += ( l*pi*pj - ri*pj - rj*pi ) * k;
    }
  }
  return opt->_h;
}

/* *** descent vector *** */

/* Steepest descent method */
zVec _zOptDMVecSD(zOptDM *opt, zVec var, zVec d, void *util)
{
  return zVecCopyNC( opt->_g, d );
}

/* Levenberg-Marquardt method */
zVec _zOptDMVecLM(zOptDM *opt, zVec var, zVec d, void *util)
{
  register int i;
  double m;

  opt->_hess( opt, var, opt->_h, util );
  m = zVecSqrNorm( opt->_g );
  for( i=0; i<zVecSizeNC(opt->_g); i++ )
    zMatElem(opt->_h,i,i) += m;
  return zLESolveGaussDST( opt->_h, opt->_g, opt->_d, opt->_idx, opt->_r );
}

/* variable metric method */
zVec _zOptDMVecVM(zOptDM *opt, zVec var, zVec d, void *util)
{
  return zMulMatVecNC( opt->_h, opt->_g, opt->_d );
}

/* conjugate gradient method */
zVec _zOptDMVecCG(zOptDM *opt, zVec var, zVec d, void *util)
{
  zVecMulDRC( opt->_d, opt->_b );
  return zVecAddNCDRC( opt->_d, opt->_g );
}

/* *** constructor & destructor *** */

#define Z_OPT_DM_SCALE ( 1.0e3 )
zOptDM *zOptDMCreate(zOptDM *opt, int dim, double scale, double (*eval)(zVec,void*), zVec (*grad)(zVec,zVec,void*), zMat (*hess)(zVec,zMat,void*))
{
  if( !( opt->eval = eval ) ){
    ZRUNERROR( ZM_ERR_OPT_NOEVAL );
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

void _zOptDMAssignStep(zOptDM *opt, const char *method)
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
  zVecFreeAO( 6, opt->_x, opt->_d, opt->_g, opt->_p, opt->_q, opt->_r );
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

int zOptDMSolve(zOptDM *opt, zVec var, void *util, double tol, int iter, double *eval)
{
  register int i;
  double _e_dummy, e;

  if( !eval ) eval = &_e_dummy;
  zMatIdent( opt->_h );
  ZITERINIT( iter );
  *eval = opt->eval( var, util );
  if( zVecIsTol( opt->_grad( opt, var, opt->_g, util ), tol ) )
    return 0; /* no need for iteration. */
  for( i=0; i<iter; i++, *eval=e ){
    opt->_vec( opt, var, opt->_d, util ); /* descent vector */
    if( !opt->_step( opt, var, util, *eval, &e ) ) break; /* line search */
    if( zVecIsTol(opt->_g,tol) || zVecIsTol(opt->_p,tol) )
      return i; /* succeed. */
    opt->_update_h( opt, i ); /* update quasi-hessian matrix */
  }
  ZITERWARN( iter );
  return -1;
}
