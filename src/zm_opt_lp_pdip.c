/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lp_pdip - optimization tools:
 * linear programming by primal-dual interior-point method.
 */

#include <zm/zm_opt.h>

/* common matrix-vector set for primal-dual interior-point method */
typedef struct{
  /* pointers */
  zMat a;
  zVec b, c;
  zVec x;
  /* working memory */
  zIndex idx;
  zMat l, u, axza;
  zVec y, z, xz, tmp;
  zVec v1, v2, v3;
} _zLP_PDIP;

static _zLP_PDIP *_zLP_PDIPAlloc(_zLP_PDIP *dat, zMat a, zVec b, zVec c, zVec x);
static void _zLP_PDIPFree(_zLP_PDIP *dat);
static bool _zLP_PDIPIni(_zLP_PDIP *dat);
static bool _zLP_PDIPEqLUDecomp(_zLP_PDIP *dat);
static void _zLP_PDIPEqLU(_zLP_PDIP *dat, zVec v1, zVec v2, zVec v3, zVec dx, zVec dy, zVec dz);

/* (static)
 * _zLP_PDIPAlloc
 * - allocate working memory for PD-IP.
 */
_zLP_PDIP *_zLP_PDIPAlloc(_zLP_PDIP *dat, zMat a, zVec b, zVec c, zVec x)
{
  /* bind pointers */
  dat->a = a;
  dat->b = b;
  dat->c = c;
  dat->x = x;
  /* allocate working memory */
  dat->idx = zIndexCreate( zMatRowSizeNC(a) );
  dat->l = zMatAllocSqr( zMatRowSizeNC(a) );
  dat->u = zMatAllocSqr( zMatRowSizeNC(a) );
  dat->axza = zMatAllocSqr( zMatRowSizeNC(a) );
  dat->y = zVecAlloc( zVecSizeNC(b) );
  dat->z = zVecAlloc( zVecSizeNC(c) );
  dat->xz = zVecAlloc( zVecSizeNC(c) );
  dat->v1 = zVecAlloc( zVecSizeNC(b) );
  dat->v2 = zVecAlloc( zVecSizeNC(c) );
  dat->v3 = zVecAlloc( zVecSizeNC(c) );
  dat->tmp = zVecAlloc( zVecSizeNC(b) );

  if( !dat->idx || !dat->l || !dat->u || !dat->axza ||
      !dat->y || !dat->z || !dat->xz || !dat->tmp ||
      !dat->v1 || !dat->v2 || !dat->v3 ){
    ZALLOCERROR();
    _zLP_PDIPFree( dat );
    return NULL;
  }
  return dat;
}

/* (static)
 * _zLP_PDIPFree
 * - free working memory for PD-IP.
 */
void _zLP_PDIPFree(_zLP_PDIP *dat)
{
  /* unbind pointers */
  dat->a = NULL;
  dat->b = dat->c = dat->x = NULL;
  /* allocate working memory */
  zIndexFree( dat->idx );
  zMatFreeAO( 3, dat->l, dat->u, dat->axza );
  zVecFreeAO( 7,
    dat->y, dat->z, dat->xz, dat->tmp, dat->v1, dat->v2, dat->v3 );
}

/* (static)
 * _zLP_PDIPIni
 * - initial interior point for PD-IP.
 */
bool _zLP_PDIPIni(_zLP_PDIP *dat)
{
  double m;

  if( !zLESolveNormMin( dat->a, dat->b, NULL, dat->x ) ||
      !zLESolveMP( dat->a, dat->b, NULL, NULL, dat->x ) ){
    ZRUNERROR( ZM_ERR_OPT_INI );
    return false;
  }
  zVecSetAll( dat->y, 1.0 );
  zMulMatTVecNC( dat->a, dat->y, dat->z );
  zVecSubNC( dat->c, dat->z, dat->z );
  if( ( m = zVecMin(dat->x,NULL) ) <= zTOL )
    zVecShift( dat->x, -m+zTOL );
  if( ( m = zVecMin(dat->z,NULL) ) <= zTOL )
    zVecShift( dat->z, -m+zTOL );
  return true;
}

/* (static)
 * _zLP_PDIPEqLUDecomp
 * - LU decomposition of gradient matrix for PD-IP.
 */
bool _zLP_PDIPEqLUDecomp(_zLP_PDIP *dat)
{
  zVecDemNC( dat->x, dat->z, dat->xz ); /* xz = x / z */
  zMatQuadNC( dat->a, dat->xz, dat->axza ); /* Ax/zA^T */
  if( zLUDecomp( dat->axza, dat->l, dat->u, dat->idx ) == zMatRowSizeNC(dat->a) )
    return true;
  ZRUNERROR( ZM_ERR_LE_SINGULAR );
  return false;
}

/* (static)
 * _zLP_PDIPEqLU
 * - solve gradient equation by LU factorization for PD-IP.
 */
void _zLP_PDIPEqLU(_zLP_PDIP *dat, zVec v1, zVec v2, zVec v3, zVec dx, zVec dy, zVec dz)
{
  /* dy (dz for temporary working space) */
  v2 ? zVecAmpNC( v2, dat->x, dz ) : zVecClear( dz );
  zVecSubNCDRC( dz, v3 );
  zVecDemNCDRC( dz, dat->z );
  zMulMatVecNC( dat->a, dz, dat->tmp );
  if( v1 ) zVecAddNCDRC( dat->tmp, v1 );
  zVecRevNCDRC( dat->tmp );
  zLESolve_L_U( dat->l, dat->u, dat->tmp, dy, dat->idx );
  /* dz */
  zMulMatTVecNC( dat->a, dy, dz );
  if( v2 ) zVecAddNCDRC( dz, v2 );
  zVecRevNCDRC( dz );
  /* dx */
  zVecAmpNC( dat->x, dz, dx );
  zVecAddNCDRC( dx, v3 );
  zVecDemNCDRC( dx, dat->z );
  zVecRevNCDRC( dx );
}

/* Mehrotra's predictor-corrector method */
typedef struct{
  zVec dx, dy, dz;
  zVec dx2, dy2, dz2;
} _zLP_PDIP_PC;

static _zLP_PDIP_PC *_zLP_PDIP_PCAlloc(_zLP_PDIP_PC *pc, zVec b, zVec c);
static void _zLP_PDIP_PCFree(_zLP_PDIP_PC *pc);
static double _zLP_PDIP_PCErr(_zLP_PDIP *dat);
static double _zLP_PDIP_PCStep(zVec x, zVec dx);
static void _zLP_PDIP_PCCorrect(_zLP_PDIP *dat, _zLP_PDIP_PC *pc, double ap, double ad, double e);

/* (static)
 * _zLP_PDIP_PCAlloc
 * - allocate working memory for PD-IP-PC.
 */
_zLP_PDIP_PC *_zLP_PDIP_PCAlloc(_zLP_PDIP_PC *pc, zVec b, zVec c)
{
  pc->dx = zVecAlloc( zVecSizeNC(c) );
  pc->dy = zVecAlloc( zVecSizeNC(b) );
  pc->dz = zVecAlloc( zVecSizeNC(c) );
  pc->dx2 = zVecAlloc( zVecSizeNC(c) );
  pc->dy2 = zVecAlloc( zVecSizeNC(b) );
  pc->dz2 = zVecAlloc( zVecSizeNC(c) );
  if( !pc->dx || !pc->dy || !pc->dz || !pc->dx2 || !pc->dy2 || !pc->dz2 ){
    ZALLOCERROR();
    _zLP_PDIP_PCFree( pc );
    return NULL;
  }
  return pc;
}

/* (static)
 * _zLP_PDIP_PCFree
 * - free working memory for PD-IP-PC.
 */
void _zLP_PDIP_PCFree(_zLP_PDIP_PC *pc)
{
  zVecFreeAO( 6, pc->dx, pc->dy, pc->dz, pc->dx2, pc->dy2, pc->dz2 );
}

/* (static)
 * _zLP_PDIP_PCErr
 * - error vector of PD-IP-PC.
 */
double _zLP_PDIP_PCErr(_zLP_PDIP *dat)
{
  zMulMatVecNC( dat->a, dat->x, dat->v1 );
  zVecSubNCDRC( dat->v1, dat->b );       /* v1 = Ax - b */
  zMulMatTVecNC( dat->a, dat->y, dat->v2 );
  zVecAddNCDRC( dat->v2, dat->z );
  zVecSubNCDRC( dat->v2, dat->c );       /* v2 = A^T y + z - c */
  zVecAmpNC( dat->x, dat->z, dat->v3 );  /* v3 = x * z */
  return zVecSum( dat->v3 );             /* x^T z */
}

/* (static)
 * _zLP_PDIP_PCStep
 * - updating step of PD-IP-PC.
 */
double _zLP_PDIP_PCStep(zVec x, zVec dx)
{
  register int i;
  bool max_ok = false, min_ok = false;
  double max, min, val, d;

  for( max=min=0, i=0; i<zVecSizeNC(dx); i++ ){
    d = zVecElem(dx,i);
    if( zIsTiny(d) && zVecElem(x,i) < -zTOL ){
      ZRUNWARN( ZM_ERR_OPT_UNSOLVE );
      return 0;
    }
    val = -zVecElem(x,i)/d;
    if( d < -zTOL ){
      if( !max_ok || val < max ){
        max = val;
        max_ok = true;
      }
    } else
    if( d > zTOL ){
      if( !min_ok || val > min ){
        min = val;
        min_ok = true;
      }
    }
  }
  if( min > max ){
    ZRUNWARN( ZM_ERR_OPT_UNSOLVE );
    return 0;
  }
  return max;
}

/* (static)
 * _zLP_PDIP_PCCorrect
 * - corrector magnitude of PD-IP-PC.
 */
void _zLP_PDIP_PCCorrect(_zLP_PDIP *dat, _zLP_PDIP_PC *pc, double ap, double ad, double e)
{
  double g;

  zVecAmpNC( pc->dx, pc->dz, dat->v3 );
  zVecCatNC( dat->x, ap, pc->dx, pc->dx2 );
  zVecCatNC( dat->z, ad, pc->dz, pc->dz2 );
  g = zVecInnerProd(pc->dx2,pc->dz2) / e;
  zVecShift( dat->v3, -zCube(g) * e / zVecSizeNC(dat->x) );
}

/* zLPSolvePDIP_PC
 * - linear programming solver based on primal-dual interior-point
 *   method with Mehrotra s predictor-corrector.
 */
#define ZM_LPSOLVE_PDIP_PCSTEP_MUL 0.99
bool zLPSolvePDIP_PC(zMat a, zVec b, zVec c, zVec x, double *cost)
{
  _zLP_PDIP dat;
  _zLP_PDIP_PC pc;
  double e, ap, ad;
  bool ret = false;
  int iter = 0;
  register int i;

  if( !_zLP_PDIPAlloc( &dat, a, b, c, x ) ) goto TERMINATE2;
  if( !_zLP_PDIP_PCAlloc( &pc, b, c ) ) goto TERMINATE1;
  _zLP_PDIPIni( &dat );
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    e = _zLP_PDIP_PCErr( &dat );
    if( zVecIsTiny(dat.v1) && zVecIsTiny(dat.v2) && zIsTiny(e) ) goto TERMINATE0;
    if( !_zLP_PDIPEqLUDecomp( &dat ) ) goto TERMINATE1;
    /* predictor */
    _zLP_PDIPEqLU( &dat, dat.v1, dat.v2, dat.v3, pc.dx, pc.dy, pc.dz );
    ap = _zLP_PDIP_PCStep( dat.x, pc.dx );
    ad = _zLP_PDIP_PCStep( dat.z, pc.dz );
    if( zIsTiny(ap) && zIsTiny(ad) ) goto TERMINATE0;
    _zLP_PDIP_PCCorrect( &dat, &pc, ap, ad, e );
    /* corrector */
    _zLP_PDIPEqLU( &dat, NULL, NULL, dat.v3, pc.dx2, pc.dy2, pc.dz2 );
    zVecAddNCDRC( pc.dx, pc.dx2 );
    zVecAddNCDRC( pc.dy, pc.dy2 );
    zVecAddNCDRC( pc.dz, pc.dz2 );
    if( ( ap = ZM_LPSOLVE_PDIP_PCSTEP_MUL * _zLP_PDIP_PCStep( dat.x, pc.dx ) ) > 1 ) ap = 1;
    if( ( ad = ZM_LPSOLVE_PDIP_PCSTEP_MUL * _zLP_PDIP_PCStep( dat.z, pc.dz ) ) > 1 ) ad = 1;
    if( zIsTiny(ap) && zIsTiny(ad) ) goto TERMINATE0;
    zVecCatNCDRC( dat.x, ap, pc.dx );
    zVecCatNCDRC( dat.y, ad, pc.dy );
    zVecCatNCDRC( dat.z, ad, pc.dz );
  }
  ZITERWARN( iter );
 TERMINATE0:
  zVecTouchup( x );
  ret = true;
 TERMINATE1:
  _zLP_PDIP_PCFree( &pc );
 TERMINATE2:
  _zLP_PDIPFree( &dat );
  if( cost ) *cost = zVecInnerProd( x, c );
  return ret;
}
