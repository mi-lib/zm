/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_nm - optimization tools: Nelder-Mead's downhill simplex method.
 */
#include <zm/zm_opt.h>

static void _zOptNMInit(zOptNM *opt, zVec var, void *util);
static void _zOptNMEvalAll(zOptNM *opt, void *util);
static void _zOptNMReord(zOptNM *opt);
static zVec _zOptNMPin(zOptNM *opt);

static zVec _zOptNMRefl(zOptNM *opt);
static zVec _zOptNMExp(zOptNM *opt);
static zVec _zOptNMShrink(zOptNM *opt);
static void _zOptNMCrunch(zOptNM *opt);

static bool _zOptNMCheck(zOptNM *opt, double tol);
static int _zOptNMTry(zOptNM *opt, zVec var, void *util, double tol, int iter, double *eval);

/* _zOptNMCreate
 * - create solver, allocating internal working space.
 */
zOptNM *zOptNMCreate(zOptNM *opt, int dim, double (*eval)(zVec,void*))
{
  register int i;
  bool check = true;

  opt->eval = eval;
  opt->num = dim+1;
  if( !( opt->f = zVecAlloc(opt->num) ) ||
      !( opt->pin = zVecAlloc(dim) ) ||
      !( opt->test = zVecAlloc(dim) ) ||
      !( opt->index = zIndexCreate(opt->num) ) ||
      !( opt->e = zAlloc( zVec, opt->num ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  for( i=0; i<opt->num; i++ )
    if( !( opt->e[i] = zVecAlloc(dim) ) ) check = false;
  if( !check ){
    ZALLOCERROR();
    return NULL;
  }
  return opt;
}

/* zOptNMDestroy
 * - destroy internal working space.
 */
void zOptNMDestroy(zOptNM *opt)
{
  register int i;

  if( opt->e ){
    for( i=0; i<opt->num; i++ )
      zVecFree( opt->e[i] );
    zFree( opt->e );
  }
  zVecFreeAO( 3, opt->f, opt->pin, opt->test );
  zIndexFree( opt->index );
  opt->eval = NULL;
}

/* (static)
 * _zOptNMInit
 * - initialize internal workspace.
 */
void _zOptNMInit(zOptNM *opt, zVec var, void *util)
{
  register int i;

  zVecCopyNC( var, opt->e[0] );
  for( i=1; i<opt->num; i++ ){
    zVecCopyNC( var, opt->e[i] );
    zVecElem(opt->e[i],i-1) += 1.0;
  }
}

/* (static)
 * _zOptNMEvalAll
 * - evaluate all vertices.
 */
void _zOptNMEvalAll(zOptNM *opt, void *util)
{
  register int i;

  for( i=0; i<opt->num; i++ )
    zVecSetElem( opt->f, i, opt->eval( opt->e[i], util ) );
  zVecSort( opt->f, opt->index );
}

/* (static)
 * _zOptNMReord
 * - reorder the working vertex according to evaluation.
 */
void _zOptNMReord(zOptNM *opt)
{
  register int i;

  for( i=0; i<opt->num; i++ )
    if( zVecElem(opt->f,zIndexElem(opt->index,i)) > zVecElem(opt->f,zIndexHead(opt->index)) ){
      zIndexMove( opt->index, zArrayNum(opt->index)-1, i );
      break;
    }
}

/* (static)
 * _zOptNMPin
 * - compute pin (COG of the bottom face of simplex), the point
 *   to manipulate working vertex.
 */
zVec _zOptNMPin(zOptNM *opt)
{
  register int i;

  zVecClear( opt->pin );
  for( i=1; i<opt->num; i++ )
    zVecAddNCDRC( opt->pin, opt->e[zIndexElem(opt->index,i-1)] );
  return zVecDivNCDRC( opt->pin, opt->num-1 );
}

/* (static)
 * _zOptNMRefl
 * - reflect the working vertex with respect to the pin.
 */
zVec _zOptNMRefl(zOptNM *opt)
{
  zVecMulNC( opt->pin, 2, opt->test );
  return zVecSubNCDRC( opt->test, opt->e[zIndexHead(opt->index)] );
}

/* (static)
 * _zOptNMExp
 * - expand working vertex with respect to the pin.
 */
zVec _zOptNMExp(zOptNM *opt)
{
  zVecMulNC( opt->e[zIndexHead(opt->index)], 2, opt->test );
  return zVecSubNCDRC( opt->test, opt->pin );
}

/* (static)
 * _zOptNMShrink
 * - shrink the working vertex with respect to the pin.
 */
zVec _zOptNMShrink(zOptNM *opt)
{
  zVecAddNC( opt->e[zIndexHead(opt->index)], opt->pin, opt->test );
  return zVecMulNCDRC( opt->test, 0.5 );
}

/* (static)
 * _zOptNMCrunch
 * - crunch the working vertex towards the best vertex.
 */
void _zOptNMCrunch(zOptNM *opt)
{
  register int i;

  for( i=0; i<opt->num; i++ )
    if( i != zIndexTail(opt->index) ){
      zVecAddNCDRC( opt->e[i], opt->e[zIndexTail(opt->index)] );
      zVecMulNCDRC( opt->e[i], 0.5 );
    }
}

/* (static)
 * _zOptNMCheck
 * - check if the volume of simplex is tiny.
 */
bool _zOptNMCheck(zOptNM *opt, double tol)
{
  register int i;

  for( i=1; i<opt->num; i++ )
    if( !zIsTol( zVecDist(opt->e[0],opt->e[i]), tol ) )
      return false;
  return true;
}

/* (static)
 * _zOptNMTry
 * - amoeba-like iteration procedure of polytope method.
 */
int _zOptNMTry(zOptNM *opt, zVec var, void *util, double tol, int iter, double *eval)
{
  register int i;
  double newval, bestval;

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    if( _zOptNMCheck( opt, tol ) ){
      zVecCopyNC( opt->e[zIndexTail(opt->index)], var );
      if( eval ) *eval = zVecElem( opt->f, zIndexTail(opt->index) );
      return i; /* succeed. */
    }
    _zOptNMPin( opt );
    _zOptNMRefl( opt );
    bestval = zVecElem( opt->f, zIndexTail(opt->index) );
    newval = opt->eval( opt->test, util );
    if( newval < bestval ){
      do{
        bestval = newval;
        zVecCopyNC( opt->test, opt->e[zIndexHead(opt->index)] );
        _zOptNMExp( opt );
        newval = opt->eval( opt->test, util );
      } while( newval < bestval );
      zVecSetElem( opt->f, zIndexHead(opt->index), bestval );
      _zOptNMReord( opt );
      continue;
    }
    if( newval > zVecElem(opt->f,zIndexNeck(opt->index)) ){
      _zOptNMShrink( opt );
      newval = opt->eval( opt->test, util );
      if( newval > zVecElem(opt->f,zIndexNeck(opt->index)) ){
        _zOptNMCrunch( opt );
        _zOptNMEvalAll( opt, util );
        continue;
      }
    }
    zVecCopyNC( opt->test, opt->e[zIndexHead(opt->index)] );
    zVecSetElem( opt->f, zIndexHead(opt->index), newval );
    _zOptNMReord( opt );
  }
  ZITERWARN( iter );
  return -1;
}

/* zOptNMSolve
 * - execute optimization.
 */
int zOptNMSolve(zOptNM *opt, zVec var, void *util, double tol, int iter, double *eval)
{
  _zOptNMInit( opt, var, util );
  _zOptNMEvalAll( opt, util );
  return _zOptNMTry( opt, var, util, tol, iter, eval );
}
