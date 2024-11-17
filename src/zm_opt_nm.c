/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_nm - optimization tools: Nelder-Mead's downhill simplex method.
 */
#include <zm/zm_opt.h>

/* data set for the Nelder-Mead (downhill simplex/polytope) method. */
typedef struct{
  int num;   /* number of vertices of simplex */
  zVec *e;   /* bases */
  zVec eval; /* evaluation values */
  zVec pin;  /* pin (COG of the bottom face of simplex) */
  zVec test; /* test stick */
  zIndex index;
} zOptNM;

/* create solver, allocating internal working space. */
static zOptNM *_zOptNMCreate(zOptNM *opt, int dim)
{
  int i;
  bool check = true;

  opt->num = dim+1;
  if( !( opt->eval = zVecAlloc(opt->num) ) ||
      !( opt->pin  = zVecAlloc(dim) ) ||
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

/* destroy internal workspace. */
static void _zOptNMDestroy(zOptNM *opt)
{
  int i;

  if( opt->e ){
    for( i=0; i<opt->num; i++ )
      zVecFree( opt->e[i] );
    zFree( opt->e );
  }
  zVecFreeAtOnce( 3, opt->eval, opt->pin, opt->test );
  zIndexFree( opt->index );
}

/* initialize internal workspace. */
static void _zOptNMInit(zOptNM *opt, zVec min, zVec max)
{
  int i;

  zVecMid( min, max, opt->e[0] );
  for( i=1; i<opt->num; i++ ){
    zVecCopyNC( opt->e[0], opt->e[i] );
    zVecElemNC(opt->e[i],i-1) += 1.0;
  }
}

/* evaluate all vertices. */
static void _zOptNMEvalAll(zOptNM *opt, double (* f)(zVec,void*), void *util)
{
  int i;

  for( i=0; i<opt->num; i++ )
    zVecSetElemNC( opt->eval, i, f( opt->e[i], util ) );
  zVecSort( opt->eval, opt->index );
}

/* reorder the working vertex according to evaluation. */
static void _zOptNMReord(zOptNM *opt)
{
  int i;

  for( i=0; i<opt->num; i++ )
    if( zVecElemNC(opt->eval,zIndexElemNC(opt->index,i)) > zVecElemNC(opt->eval,zIndexHead(opt->index)) ){
      zIndexMove( opt->index, zIndexSizeNC(opt->index)-1, i );
      break;
    }
}

/* compute pin (centroid of the bottom face of simplex), the point
 * to manipulate working vertex. */
static zVec _zOptNMPin(zOptNM *opt)
{
  int i;

  zVecZero( opt->pin );
  for( i=1; i<opt->num; i++ )
    zVecAddNCDRC( opt->pin, opt->e[zIndexElemNC(opt->index,i-1)] );
  return zVecDivNCDRC( opt->pin, opt->num-1 );
}

/* reflect the working vertex with respect to the pin. */
static zVec _zOptNMRefl(zOptNM *opt)
{
  zVecMulNC( opt->pin, 2, opt->test );
  return zVecSubNCDRC( opt->test, opt->e[zIndexHead(opt->index)] );
}

/* expand working vertex with respect to the pin. */
static zVec _zOptNMExp(zOptNM *opt)
{
  zVecMulNC( opt->e[zIndexHead(opt->index)], 2, opt->test );
  return zVecSubNCDRC( opt->test, opt->pin );
}

/* shrink the working vertex with respect to the pin. */
static zVec _zOptNMShrink(zOptNM *opt)
{
  zVecAddNC( opt->e[zIndexHead(opt->index)], opt->pin, opt->test );
  return zVecMulNCDRC( opt->test, 0.5 );
}

/* crunch the working vertex towards the best vertex. */
static void _zOptNMCrunch(zOptNM *opt)
{
  int i;

  for( i=0; i<opt->num; i++ )
    if( i != zIndexTail(opt->index) ){
      zVecAddNCDRC( opt->e[i], opt->e[zIndexTail(opt->index)] );
      zVecMulNCDRC( opt->e[i], 0.5 );
    }
}

/* check if the volume of simplex is tiny. */
static bool _zOptNMCheck(zOptNM *opt, double tol)
{
  int i;

  for( i=1; i<opt->num; i++ )
    if( !zIsTol( zVecDist(opt->e[0],opt->e[i]), tol ) )
      return false;
  return true;
}

/* amoeba-like iteration procedure of polytope method. */
static int _zOptNMTry(zOptNM *opt, zVec var, double (* f)(zVec,void*), void *util, int iter, double tol, double *eval)
{
  int i;
  double newval, bestval;

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    if( _zOptNMCheck( opt, tol ) ){
      zVecCopyNC( opt->e[zIndexTail(opt->index)], var );
      if( eval ) *eval = zVecElemNC( opt->eval, zIndexTail(opt->index) );
      return i; /* succeed. */
    }
    _zOptNMPin( opt );
    _zOptNMRefl( opt );
    bestval = zVecElemNC( opt->eval, zIndexTail(opt->index) );
    newval = f( opt->test, util );
    if( newval < bestval ){
      do{
        bestval = newval;
        zVecCopyNC( opt->test, opt->e[zIndexHead(opt->index)] );
        _zOptNMExp( opt );
        newval = f( opt->test, util );
      } while( newval < bestval );
      zVecSetElemNC( opt->eval, zIndexHead(opt->index), bestval );
      _zOptNMReord( opt );
      continue;
    }
    if( newval > zVecElemNC(opt->eval,zIndexNeck(opt->index)) ){
      _zOptNMShrink( opt );
      newval = f( opt->test, util );
      if( newval > zVecElemNC(opt->eval,zIndexNeck(opt->index)) ){
        _zOptNMCrunch( opt );
        _zOptNMEvalAll( opt, f, util );
        continue;
      }
    }
    zVecCopyNC( opt->test, opt->e[zIndexHead(opt->index)] );
    zVecSetElemNC( opt->eval, zIndexHead(opt->index), newval );
    _zOptNMReord( opt );
  }
  ZITERWARN( iter );
  return -1;
}

/* solve an optimization problem by Nelder-Mead method. */
int zOptSolveNM(double (* f)(zVec,void*), void *util, zVec min, zVec max, int iter, double tol, zVec ans, double *eval)
{
  zOptNM opt;
  int ret;

  if( !zVecSizeIsEqual( min, ans ) || !zVecSizeIsEqual( max, ans ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return -1;
  }
  if( !_zOptNMCreate( &opt, zVecSizeNC(ans) ) ) return -1;
  _zOptNMInit( &opt, min, max );
  _zOptNMEvalAll( &opt, f, util );
  ret = _zOptNMTry( &opt, ans, f, util, iter, tol, eval );
  _zOptNMDestroy( &opt );
  return ret;
}
