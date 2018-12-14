/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_rrt: RRT.
 */

#include <zm/zm_rrt.h>

/* return type of 'extend-tree' operation */
enum{
  ZRRT_EXT_OK, ZRRT_EXT_EPS, ZRRT_EXT_NONEED, ZRRT_EXT_UNACCEPT, ZRRT_EXT_ALLOCERR,
};

static void _zRRTListDestroy(zRRTList *list);
static bool _zRRTListAdd(zRRTList *list, zVec v, zRRTNode *parent);
static zRRTNode *_zRRTListNN(zRRTList *list, zVec vr, void *util, double (* dist)(zVec,zVec,void*));
static int _zRRTListExtend(zRRTList *list, zVec vr, double eps, void *util, double (* dist)(zVec,zVec,void*), zVec (* ext)(zVec,zVec,double,zVec,void*), bool (* chk)(zVec,void*));

static double _zRRTDistDefault(zVec v1, zVec v2, void *util);
static zVec _zRRTExtDefault(zVec vfrom, zVec vto, double eps, zVec v, void *util);

static void _zRRTPath(zRRT *rrt, zVecList *path);
static bool _zRRTPathShortcutCheck(zVec v1, zVec v2, zVec vm, double s1, double s2, double d, double eps, void *util, zVec (* ext)(zVec,zVec,double,zVec,void*), bool (* chk)(zVec,void*));

static int _zRRTListExtendEsc(zRRTList *list, zVec vr, double eps, void *util, double (* dist)(zVec,zVec,void*), zVec (* ext)(zVec,zVec,double,zVec,void*), bool (* chk)(zVec,void*), zVec vg);

/* (static)
 * _zRRTListDestroy
 * - destroy a node tree.
 */
void _zRRTListDestroy(zRRTList *list)
{
  zRRTListCell *cp;

  while( !zListIsEmpty( list ) ){
    zListDeleteTail( list, &cp );
    zVecFree( cp->data.v );
    zFree( cp );
  }
}

/* (static)
 * _zRRTListAdd
 * - add a node to a tree.
 */
bool _zRRTListAdd(zRRTList *list, zVec v, zRRTNode *parent)
{
  zRRTListCell *cell;

  if( !( cell = zAlloc( zRRTListCell, 1 ) ) ){
    ZALLOCERROR();
    return false;
  }
  cell->data.parent = parent;
  cell->data.v = v;
  zListInsertHead( list, cell );
  return true;
}

/* (static)
 * _zRRTListNN
 * - nearest neighbor search by a naive linear algorithm.
 */
zRRTNode *_zRRTListNN(zRRTList *list, zVec vr, void *util, double (* dist)(zVec,zVec,void*))
{
  zRRTListCell *cp, *nn = NULL;
  double d, min = HUGE_VAL;

  zListForEach( list, cp ){
    if( ( d = dist( vr, cp->data.v, util ) ) < min ){
      min = d;
      nn = cp;
    }
  }
  return &nn->data;
}

/* (static)
 * _zRRTListExtend
 * - extend-tree operation.
 * RETURN VALUE:
 * ZRRT_EXT_OK       success to extend tree.
 * ZRRT_EXT_EPS      success to directly connect to the node.
 * ZRRT_EXT_NONEED   the node is already on the tree.
 * ZRRT_EXT_UNACCEPT failure to extend tree due to the violation of constraints.
 * ZRRT_EXT_ALLOCERR failure to allocate new node.
 */
int _zRRTListExtend(zRRTList *list, zVec vr, double eps, void *util, double (* dist)(zVec,zVec,void*), zVec (* ext)(zVec,zVec,double,zVec,void*), bool (* chk)(zVec,void*))
{
  zRRTNode *nn;
  double norm;
  zVec ve;
  int ret;

  nn = _zRRTListNN( list, vr, util, dist );
  if( zIsTiny( ( norm = dist( vr, nn->v, util ) ) ) )
    return ZRRT_EXT_NONEED; /* no need to extend the tree */
  if( !( ve = zVecAlloc( zVecSizeNC(vr) ) ) ){
    ZALLOCERROR();
    return ZRRT_EXT_ALLOCERR;
  }
  if( zIsTol( norm, eps ) ){
    zVecCopyNC( vr, ve );
    ret = ZRRT_EXT_EPS;
  } else{
    ext( nn->v, vr, eps/norm, ve, util );
    ret = ZRRT_EXT_OK;
  }
  if( chk( ve, util ) ){
    zVecFree( ve );
    return ZRRT_EXT_UNACCEPT;
  }
  if( !_zRRTListAdd( list, ve, nn ) ){
    zVecFree( ve );
    return ZRRT_EXT_ALLOCERR;
  }
  return ret;
}

/* (static)
 * _zRRTDistDefault
 * - default distance function by a squareroot norm.
 */
double _zRRTDistDefault(zVec v1, zVec v2, void *util)
{
  return zVecDist( v1, v2 );
}

/* (static)
 * _zRRTExtDefault
 * - default extention function by a linear division.
 */
zVec _zRRTExtDefault(zVec vfrom, zVec vto, double eps, zVec v, void *util)
{
  zVecMulNC( vfrom, 1-eps, v );
  zVecCatNCDRC( v, eps, vto );
  return v;
}

/* zRRTInit
 * - initialize RRT solver.
 */
void zRRTInit(zRRT *rrt, zVec min, zVec max, double eps, double (* dist)(zVec,zVec,void*), zVec (* ext)(zVec,zVec,double,zVec,void*), bool (* chk)(zVec,void*))
{
  rrt->min = min;
  rrt->max = max;
  rrt->eps = eps;
  rrt->dist = dist ? dist : _zRRTDistDefault;
  rrt->ext = ext ? ext : _zRRTExtDefault;
  rrt->chk = chk;
  zListInit( &rrt->slist );
  zListInit( &rrt->glist );
}

/* zRRTDestroy
 * - destroy RRTs.
 */
void zRRTDestroy(zRRT *rrt)
{
  _zRRTListDestroy( &rrt->slist );
  _zRRTListDestroy( &rrt->glist );
}

/* (static)
 * _zRRTPath
 * - pick up the path from two RRTs.
 */
void _zRRTPath(zRRT *rrt, zVecList *path)
{
  zRRTNode *node;

  zListInit( path );
  for( node=&zListHead(&rrt->slist)->data; ; node=node->parent ){
    zVecListInsertTail( path, node->v, true );
    if( !node->parent ) break;
  }
  for( node=&zListHead(&rrt->glist)->data; ; node=node->parent ){
    zVecListInsertHead( path, node->v, true );
    if( !node->parent ) break;
  }
}

/* (static)
 * _zRRTPathShortcutCheck
 * - check if shortcut between two nodes on a path is available.
 */
bool _zRRTPathShortcutCheck(zVec v1, zVec v2, zVec vm, double s1, double s2, double d, double eps, void *util, zVec (* ext)(zVec,zVec,double,zVec,void*), bool (* chk)(zVec,void*))
{
  double s;

  if( d < eps ) return true;
  ext( v1, v2, ( s = 0.5 * ( s1 + s2 ) ), vm, util );
  if( chk( vm, util ) ) return false;
  d *= 0.5;
  return _zRRTPathShortcutCheck( v1, v2, vm, s1, s, d, eps, util, ext, chk ) &&
         _zRRTPathShortcutCheck( v1, v2, vm, s, s2, d, eps, util, ext, chk ) ?
    true : false;
}

/* zRRTConnect
 * - RRT-connect algorithm to find a path.
 */
bool zRRTConnect(zRRT *rrt, zVec start, zVec goal, int iter, void *util, zVecList *path)
{
  zVec vr;
  int i;
  bool ret = true;

  /* initialize trees from both start and goal */
  _zRRTListDestroy( &rrt->slist );
  _zRRTListDestroy( &rrt->glist );
  _zRRTListAdd( &rrt->slist, start, NULL );
  _zRRTListAdd( &rrt->glist, goal, NULL );
  vr = zVecAlloc( zVecSizeNC(start) );
  /* loop */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecRand( vr, rrt->min, rrt->max );
    _zRRTListExtend( &rrt->slist, vr, rrt->eps, util, rrt->dist, rrt->ext, rrt->chk );
    _zRRTListExtend( &rrt->glist, vr, rrt->eps, util, rrt->dist, rrt->ext, rrt->chk );
    if( _zRRTListExtend( &rrt->slist, zListHead(&rrt->glist)->data.v, rrt->eps, util, rrt->dist, rrt->ext, rrt->chk ) == ZRRT_EXT_EPS ||
        _zRRTListExtend( &rrt->glist, zListHead(&rrt->slist)->data.v, rrt->eps, util, rrt->dist, rrt->ext, rrt->chk ) == ZRRT_EXT_EPS )
      goto TERMINATE;
  }
  ZITERWARN( iter );
  ret = false;

 TERMINATE:
  _zRRTPath( rrt, path );
  zVecFree( vr );
  return ret;
}

/* zRRTPathShortcut
 * - shortcut a path.
 */
void zRRTPathShortcut(zRRT *rrt, void *util, zVecList *path)
{
  zVecListCell *cp1, *cp2, *tmp;
  zVec vm;

  if( zListNum(path) < 3 ) return;
  vm = zVecAlloc( zVecSizeNC( zListHead(path)->data ) );
  cp1 = zListTail(path);
  cp2 = zListCellNext( zListCellNext(cp1) );
  while( cp1!=zListHead(path) && cp2!=zListRoot(path) ){
    while( cp2!=zListRoot(path) ){
      if( !_zRRTPathShortcutCheck( cp1->data, cp2->data, vm, 0, 1, zVecDist(cp1->data,cp2->data), rrt->eps, util, rrt->ext, rrt->chk ) )
        break;
      zListDeletePrev( path, cp2, &tmp );
      zVecFree( tmp->data );
      zFree( tmp );
      cp2 = zListCellNext(cp2);
    }
    cp1 = zListCellPrev(cp2);
    cp2 = zListCellNext(cp2);
  }
  zVecFree( vm );
}

/* (static)
 * _zRRTListExtendEsc
 * - extend-tree operation for RRT-escapement.
 * RETURN VALUE:
 * ZRRT_EXT_OK       success to escape.
 * ZRRT_EXT_EPS      success to escape by directly connecting to the node.
 * ZRRT_EXT_NONEED   the node is already on the tree.
 * ZRRT_EXT_UNACCEPT yet to escape.
 * ZRRT_EXT_ALLOCERR failure to allocate new node.
 */
int _zRRTListExtendEsc(zRRTList *list, zVec vr, double eps, void *util, double (* dist)(zVec,zVec,void*), zVec (* ext)(zVec,zVec,double,zVec,void*), bool (* chk)(zVec,void*), zVec vg)
{
  zRRTNode *nn;
  double norm;
  zVec ve;

  nn = _zRRTListNN( list, vr, util, dist );
  if( zIsTiny( ( norm = dist( vr, nn->v, util ) ) ) )
    return ZRRT_EXT_NONEED; /* no need to extend the tree */
  if( !( ve = zVecAlloc( zVecSizeNC(vr) ) ) ){
    ZALLOCERROR();
    return ZRRT_EXT_ALLOCERR;
  }
  if( zIsTol( norm, eps ) ){
    zVecCopyNC( vr, ve );
  } else{
    ext( nn->v, vr, eps/norm, ve, util );
  }
  if( !_zRRTListAdd( list, ve, nn ) ){
    zVecFree( ve );
    return ZRRT_EXT_ALLOCERR;
  }
  if( chk( ve, util ) )
    return ZRRT_EXT_UNACCEPT;
  zVecCopy( ve, vg );
  return ZRRT_EXT_OK;
}

/* zRRTEsc
 * - RRT-Escapement algorithm to find a collision-free point
 *   proposed by Y. Shimizu in 2012.
 */
bool zRRTEsc(zRRT *rrt, zVec start, int iter, void *util, zVec goal)
{
  zVec vr;
  int i;
  bool ret = true;

  if( !rrt->chk( start, util ) ){
    zVecCopy( start, goal ); /* already escaped */
    return true;
  }
  /* initialize trees from both start and goal */
  _zRRTListDestroy( &rrt->slist );
  _zRRTListDestroy( &rrt->glist );
  _zRRTListAdd( &rrt->slist, start, NULL );
  vr = zVecAlloc( zVecSizeNC(start) );
  /* loop */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecRand( vr, rrt->min, rrt->max );
    if( _zRRTListExtendEsc( &rrt->slist, vr, rrt->eps, util, rrt->dist, rrt->ext, rrt->chk, goal ) != ZRRT_EXT_UNACCEPT )
      goto TERMINATE;
  }
  ZITERWARN( iter );
  ret = false;

 TERMINATE:
  zVecFree( vr );
  return ret;
}
