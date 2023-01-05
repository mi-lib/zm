/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_rrt: Rapidly-explored Random Tree algorithm and its family.
 */

#include <zm/zm_rrt.h>

/* destroy a tree. */
static void _zRRTListDestroy(zRRTList *list)
{
  zRRTListCell *cp;

  while( !zListIsEmpty( list ) ){
    zListDeleteHead( list, &cp );
    zVecFree( cp->data.v );
    zFree( cp );
  }
}

/* add a node to a tree. */
static bool _zRRTListAddNode(zRRTList *list, zVec v, zRRTNode *parent)
{
  zRRTListCell *cell;

  if( !( cell = zAlloc( zRRTListCell, 1 ) ) ){
    ZALLOCERROR();
    return false;
  }
  cell->data.v = v;
  cell->data.parent = parent;
  zListInsertHead( list, cell );
  return true;
}

/* clone a node and add to a tree. */
static bool _zRRTListCloneNode(zRRTList *list, zVec v, zRRTNode *parent)
{
  zVec vc;

  if( !( vc = zVecClone( v ) ) ) return false;
  if( !_zRRTListAddNode( list, vc, parent ) ){
    zVecFree( vc );
    return false;
  }
  return true;
}

/* nearest neighbor search by a naive linear algorithm. */
static zRRTNode *_zRRTListNN(zRRTList *list, zVec vr, void *util, double (* distance)(zVec,zVec,void*), double *dmin)
{
  zRRTListCell *cp, *nn = NULL;
  double d;

  *dmin = HUGE_VAL;
  zListForEach( list, cp ){
    if( ( d = distance( vr, cp->data.v, util ) ) < *dmin ){
      *dmin = d;
      nn = cp;
    }
  }
  return &nn->data;
}

/* default distance function by a squareroot norm. */
static double _zRRTDistanceFuncDefault(zVec v1, zVec v2, void *dummy)
{
  return zVecDist( v1, v2 );
}

/* default extention function by a linear division. */
static zVec _zRRTExtendFuncDefault(zVec vfrom, zVec vto, double eps, zVec v, void *dummy)
{
  zVecMulNC( vfrom, 1-eps, v );
  zVecCatNCDRC( v, eps, vto );
  return v;
}

/* default (dummy) collision-check function. */
static bool _zRRTCheckCollisionFuncDefault(zVec v, void *dummy){ return true; }

/* default (dummy) goal-check function. */
static bool _zRRTCheckGoalFuncDefault(zVec v, void *dummy){ return true; }

/* return type of extend-tree operation */
enum{
  ZRRT_EXT_NOP,      /* the node already on the tree */
  ZRRT_EXT_ADVANCED, /* succeed to extend tree */
  ZRRT_EXT_REACHED,  /* reached to the node */
  ZRRT_EXT_TRAPPED,  /* failure to extend tree due to the violation of constraints */
  ZRRT_EXT_ALLOCERR, /* failure to allocate new node */
};

/* extend-tree operation. */
static int _zRRTListExtend(zRRT *rrt, zRRTList *list, zVec vr, void *util)
{
  zRRTNode *nn;
  double norm;
  zVec ve;
  int ret = ZRRT_EXT_ADVANCED;

  nn = _zRRTListNN( list, vr, util, rrt->_distance, &norm );
  if( zIsTiny( norm ) ) return ZRRT_EXT_NOP; /* no need to extend the tree */
  if( !( ve = zVecAlloc( zVecSizeNC(vr) ) ) ){
    ZALLOCERROR();
    return ZRRT_EXT_ALLOCERR;
  }
  if( zIsTol( norm, rrt->eps ) ){
    zVecCopyNC( vr, ve );
    ret = ZRRT_EXT_REACHED;
  } else{
    rrt->_extend( nn->v, vr, rrt->eps/norm, ve, util );
  }
  if( rrt->_check_collision( ve, util ) ){
    zVecFree( ve );
    return ZRRT_EXT_TRAPPED;
  }
  if( !_zRRTListAddNode( list, ve, nn ) ){
    zVecFree( ve );
    return ZRRT_EXT_ALLOCERR;
  }
  return ret;
}

/* greedy-extend-tree operation for RRT-connect. */
static int _zRRTListExtendGreedy(zRRT *rrt, zRRTList *list, zVec vr, void *util)
{
  int ret;

  do{
    ret = _zRRTListExtend( rrt, list, vr, util );
  } while( ret == ZRRT_EXT_ADVANCED );
  return ret;
}

/* initialize RRTs. */
void zRRTInit(zRRT *rrt, zVec min, zVec max, double eps, double (* distance)(zVec,zVec,void*), zVec (* extend)(zVec,zVec,double,zVec,void*), bool (* chk_collision)(zVec,void*), bool (* chk_goal)(zVec,void*))
{
  rrt->min = min;
  rrt->max = max;
  rrt->eps = eps;
  zRRTSetDistanceFunc( rrt, distance ? distance : _zRRTDistanceFuncDefault );
  zRRTSetExtendFunc( rrt, extend ? extend : _zRRTExtendFuncDefault );
  zRRTSetCheckCollisionFunc( rrt, chk_collision ? chk_collision : _zRRTCheckCollisionFuncDefault );
  zRRTSetCheckGoalFunc( rrt, chk_goal ? chk_goal : _zRRTCheckGoalFuncDefault );
  zListInit( &rrt->slist );
  zListInit( &rrt->glist );
}

/* destroy RRTs. */
void zRRTDestroy(zRRT *rrt)
{
  _zRRTListDestroy( &rrt->slist );
  _zRRTListDestroy( &rrt->glist );
}

/* pick up the path from an RRT. */
static void _zRRTPath(zRRT *rrt, zVecList *path)
{
  zRRTNode *node;

  zListInit( path );
  for( node=&zListHead(&rrt->slist)->data; ; node=node->parent ){
    zVecListInsertTail( path, node->v );
    if( !node->parent ) break;
  }
}

/* find a path based on the RRT algorithm. */
bool zRRTFindPath(zRRT *rrt, zVec start, int iter, void *util, zVecList *path)
{
  zVec vr;
  int i, ext_ret;
  bool ret = true;

  /* initialize trees from both start and goal */
  _zRRTListDestroy( &rrt->slist );
  _zRRTListCloneNode( &rrt->slist, start, NULL );
  vr = zVecAlloc( zVecSizeNC(start) );
  /* loop */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecRand( vr, rrt->min, rrt->max );
    if( ( ext_ret = _zRRTListExtend( rrt, &rrt->slist, vr, util ) ) != ZRRT_EXT_TRAPPED )
      continue;
    if( ext_ret == ZRRT_EXT_ALLOCERR ) break;
    if( rrt->_check_goal( zListHead(&rrt->slist)->data.v, util ) )
      goto TERMINATE;
  }
  ZITERWARN( iter );
  ret = false;

 TERMINATE:
  _zRRTPath( rrt, path );
  zVecFree( vr );
  return ret;
}

/* pick up the path from two RRTs (for RRT-connect). */
static void _zRRTPathDual(zRRT *rrt, zVecList *path)
{
  zRRTNode *node;

  zListInit( path );
  for( node=&zListHead(&rrt->slist)->data; ; node=node->parent ){
    zVecListInsertTail( path, node->v );
    if( !node->parent ) break;
  }
  for( node=&zListHead(&rrt->glist)->data; ; node=node->parent ){
    zVecListInsertHead( path, node->v );
    if( !node->parent ) break;
  }
}

/* find a path based on the RRT-connect algorithm. */
bool zRRTFindPathDual(zRRT *rrt, zVec start, zVec goal, int iter, void *util, zVecList *path)
{
  zVec vr;
  int i, ext_ret;
  bool ret = true;
  zRRTList *list1, *list2;

  /* initialize trees from both start and goal */
  _zRRTListDestroy( &rrt->slist );
  _zRRTListDestroy( &rrt->glist );
  _zRRTListCloneNode( &rrt->slist, start, NULL );
  _zRRTListCloneNode( &rrt->glist, goal, NULL );
  vr = zVecAlloc( zVecSizeNC(start) );
  list1 = &rrt->slist;
  list2 = &rrt->glist;
  /* loop */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecRand( vr, rrt->min, rrt->max );
    if( ( ext_ret = _zRRTListExtend( rrt, list1, vr, util ) ) != ZRRT_EXT_TRAPPED )
      continue;
    if( ext_ret == ZRRT_EXT_ALLOCERR ) break;
    if( _zRRTListExtendGreedy( rrt, list2, zListHead(list1)->data.v, util ) == ZRRT_EXT_REACHED )
      goto TERMINATE;
    zSwap( zRRTList*, list1, list2 );
  }
  ZITERWARN( iter );
  ret = false;

 TERMINATE:
  _zRRTPathDual( rrt, path );
  zVecFree( vr );
  return ret;
}

/* check if shortcut between two nodes on a path is collision-free. */
static bool _zRRTShortcutPathCheck(zRRT *rrt, zVec v1, zVec v2, zVec vm, double s1, double s2, double d, void *util)
{
  double s;

  if( d < rrt->eps ) return true;
  rrt->_extend( v1, v2, ( s = 0.5 * ( s1 + s2 ) ), vm, util );
  if( rrt->_check_collision( vm, util ) ) return false;
  d *= 0.5;
  return _zRRTShortcutPathCheck( rrt, v1, v2, vm, s1, s, d, util ) &&
         _zRRTShortcutPathCheck( rrt, v1, v2, vm, s, s2, d, util ) ?
    true : false;
}

/* a postprocess for RRT family to shortcut a path. */
void zRRTShortcutPath(zRRT *rrt, void *util, zVecList *path)
{
  zVecListCell *cp1, *cp2, *tmp;
  zVec vm;

  if( zListSize(path) < 3 ) return;
  vm = zVecAlloc( zVecSizeNC( zListHead(path)->data ) );
  cp1 = zListTail(path);
  cp2 = zListCellNext( zListCellNext(cp1) );
  while( cp1!=zListHead(path) && cp2!=zListRoot(path) ){
    while( cp2!=zListRoot(path) ){
      if( !_zRRTShortcutPathCheck( rrt, cp1->data, cp2->data, vm, 0, 1, rrt->_distance(cp1->data,cp2->data,util), util ) )
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

/* extend-tree operation for RRT-escapement. */
static int _zRRTListExtendEscape(zRRT *rrt, zRRTList *list, zVec vr, void *util)
{
  zRRTNode *nn;
  double norm;
  zVec ve;
  int ret = ZRRT_EXT_TRAPPED;

  nn = _zRRTListNN( list, vr, util, rrt->_distance, &norm );
  if( zIsTiny( norm ) ) return ZRRT_EXT_NOP; /* no need to extend the tree */
  if( !( ve = zVecAlloc( zVecSizeNC(vr) ) ) ){
    ZALLOCERROR();
    return ZRRT_EXT_ALLOCERR;
  }
  if( zIsTol( norm, rrt->eps ) ){
    zVecCopyNC( vr, ve );
  } else{
    rrt->_extend( nn->v, vr, rrt->eps/norm, ve, util );
  }
  if( !rrt->_check_collision( ve, util ) ){
    ret = ZRRT_EXT_REACHED;
  }
  if( !_zRRTListAddNode( list, ve, nn ) ){
    zVecFree( ve );
    return ZRRT_EXT_ALLOCERR;
  }
  return ret;
}

/* RRT-escapement algorithm to find a collision-free point proposed by Y. Shimizu in 2012. */
bool zRRTEscape(zRRT *rrt, zVec start, int iter, void *util, zVec goal)
{
  zVec vr;
  int i, ext_ret;
  bool ret = true;

  if( !rrt->_check_collision( start, util ) ){
    zVecCopy( start, goal ); /* already escaped */
    return true;
  }
  /* initialize trees from both start and goal */
  _zRRTListDestroy( &rrt->slist );
  _zRRTListCloneNode( &rrt->slist, start, NULL );
  vr = zVecAlloc( zVecSizeNC(start) );
  /* loop */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecRand( vr, rrt->min, rrt->max );
    if( ( ext_ret = _zRRTListExtendEscape( rrt, &rrt->slist, vr, util ) ) == ZRRT_EXT_REACHED )
      goto TERMINATE;
    if( ext_ret == ZRRT_EXT_ALLOCERR ) break;
  }
  ZITERWARN( iter );
  ret = false;

 TERMINATE:
  zVecFree( vr );
  zVecCopy( zListHead(&rrt->slist)->data.v, goal );
  return ret;
}
