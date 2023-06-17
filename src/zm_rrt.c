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
static bool _zRRTListAddNode(zRRTList *list, zVec v, double cost, zRRTNode *parent)
{
  zRRTListCell *cell;

  if( !( cell = zAlloc( zRRTListCell, 1 ) ) ){
    ZALLOCERROR();
    return false;
  }
  cell->data.v = v;
  cell->data.cost_edge = cost;
  cell->data.cost_total = ( cell->data.parent = parent ) ? parent->cost_total + cost : cost;
  zListInsertHead( list, cell );
  return true;
}

/* clone a node and add to a tree. */
static bool _zRRTListCloneNode(zRRTList *list, zVec v, double cost, zRRTNode *parent)
{
  zVec vc;

  if( !( vc = zVecClone( v ) ) ) return false;
  if( !_zRRTListAddNode( list, vc, cost, parent ) ){
    zVecFree( vc );
    return false;
  }
  return true;
}

/* nearest neighbor search by a naive linear algorithm. */
static double _zRRTListNN(zRRTList *list, zVec vr, void *util, double (* distance)(zVec,zVec,void*), zRRTNode **nn)
{
  zRRTListCell *cp;
  double d, dmin = HUGE_VAL;

  zListForEach( list, cp ){
    if( ( d = distance( vr, cp->data.v, util ) ) < dmin ){
      dmin = d;
      *nn = &cp->data;
    }
  }
  return dmin;
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
  ZRRT_EXT_NOP,        /* the node already on the tree */
  ZRRT_EXT_ADVANCED,   /* succeed to extend tree */
  ZRRT_EXT_REACHED,    /* reached to the node */
  ZRRT_EXT_TRAPPED,    /* failure to extend tree due to the violation of constraints */
  ZRRT_EXT_ALLOCERROR, /* failure to allocate new node */
};

/* create a new node and test if it is collision-free. */
static int _zRRTListExtendTest(zRRT *rrt, zRRTNode *nn, zVec vr, void *util, double dist, zVec ve, double *cost)
{
  int ret = ZRRT_EXT_ADVANCED;

  if( dist < rrt->eps ){
    zVecCopyNC( vr, ve );
    *cost = dist;
    ret = ZRRT_EXT_REACHED;
  } else{
    rrt->_extend( nn->v, vr, rrt->eps/dist, ve, util );
    *cost = rrt->eps;
  }
  return rrt->_check_collision( ve, util ) ? ZRRT_EXT_TRAPPED : ret;
}

/* extend-tree operation. */
static int _zRRTListExtend(zRRT *rrt, zRRTList *list, zVec vr, void *util)
{
  zRRTNode *nn = NULL;
  double dist, cost = HUGE_VAL;
  zVec ve;
  int ret;

  dist = _zRRTListNN( list, vr, util, rrt->_distance, &nn );
  if( zIsTiny( dist ) ) return ZRRT_EXT_NOP; /* no need to extend the tree */
  if( !( ve = zVecAlloc( zVecSizeNC(vr) ) ) ) return ZRRT_EXT_ALLOCERROR;
  if( ( ret = _zRRTListExtendTest( rrt, nn, vr, util, dist, ve, &cost ) ) == ZRRT_EXT_TRAPPED ){
    zVecFree( ve );
  } else{
    if( !_zRRTListAddNode( list, ve, cost, nn ) ){
      zVecFree( ve );
      return ZRRT_EXT_ALLOCERROR;
    }
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
static double _zRRTPath(zRRT *rrt, zVecList *path)
{
  zRRTNode *node;
  double cost = 0;

  zListInit( path );
  for( node=&zListHead(&rrt->slist)->data; ; node=node->parent ){
    cost += node->cost_edge;
    zVecListInsertTail( path, node->v );
    if( !node->parent ) break;
  }
  return cost;
}

/* skeleton of the RRT algorithm family. */
static bool _zRRTFindPath(zRRT *rrt, zVec start, int iter, void *util, zVecList *path, double *cost, int (* extendfunc)(zRRT*,zRRTList*,zVec,void*))
{
  zVec vr;
  int i, ext_ret;
  bool ret = true;

  /* initialize trees from both start and goal */
  _zRRTListDestroy( &rrt->slist );
  if( !_zRRTListCloneNode( &rrt->slist, start, 0, NULL ) ||
      !( vr = zVecAlloc( zVecSizeNC(start) ) ) ) return false;
  /* loop */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecRand( vr, rrt->min, rrt->max );
    if( ( ext_ret = extendfunc( rrt, &rrt->slist, vr, util ) ) == ZRRT_EXT_ALLOCERROR ) break;
    if( ext_ret == ZRRT_EXT_TRAPPED ) continue;
    if( rrt->_check_goal( zListHead(&rrt->slist)->data.v, util ) )
      goto TERMINATE;
  }
  ZITERWARN( iter );
  ret = false;

 TERMINATE:
  if( cost ) *cost = _zRRTPath( rrt, path );
  zVecFree( vr );
  return ret;
}

/* find a path based on the RRT algorithm. */
bool zRRTFindPath(zRRT *rrt, zVec start, int iter, void *util, zVecList *path, double *cost)
{
  return _zRRTFindPath( rrt, start, iter, util, path, cost, _zRRTListExtend );
}

/* node list for RRT* algorithm */
typedef struct{
  zRRTNode *node;
  double cost;
} zRRTNearNode;
zListClass( zRRTNearList, zRRTNearListCell, zRRTNearNode );

/* check if a pair of nodes is rewirable. */
bool _zRRTListCheckNodeRewirable(zRRT *rrt, zVec v1, zVec v2, double dist, void *util)
{
  zVec vm;
  int i, n;
  bool ret = true;

  if( !( vm = zVecAlloc( zVecSizeNC( v1 ) ) ) ) return false;
  n = dist / rrt->eps + 1;
  for( i=1; i<n; i++ ){
    rrt->_extend( v1, v2, (double)i/n, vm, util );
    if( rrt->_check_collision( vm, util ) ){
      ret = false;
      break;
    }
  }
  zVecFree( vm );
  return ret;
}

/* add a node to the near-list. */
static zRRTNearListCell *_zRRTNearListAdd(zRRTNearList *nearlist, zRRTNode *node, double cost)
{
  zRRTNearListCell *nearcp;

  if( !( nearcp = zAlloc( zRRTNearListCell, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  nearcp->data.node = node;
  nearcp->data.cost = cost;
  zListInsertHead( nearlist, nearcp );
  return nearcp;
}

/* add a new node to an RRT and rewire the near-nodes if possible. */
static bool _zRRTListAddAndRewireNode(zRRT *rrt, zRRTList *list, zVec ve, void *util, zRRTNode *nn, double dn, double r)
{
  zRRTListCell *cp;
  zRRTNearList nearlist;
  zRRTNearListCell *nearcp, *mincp = NULL;
  double d, cost_total_min = HUGE_VAL, cost_total;
  bool ret = true;

  zListInit( &nearlist );
  zListForEach( list, cp ){
    if( ( d = rrt->_distance( cp->data.v, ve, util ) ) < r ){
      if( !_zRRTListCheckNodeRewirable( rrt, ve, cp->data.v, d, util ) ) continue;
      if( !( nearcp = _zRRTNearListAdd( &nearlist, &cp->data, d ) ) )
        goto TERMINATE;
      if( ( cost_total = cp->data.cost_total + d ) < cost_total_min ){
        cost_total_min = cost_total;
        mincp = nearcp;
      }
    }
  }
  if( !mincp ){
    if( !( mincp = _zRRTNearListAdd( &nearlist, nn, dn ) ) ) goto TERMINATE;
  }
  if( !_zRRTListAddNode( list, ve, mincp->data.cost, mincp->data.node ) ){
    ret = false;
    goto TERMINATE;
  }
  zListPurge( &nearlist, mincp );
  zFree( mincp );
  zListForEach( &nearlist, nearcp ){
    if( ( cost_total = zListHead(list)->data.cost_total + nearcp->data.cost ) < nearcp->data.node->cost_total ){
      /* rewire */
      nearcp->data.node->parent = &zListHead(list)->data;
      nearcp->data.node->cost_edge = nearcp->data.cost;
      nearcp->data.node->cost_total = cost_total;
    }
  }

 TERMINATE:
  zListDestroy( zRRTNearListCell, &nearlist );
  return ret;
}

/* extend-tree operation for the RRT* algorithm. */
#define ZRRT_STAR_R_RATIO 100
static int _zRRTListExtendOpt(zRRT *rrt, zRRTList *list, zVec vr, void *util)
{
  zRRTNode *nn = NULL;
  double r, dist, dn = HUGE_VAL;
  zVec ve;
  int ret = ZRRT_EXT_ADVANCED;

  dist = _zRRTListNN( list, vr, util, rrt->_distance, &nn );
  if( zIsTiny( dist ) ) return ZRRT_EXT_NOP; /* no need to extend the tree */
  if( !( ve = zVecAlloc( zVecSizeNC(vr) ) ) ) return ZRRT_EXT_ALLOCERROR;
  if( ( ret = _zRRTListExtendTest( rrt, nn, vr, util, dist, ve, &dn ) ) == ZRRT_EXT_TRAPPED ){
    zVecFree( ve );
  } else{
    /* the following definition of near-radius is slightly different from the original paper
       because it should converge to eps as the number of nodes increases. */
    r = rrt->eps * ZRRT_STAR_R_RATIO *
      ( zListSize(list) == 1 ? 1 : pow( log(zListSize(list)) / zListSize(list), 1.0/zVecSizeNC(ve) ) ) + rrt->eps;
    if( !_zRRTListAddAndRewireNode( rrt, list, ve, util, nn, dn, r ) ){
      zVecFree( ve );
      return ZRRT_EXT_ALLOCERROR;
    }
  }
  return ret;
}

/* find the optimum path based on the RRT* algorithm. */
bool zRRTFindPathOpt(zRRT *rrt, zVec start, int iter, void *util, zVecList *path, double *cost)
{
  return _zRRTFindPath( rrt, start, iter, util, path, cost, _zRRTListExtendOpt );
}

/* pick up the path from two RRTs (for RRT-connect). */
static double _zRRTPathDual(zRRT *rrt, zVecList *path, void *util)
{
  zRRTNode *node;
  double cost = 0;

  zListInit( path );
  for( node=&zListHead(&rrt->slist)->data; ; node=node->parent ){
    cost += node->cost_edge;
    zVecListInsertTail( path, node->v );
    if( !node->parent ) break;
  }
  for( node=&zListHead(&rrt->glist)->data; ; node=node->parent ){
    cost += node->cost_edge;
    zVecListInsertHead( path, node->v );
    if( !node->parent ) break;
  }
  return cost + rrt->_distance( zListHead(&rrt->slist)->data.v, zListHead(&rrt->glist)->data.v, util );
}

/* find a path based on the RRT-connect algorithm. */
bool zRRTFindPathDual(zRRT *rrt, zVec start, zVec goal, int iter, void *util, zVecList *path, double *cost)
{
  zVec vr;
  int i, ext_ret;
  bool ret = true;
  zRRTList *list1, *list2;

  /* initialize trees from both start and goal */
  _zRRTListDestroy( &rrt->slist );
  _zRRTListDestroy( &rrt->glist );
  if( !_zRRTListCloneNode( &rrt->slist, start, 0, NULL ) ||
      !_zRRTListCloneNode( &rrt->glist, goal,  0, NULL ) ||
      !( vr = zVecAlloc( zVecSizeNC(start) ) ) ) return false;
  list1 = &rrt->slist;
  list2 = &rrt->glist;
  /* loop */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecRand( vr, rrt->min, rrt->max );
    if( ( ext_ret = _zRRTListExtend( rrt, list1, vr, util ) ) == ZRRT_EXT_ALLOCERROR ) break;
    if( ext_ret == ZRRT_EXT_TRAPPED ) continue;
    if( _zRRTListExtendGreedy( rrt, list2, zListHead(list1)->data.v, util ) == ZRRT_EXT_REACHED )
      goto TERMINATE;
    zSwap( zRRTList*, list1, list2 );
  }
  ZITERWARN( iter );
  ret = false;

 TERMINATE:
  if( cost ) *cost = _zRRTPathDual( rrt, path, util );
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
bool zRRTShortcutPath(zRRT *rrt, void *util, zVecList *path, double *cost)
{
  zVecListCell *cp1, *cp2, *tmp;
  zVec vm;
  double __cost;

  if( zListSize(path) < 3 ) return true;
  if( !( vm = zVecAlloc( zVecSizeNC( zListHead(path)->data ) ) ) ) return false;
  if( !cost ) cost = &__cost;
  *cost = 0;
  cp1 = zListTail(path);
  cp2 = zListCellNext( zListCellNext(cp1) );
  while( cp1 != zListHead(path) && cp2 != zListRoot(path) ){
    while( cp2 != zListRoot(path) ){
      if( !_zRRTShortcutPathCheck( rrt, cp1->data, cp2->data, vm, 0, 1, rrt->_distance(cp1->data,cp2->data,util), util ) )
        break;
      zListDeletePrev( path, cp2, &tmp );
      zVecFree( tmp->data );
      zFree( tmp );
      cp2 = zListCellNext(cp2);
    }
    *cost += rrt->_distance( cp1->data, zListCellPrev(cp2)->data, util );
    cp1 = zListCellPrev(cp2);
    cp2 = zListCellNext(cp2);
  }
  zVecFree( vm );
  return true;
}

/* extend-tree operation for RRT-escapement. */
static int _zRRTListExtendEscape(zRRT *rrt, zRRTList *list, zVec vr, void *util)
{
  zRRTNode *nn = NULL;
  double dist;
  zVec ve;
  int ret = ZRRT_EXT_TRAPPED;

  dist = _zRRTListNN( list, vr, util, rrt->_distance, &nn );
  if( zIsTiny( dist ) ) return ZRRT_EXT_NOP; /* no need to extend the tree */
  if( !( ve = zVecAlloc( zVecSizeNC(vr) ) ) ) return ZRRT_EXT_ALLOCERROR;
  if( dist < rrt->eps ){
    zVecCopyNC( vr, ve );
  } else{
    rrt->_extend( nn->v, vr, rrt->eps/dist, ve, util );
  }
  if( !rrt->_check_collision( ve, util ) ){
    ret = ZRRT_EXT_REACHED;
  }
  if( !_zRRTListAddNode( list, ve, 0, nn ) ){
    zVecFree( ve );
    return ZRRT_EXT_ALLOCERROR;
  }
  return ret;
}

/* RRT-escapement algorithm to find a collision-free point (proposed by Y. Shimizu in 2012). */
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
  if( !_zRRTListCloneNode( &rrt->slist, start, 0, NULL ) ) return false;
  if( !( vr = zVecAlloc( zVecSizeNC(start) ) ) ) return false;
  /* loop */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecRand( vr, rrt->min, rrt->max );
    if( ( ext_ret = _zRRTListExtendEscape( rrt, &rrt->slist, vr, util ) ) == ZRRT_EXT_REACHED )
      goto TERMINATE;
    if( ext_ret == ZRRT_EXT_ALLOCERROR ) break;
  }
  ZITERWARN( iter );
  ret = false;

 TERMINATE:
  zVecFree( vr );
  zVecCopy( zListHead(&rrt->slist)->data.v, goal );
  return ret;
}
