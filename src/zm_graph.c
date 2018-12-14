/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_graph - graph class.
 */

#include <zm/zm_graph.h>

/* ********************************************************** */
/* CLASS: zGraphNode
 * graph node class
 * ********************************************************** */

/* zGraphNodeInit
 * - initialize a graph node.
 */
void zGraphNodeInit(zGraphNode *node)
{
  node->val = node->hval = HUGE_VAL;
  zListInit( &node->arc );
  node->to = NULL;
}

/* zGraphNodeConnect
 * - connect two graph nodes.
 */
bool zGraphNodeConnect(zGraphNode *from, zGraphNode *to, double cost)
{
  zGraphArcListCell *ac;

  zListForEach( &from->arc, ac ){
    if( ac->data.node == to ){
      ZRUNWARN( ZM_WARN_GRAPH_DUPCON );
      ac->data.cost = cost;
      return true;
    }
  }
  if( !( ac = zAlloc( zGraphArcListCell, 1 ) ) ){
    ZALLOCERROR();
    return false;
  }
  ac->data.node = to;
  ac->data.cost = cost;
  zListInsertHead( &from->arc, ac );
  return true;
}

/* zGraphNodeBiconnect
 * - connect two graph nodes bidirectionally.
 */
bool zGraphNodeBiconnect(zGraphNode *n1, zGraphNode *n2, double cost)
{
  return zGraphNodeConnect( n1, n2, cost ) &&
         zGraphNodeConnect( n2, n1, cost ) ? true : false;
}

/* ********************************************************** */
/* CLASS: zGraph
 * graph class
 * ********************************************************** */

static void *_zGraphDupDummy(void *data);
static bool _zGraphEqualDummy(void *d1, void *d2);
static void _zGraphFWriteDummy(FILE *fp, void *data);
static void _zGraphDestroyDummy(void *data);
static double _zGraphEvalHDummy(void *d1, void *d2, void *util);

void *_zGraphDupDummy(void *data){ return data; }
bool _zGraphEqualDummy(void *d1, void *d2){ return true; }
void _zGraphFWriteDummy(FILE *fp, void *data){}
void _zGraphDestroyDummy(void *data){ zFree( data ); }
double _zGraphEvalHDummy(void *d1, void *d2, void *util){ return 0; }

/* zGraphInit
 * - initialize a graph.
 */
void zGraphInit(zGraph *graph)
{
  zListInit( &graph->list );
  graph->dup = _zGraphDupDummy;
  graph->equal = _zGraphEqualDummy;
  graph->fwrite = _zGraphFWriteDummy;
  graph->destroy = _zGraphDestroyDummy;
  graph->h = _zGraphEvalHDummy;
}

/* zGraphDestroy
 * - destroy a graph.
 */
void zGraphDestroy(zGraph *graph)
{
  zGraphCell *gc;

  while( !zListIsEmpty(&graph->list) ){
    zListDeleteHead( &graph->list, &gc );
    graph->destroy( gc->data.data );
    zListDestroy( zGraphArcListCell, &gc->data.arc );
    zFree( gc );
  }
}

/* zGraphAddNode
 * - add a node to a graph.
 */
bool zGraphAddNode(zGraph *graph, void *data)
{
  zGraphCell *gc;

  if( !zListIsEmpty( &graph->list ) )
    zListForEach( &graph->list, gc )
      if( graph->equal( gc->data.data, data ) ){
        ZRUNWARN( ZM_WARN_GRAPH_DUPSPC );
        return true;
      }
  if( !( gc = zAlloc( zGraphCell, 1 ) ) ){
    ZALLOCERROR();
    return false;
  }
  zGraphNodeInit( &gc->data );
  if( !( gc->data.data = graph->dup( data ) ) ){
    ZALLOCERROR();
    zFree( gc );
    return false;
  }
  zListInsertHead( &graph->list, gc );
  return true;
}

/* zGraphFindNode
 * - find a node from a graph.
 */
zGraphNode *zGraphFindNode(zGraph *graph, void *data)
{
  zGraphCell *gc;

  zListForEach( &graph->list, gc )
    if( graph->equal( gc->data.data, data ) )
      return &gc->data;
  return NULL;
}

/* zGraphConnect
 * - connect two graph nodes.
 */
bool zGraphConnect(zGraph *graph, void *from, void *to, double cost)
{
  zGraphNode *fromnode, *tonode;

  if( !( fromnode = zGraphFindNode( graph, from ) ) ||
      !( tonode = zGraphFindNode( graph, to ) ) )
    return false;
  return zGraphNodeConnect( fromnode, tonode, cost );
}

/* zGraphBiconnect
 * - connect two graph nodes bidirectionally.
 */
bool zGraphBiconnect(zGraph *graph, void *n1, void *n2, double cost)
{
  zGraphNode *node1, *node2;

  if( !( node1 = zGraphFindNode( graph, n1 ) ) ||
      !( node2 = zGraphFindNode( graph, n2 ) ) )
    return false;
  return zGraphNodeBiconnect( node1, node2, cost );
}

/* zGraphFWrite
 * - output information of a graph.
 */
void zGraphFWrite(FILE *fp, zGraph *graph)
{
  zGraphCell *gc;
  zGraphArcListCell *ac;

  fprintf( fp, "number of nodes = %d\n", zListNum(&graph->list) );
  zListForEach( &graph->list, gc ){
    fprintf( fp, "node" );
    graph->fwrite( fp, gc->data.data );
    fprintf( fp, " ... value=%g (H=%g)\n", gc->data.val, gc->data.hval );
    zListForEach( &gc->data.arc, ac ){
      fprintf( fp, " --> node" );
      graph->fwrite( fp, ac->data.node->data );
      fprintf( fp, " ... cost=%g\n", ac->data.cost );
    }
    if( gc->data.to ){
      fprintf( fp, " *min. cost providing to node" );
      graph->fwrite( fp, gc->data.to->data );
      fprintf( fp, "\n" );
    }
  }
}

/* ********************************************************** */
/* CLASS: zGraphPath
 * graph path class
 * ********************************************************** */

/* zGraphNodeListAdd
 * - add a node to a node list.
 */
bool zGraphNodeListAdd(zGraphNodeList *list, zGraphNode *node)
{
  zGraphNodeListCell *cp;

  zListForEach( list, cp )
    if( cp->data == node ) return true;
  if( !( cp = zAlloc( zGraphNodeListCell, 1 ) ) ){
    ZALLOCERROR();
    return false;
  }
  cp->data = node;
  zListInsertHead( list, cp );
  return true;
}
