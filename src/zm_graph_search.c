/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_graph_search - graph class : shortest-path-finding methods.
 */

#include <zm/zm_graph.h>

/* graph node pool */

static double _zGraphEvalHDummy(void *d1, void *d2, void *util);

/* (static)
 * _zGraphEvalHDummy
 * - zero cost function (for Dijkstra method).
 */
double _zGraphEvalHDummy(void *d1, void *d2, void *util)
{
  return 0;
}

/* zGraphAStar
 * - find the shortest path of a graph by A* method.
 */
double zGraphAStar(zGraph *graph, void *start, void *goal, void *util, zGraphNodeList *path)
{
  zGraphNode *s=NULL, *g=NULL, *n;
  zGraphCell *gc;
  zGraphNodeListCell *nc, *cc;
  zGraphArcListCell *ac;
  double testcost, testcost_h;
  zGraphNodeList list_done, list_yet;

  if( !graph->h ) graph->h = _zGraphEvalHDummy;
  /* initialize */
  zListForEach( &graph->list, gc ){
    if( graph->equal( gc->data.data, goal ) ){
      g = &gc->data;
    } else{
      gc->data.val = gc->data.hval = HUGE_VAL;
      if( graph->equal( gc->data.data, start ) )
        s = &gc->data;
    }
  }
  if( !s || !g ){
    ZRUNERROR( "an invalid start or goal node specified" );
    return HUGE_VAL;
  }
  g->val = 0;
  g->hval = graph->h( g->data, s->data, util );

  zListInit( &list_done );
  zListInit( &list_yet );
  zGraphNodeListAdd( &list_yet, g );
  /* start updating */
  while( !zListIsEmpty( &list_yet ) ){
    cc = NULL;
    zListForEach( &list_yet, nc )
      if( !cc || nc->data->hval < cc->data->hval ) cc = nc;
    zListPurge( &list_yet, cc );
    zListInsertTail( &list_done, cc );
    if( cc->data == s ) break;
    zListForEach( &cc->data->arc, ac ){
      testcost = cc->data->val + ac->data.cost;
      testcost_h = testcost + graph->h( ac->data.node->data, s->data, util );
      if( ac->data.node->hval > testcost_h ){
        ac->data.node->val = testcost;
        ac->data.node->hval = testcost_h;
        ac->data.node->to = cc->data;
        zGraphNodeListAdd( &list_yet, ac->data.node );
      }
    }
  }
  /* path */
  zListInit( path );
  for( n=s; ; n=n->to ){
    zGraphNodeListAdd( path, n );
    if( n == g ) break;
  }
  return s->val;
}

/* zGraphDijkstra
 * - find the shortest path of a graph by Dijkstra's method.
 */
double zGraphDijkstra(zGraph *graph, void *start, void *goal, zGraphNodeList *path)
{
  graph->h = NULL;
  return zGraphAStar( graph, start, goal, NULL, path );
}
