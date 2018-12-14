#include <zm/zm_graph.h>

#define TEST 0

typedef struct{
  int id;
  double x;
  double y;
} data_t;

void *test_dup(void *data)
{
  data_t *d;

  if( !( d = zAlloc( data_t, 1 ) ) ) return NULL;
  d->id = ((data_t *)data)->id;
  d->x  = ((data_t *)data)->x;
  d->y  = ((data_t *)data)->y;
  return d;
}

bool test_equal(void *d1, void *d2)
{
  return ((data_t *)d1)->id == ((data_t *)d2)->id;
}

void test_fwrite(FILE *fp, void *data)
{
  fprintf( fp, "[%d](%g,%g)", ((data_t *)data)->id, ((data_t *)data)->x, ((data_t *)data)->y );
}

void test_destroy(void *data)
{
  zFree( data );
}

double test_h(void *d1, void *d2, void *dummy)
{
  data_t *v1, *v2;

  v1 = d1;
  v2 = d2;
  return fabs( v1->x - v2->x ) + fabs( v1->y - v2->y );
}

bool test_connect(zGraph *graph, int i1, int i2, double cost)
{
  data_t n1, n2;

  n1.id = i1;
  n2.id = i2;
  return zGraphBiconnect( graph, &n1, &n2, cost );
}


#if TEST == 1
/* TEST1 maze:
 4 2 ooo G
   o o o o
 3 1o3 o o
   o o o o
 2 o 4o5o6
   o o   o
 1 S ooooo
   1 2 3 4
 */
static data_t maze[]={
  { 0, 1.0, 1.0 }, /* 0:S */
  { 1, 1.0, 3.0 }, /* 1 */
  { 2, 1.0, 4.0 }, /* 2 */
  { 3, 2.0, 3.0 }, /* 3 */
  { 4, 2.0, 2.0 }, /* 4 */
  { 5, 3.0, 2.0 }, /* 5 */
  { 6, 4.0, 2.0 }, /* 6 */
  { 7, 4.0, 4.0 }, /* 7:G */
};

void test(zGraph *graph, data_t *s, data_t *g)
{
  int i;

  /* nodes */
  for( i=0; i<8; i++ )
    zGraphAddNode( graph, &maze[i] );
  /* arcs */
  test_connect( graph, 0, 1, 2 );
  test_connect( graph, 1, 2, 1 );
  test_connect( graph, 1, 3, 1 );
  test_connect( graph, 3, 4, 1 );
  test_connect( graph, 3, 5, 4 );
  test_connect( graph, 4, 5, 1 );
  test_connect( graph, 4, 6, 4 );
  test_connect( graph, 5, 6, 1 );
  test_connect( graph, 6, 7, 2 );
  /* start & goal */
  s->id = 0;
  g->id = 7;
  /* answer: 0-1-3-4-5-6-7 ... 8 */
}
#else
/* TEST2 maze:
 1 S   6o8o9ooooo12
   o   o o o     o
 2 o 3o5o7 o     o
   o o o   o     o
 3 1o2o4   1011oo13
             o
 4 1415oo16oo17oo18
     o   o
 5 oo19oo22ooooGo25
   o o   o       o
 6 2021  23oooooo24
   1 2 3 4 5 6 7 8
 */
static data_t maze[]={
  { 0,  1.0, 1.0 }, /* 0:S */
  { 1,  1.0, 3.0 }, /* 1 */
  { 2,  2.0, 3.0 }, /* 2 */
  { 3,  2.0, 2.0 }, /* 3 */
  { 4,  3.0, 3.0 }, /* 4 */
  { 5,  3.0, 2.0 }, /* 5 */
  { 6,  3.0, 1.0 }, /* 6 */
  { 7,  4.0, 2.0 }, /* 7 */
  { 8,  4.0, 1.0 }, /* 8 */
  { 9,  5.0, 1.0 }, /* 9 */
  { 10, 5.0, 3.0 }, /* 10 */
  { 11, 6.0, 3.0 }, /* 11 */
  { 12, 8.0, 1.0 }, /* 12 */
  { 13, 8.0, 3.0 }, /* 13 */
  { 14, 1.0, 4.0 }, /* 14 */
  { 15, 2.0, 4.0 }, /* 15 */
  { 16, 4.0, 4.0 }, /* 16 */
  { 17, 6.0, 4.0 }, /* 17 */
  { 18, 8.0, 4.0 }, /* 18 */
  { 19, 2.0, 5.0 }, /* 19 */
  { 20, 1.0, 6.0 }, /* 20 */
  { 21, 2.0, 6.0 }, /* 21 */
  { 22, 4.0, 5.0 }, /* 22 */
  { 23, 4.0, 6.0 }, /* 23 */
  { 24, 8.0, 6.0 }, /* 24 */
  { 25, 8.0, 5.0 }, /* 25 */
  { 26, 7.0, 5.0 }, /* 26:G */
};

void test(zGraph *graph, data_t *s, data_t *g)
{
  int i;

  /* nodes */
  for( i=0; i<27; i++ )
    zGraphAddNode( graph, &maze[i] );
  /* arcs */
  test_connect( graph, 0, 1, 2 );
  test_connect( graph, 1, 2, 1 );
  test_connect( graph, 2, 3, 1 );
  test_connect( graph, 2, 4, 1 );
  test_connect( graph, 3, 5, 1 );
  test_connect( graph, 4, 5, 1 );
  test_connect( graph, 5, 6, 1 );
  test_connect( graph, 5, 7, 1 );
  test_connect( graph, 6, 8, 1 );
  test_connect( graph, 7, 8, 1 );
  test_connect( graph, 8, 9, 1 );
  test_connect( graph, 9,10, 2 );
  test_connect( graph, 9,12, 3 );
  test_connect( graph,12,13, 2 );
  test_connect( graph,10,11, 1 );
  test_connect( graph,13,11, 2 );
  test_connect( graph,11,17, 1 );
  test_connect( graph,17,16, 2 );
  test_connect( graph,17,18, 2 );
  test_connect( graph,16,15, 2 );
  test_connect( graph,16,22, 1 );
  test_connect( graph,15,14, 1 );
  test_connect( graph,15,19, 1 );
  test_connect( graph,19,22, 2 );
  test_connect( graph,19,20, 2 );
  test_connect( graph,19,21, 1 );
  test_connect( graph,22,23, 1 );
  test_connect( graph,22,26, 2 );
  test_connect( graph,23,24, 3 );
  test_connect( graph,24,25, 1 );
  test_connect( graph,26,25, 1 );
  /* start & goal */
  s->id = 0;
  g->id = 26;
  /* answer: 0-1-2-3/4-5-6/7-8-9-10-11-17-16-22-26 ... 17 */
}
#endif


void output_path(zGraph *graph, zGraphNodeList *path, double cost)
{
  zGraphNodeListCell *gc;

  zListForEach( path, gc ){
    printf( " -> node" );
    graph->fwrite( stdout, gc->data->data );
    printf( " ... %f\n", gc->data->val );
  }
  printf( "cost = %f\n", cost );
}

int main(int argc, char *argv[])
{
  zGraph graph;
  zGraphNodeList path;
  data_t s, g;
  double cost;

  zGraphInit( &graph );
  graph.dup = test_dup;
  graph.equal = test_equal;
  graph.fwrite = test_fwrite;
  graph.destroy = test_destroy;
  /* test graph */
  test( &graph, &s, &g );
  zGraphFWrite( stdout, &graph );

  graph.h = test_h;
  cost = zGraphAStar( &graph, &s, &g, NULL, &path );
  printf( ">> path (by A*)\n" );
  output_path( &graph, &path, cost );

  cost = zGraphDijkstra( &graph, &s, &g, &path );
  printf( ">> path (by Dijkstra)\n" );
  output_path( &graph, &path, cost );

  zGraphDestroy( &graph );
  zListDestroy( zGraphNodeListCell, &path );
  return 0;
}
