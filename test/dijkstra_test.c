#include <zm/zm_graph.h>

void *test_dup(void *data)
{
  int *idata;

  if( !( idata = zAlloc( int, 1 ) ) ) return NULL;
  *idata = *(int *)data;
  return idata;
}

bool test_equal(void *d1, void *d2){ return *(int *)d1 == *(int *)d2; }
void test_fprint(FILE *fp, void *data){ fprintf( fp, "%d", *(int *)data ); }
void test_destroy(void *data){ zFree( data ); }
void test_add_node(zGraph *graph, int val){ zGraphAddNode( graph, &val ); }
void test_connect(zGraph *graph, int i1, int i2, double cost){ zGraphBiconnect( graph, &i1, &i2, cost ); }

double test1(zGraph *graph, int *s, int *g, int **ap)
{
  static int answer[] = { 1, 3, 5, 6, -1 };

  test_add_node( graph, 1 );
  test_add_node( graph, 2 );
  test_add_node( graph, 3 );
  test_add_node( graph, 4 );
  test_add_node( graph, 5 );
  test_add_node( graph, 6 );
  test_connect( graph, 1, 2, 5 );
  test_connect( graph, 1, 3, 4 );
  test_connect( graph, 1, 4, 2 );
  test_connect( graph, 2, 3, 2 );
  test_connect( graph, 3, 4, 3 );
  test_connect( graph, 3, 5, 2 );
  test_connect( graph, 4, 5, 6 );
  test_connect( graph, 2, 6, 6 );
  test_connect( graph, 5, 6, 4 );

  *s = 1;
  *g = 6;
  *ap = answer;
  return 10.0;
}

double test2(zGraph *graph, int *s, int *g, int **ap)
{
  static int answer[] = { 4, 3, 1, 2, 5, -1 };

  test_add_node( graph, 1 );
  test_add_node( graph, 2 );
  test_add_node( graph, 3 );
  test_add_node( graph, 4 );
  test_add_node( graph, 5 );
  test_connect( graph, 1, 2, 2 );
  test_connect( graph, 1, 3, 1 );
  test_connect( graph, 1, 4, 7 );
  test_connect( graph, 2, 3, 5 );
  test_connect( graph, 2, 4, 6 );
  test_connect( graph, 2, 5, 1 );
  test_connect( graph, 3, 4, 1 );
  test_connect( graph, 4, 5, 8 );

  *s = 4;
  *g = 5;
  *ap = answer;
  return 5.0;
}

double test3(zGraph *graph, int *s, int *g, int **ap)
{
  static int answer[] = { 1, 2, 4, 6, -1 };

  test_add_node( graph, 1 );
  test_add_node( graph, 2 );
  test_add_node( graph, 3 );
  test_add_node( graph, 4 );
  test_add_node( graph, 5 );
  test_add_node( graph, 6 );
  test_connect( graph, 1, 2, 10 );
  test_connect( graph, 1, 3, 10 );
  test_connect( graph, 2, 3, 10 );
  test_connect( graph, 2, 4, 8 );
  test_connect( graph, 3, 5, 10 );
  test_connect( graph, 4, 6, 8 );
  test_connect( graph, 5, 6, 10 );

  *s = 1;
  *g = 6;
  *ap = answer;
  return 26.0;
}

double test4(zGraph *graph, int *s, int *g, int **ap)
{
  static int answer[] = { 1, 4, 7, 8, -1 };

  test_add_node( graph, 1 );
  test_add_node( graph, 2 );
  test_add_node( graph, 3 );
  test_add_node( graph, 4 );
  test_add_node( graph, 5 );
  test_add_node( graph, 6 );
  test_add_node( graph, 7 );
  test_add_node( graph, 8 );
  test_connect( graph, 1, 2, 1 );
  test_connect( graph, 1, 3, 7 );
  test_connect( graph, 1, 4, 2 );
  test_connect( graph, 2, 5, 2 );
  test_connect( graph, 2, 6, 4 );
  test_connect( graph, 3, 6, 2 );
  test_connect( graph, 3, 7, 3 );
  test_connect( graph, 4, 7, 5 );
  test_connect( graph, 5, 6, 1 );
  test_connect( graph, 6, 8, 6 );
  test_connect( graph, 7, 8, 2 );

  *s = 1;
  *g = 8;
  *ap = answer;
  return 9.0;
}

bool assert_dijkstra_case(double (* testfunc)(zGraph*,int*,int*,int**))
{
  zGraph graph;
  zGraphNodeList path;
  zGraphNodeListCell *gc;
  double cost, cost_a=HUGE_VAL;
  int s, g, i, *answer=NULL;

  zGraphInit( &graph );
  graph.dup = test_dup;
  graph.equal = test_equal;
  graph.fprint = test_fprint;
  graph.destroy = test_destroy;

  cost_a = testfunc( &graph, &s, &g, &answer );
  cost = zGraphSearchDijkstra( &graph, &s, &g, &path );

  printf( "\n  path: " );
  zListForEach( &path, gc ){
    graph.fprint( stdout, gc->data->data );
    printf( " - " );
  }
  printf( "cost = %f\n", cost );
  printf( "answer: " );
  for( i=0; answer[i]!=-1; i++ )
    printf( "%d - ", answer[i] );
  printf( "cost = %f\n", cost_a );

  zGraphDestroy( &graph );
  zListDestroy( zGraphNodeListCell, &path );
  return zIsEqual( cost, cost_a, zTOL );
}

void assert_dijkstra(void)
{
  zAssert( zGraphSearchDijkstra (test case1), assert_dijkstra_case( test1 ) );
  zAssert( zGraphSearchDijkstra (test case2), assert_dijkstra_case( test2 ) );
  zAssert( zGraphSearchDijkstra (test case3), assert_dijkstra_case( test3 ) );
  zAssert( zGraphSearchDijkstra (test case4), assert_dijkstra_case( test4 ) );
}

int main(int argc, char *argv[])
{
  assert_dijkstra();
  return EXIT_SUCCESS;
}
