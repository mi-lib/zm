#include <zm/zm_graph.h>

typedef struct{
  double x, y;
  bool check;
} data_t;

void *test_dup(void *data)
{
  data_t *d;

  if( !( d = zAlloc( data_t, 1 ) ) ) return NULL;
  d->x  = ((data_t *)data)->x;
  d->y  = ((data_t *)data)->y;
  d->check = false;
  return d;
}

bool test_equal(void *d1, void *d2)
{
  return ((data_t *)d1)->x == ((data_t *)d2)->x
      && ((data_t *)d1)->y == ((data_t *)d2)->y;
}

void test_fwrite(FILE *fp, void *data)
{
  fprintf( fp, "(%g,%g)", ((data_t *)data)->x, ((data_t *)data)->y );
}

void test_destroy(void *data)
{
  zFree( data );
}

void test_add(zGraph *graph, data_t *d1, double dx, double dy)
{
  data_t d;

  d.x = d1->x + dx;
  d.y = d1->y + dy;
  zGraphAddNode( graph, &d );
  zGraphBiconnect( graph, d1, &d, fabs(dx)+fabs(dy) );
}

double test_h(void *n1, void *n2, void *util)
{
  data_t *d1, *d2;
  double l;

  d1 = n1;
  d2 = n2;
  l = fabs( d1->x - d2->x ) + fabs( d1->y - d2->y );
  if( !d1->check ){
    d1->check = true;
    if( l < 1.0 ){
      zGraphBiconnect( util, n1, d2, l );
    } else{
      test_add( util, d1, 1.0, 0.0 );
      test_add( util, d1, 1.0, 1.0 );
      test_add( util, d1, 0.0, 1.0 );
      test_add( util, d1,-1.0, 1.0 );
      test_add( util, d1,-1.0, 0.0 );
      test_add( util, d1,-1.0,-1.0 );
      test_add( util, d1, 0.0,-1.0 );
      test_add( util, d1, 1.0,-1.0 );
    }
  }
  return l;
}

void node_create(zGraph *graph, data_t *node, double x, double y)
{
  node->x = x;
  node->y = y;
  node->check = false;
  zGraphAddNode( graph, node );
}

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
  graph.h = test_h;

  /* test graph */
  node_create( &graph, &s, 0, 0 );
  node_create( &graph, &g, 5, 10 );
  printf( ">> initial graph\n" );
  zGraphFWrite( stdout, &graph );

  printf( ">> start computing the shortest path...\n" );
  cost = zGraphAStar( &graph, &s, &g, &graph, &path );
  printf( "   done.\n" );
  printf( ">> final graph\n" );
  zGraphFWrite( stdout, &graph );
  printf( ">> path\n" );
  output_path( &graph, &path, cost );

  zGraphDestroy( &graph );
  zListDestroy( zGraphNodeListCell, &path );
  return 0;
}
