#include <zm/zm_graph.h>

void *test_dup(void *data)
{
  int *idata;

  if( !( idata = zAlloc( int, 1 ) ) ) return NULL;
  *idata = *(int *)data;
  return idata;
}

bool test_equal(void *d1, void *d2)
{
  return *(int *)d1 == *(int *)d2;
}

void test_fwrite(FILE *fp, void *data)
{
  fprintf( fp, "[%d]", *(int *)data );
}

void test_destroy(void *data)
{
  zFree( data );
}

int main(int argc, char *argv[])
{
  zGraph graph;
  int val, i1, i2;

  zGraphInit( &graph );
  graph.dup = test_dup;
  graph.equal = test_equal;
  graph.fwrite = test_fwrite;
  graph.destroy = test_destroy;
  val = 1; zGraphAddNode( &graph, &val );
  val = 4; zGraphAddNode( &graph, &val );
  val = 3; zGraphAddNode( &graph, &val );
  val = 2; zGraphAddNode( &graph, &val );
  i1 = 1; i2 = 2; zGraphBiconnect( &graph, &i1, &i2, 10 );
  i1 = 1; i2 = 3; zGraphBiconnect( &graph, &i1, &i2, 20 );
  i1 = 2; i2 = 4; zGraphBiconnect( &graph, &i1, &i2, 10 );
  i1 = 3; i2 = 4; zGraphConnect( &graph, &i1, &i2,  5 );
  zGraphFWrite( stdout, &graph );

  printf( ">> destroying\n" );
  zGraphDestroy( &graph );
  zGraphFWrite( stdout, &graph );
  return 0;
}
