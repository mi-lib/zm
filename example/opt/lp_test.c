#include <zm/zm_opt.h>

bool zLPFRead(FILE *fp, zMat *a, zVec *b, zVec *c, zVec *ans)
{
  
}

int main(int argc, char *argv[])
{
  zMat a;
  zVec b, c, x, ans;
  double cost;
  FILE *fp;

  if( argc < 2 ){
    eprintf( "specify a test problem file.\n" );
    return 1;
  }
  if( !( fp = fopen( argv[1], "r" ) ) ){
    ZOPENERROR( argv[1] );
    return 1;
  }

  a = zMatCloneArray( a_arr, row, col );
  b = zVecCloneArray( b_arr, row );
  c = zVecCloneArray( c_arr, col );
  x = zVecAlloc( col );

  zVecSetSize( &ans, sizeof(answer)/sizeof(double) );
  zVecBuf(&ans) = answer;

  printf( "revised two-phase simplex method\n" );
  printf( "minimize c^T x subject to Ax=b and x>=0, where:\n" );
  printf( "c: " ); zVecWrite( c );
  printf( "A: " ); zMatWrite( a );
  printf( "b: " ); zVecWrite( b );
  printf( "<result>\n" );
  if( !zLPSolveSimplex( a, b, c, x, &cost ) ){
    printf( "failed.\n" );
    return 0;
  }
  zVecWrite( x );
  printf( "cost=%f\n", cost );
  printf( "<true answer>\n" );
  zVecWrite( &ans );
  printf( "cost=%f\n", zRawVecInnerProd( answer, zVecBuf(c), zVecSizeNC(&ans) ) );

  zMatFree( a );
  zVecFree( b );
  zVecFree( c );
  zVecFree( x );
  return 0;
}
