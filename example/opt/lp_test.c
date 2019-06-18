#include <zm/zm_opt.h>

bool zLPFScan(FILE *fp, zMat *a, zVec *b, zVec *c, zVec *ans)
{
  
}

int main(int argc, char *argv[])
{
  zMat a;
  zVec b, c, x;
  zVecStruct ans;
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

  zVecSetSizeNC( &ans, sizeof(answer)/sizeof(double) );
  zVecBufNC(&ans) = answer;

  printf( "revised two-phase simplex method\n" );
  printf( "minimize c^T x subject to Ax=b and x>=0, where:\n" );
  printf( "c: " ); zVecPrint( c );
  printf( "A: " ); zMatPrint( a );
  printf( "b: " ); zVecPrint( b );
  printf( "<result>\n" );
  if( !zLPSolveSimplex( a, b, c, x, &cost ) ){
    printf( "failed.\n" );
    return 0;
  }
  zVecPrint( x );
  printf( "cost=%f\n", cost );
  printf( "<true answer>\n" );
  zVecPrint( &ans );
  printf( "cost=%f\n", zRawVecInnerProd( answer, zVecBufNC(c), zVecSizeNC(&ans) ) );

  zMatFree( a );
  zVecFree( b );
  zVecFree( c );
  zVecFree( x );
  return 0;
}
