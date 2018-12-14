#include <zm/zm_ip.h>

#define N 5
#define T 5.0

int main(void)
{
  zFerguson ferg[N];
  double v1, v2, t;
  int i;

  for( v1=-1.0, v2=1.0, i=0; i<N; i++, v1+=0.5, v2-=0.5 )
    zFergusonCreate( &ferg[i], T, 0, v1, 1.5, v2 );

  for( t=0; t<T; t+=0.01 ){
    for( i=0; i<N; i++ )
      printf( "%f ", zFergusonValue( &ferg[i], t ) );
    printf( "\n" );
  }

  return 0;
}
