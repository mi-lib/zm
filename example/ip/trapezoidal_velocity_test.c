#include <zm/zm_ip.h>

#define DT 0.01

int main(int argc, char *argv[])
{
  zTRVelProf trvelprof;
  double t;
  int i, step;

#ifdef __cplusplus
  trvelprof.create( -1, -2, 5, 4, argc > 1 ? atof(argv[1]) : 5 );
  step = trvelprof.term() / DT + 1;
  for( i=0; i<=step; i++ ){
    t = i * DT;
    printf( "%g %g %g\n", trvelprof.dist(t), trvelprof.vel(t), trvelprof.acc(t) );
  }
#else
  zTRVelProfCreate( &trvelprof, -1, -2, 5, 4, argc > 1 ? atof(argv[1]) : 5 );
  step = zTRVelProfTerm( &trvelprof ) / DT + 1;
  for( i=0; i<=step; i++ ){
    t = i * DT;
    printf( "%g %g %g\n", zTRVelProfDist(&trvelprof,t), zTRVelProfVel(&trvelprof,t), zTRVelProfAcc(&trvelprof,t) );
  }
#endif
  return 0;
}
