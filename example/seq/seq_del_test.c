#include <zm/zm_seq.h>

#define DT 0.1

int main(void)
{
  zSeq seq;
  zVec v;
  zSeqCell *cp;
  double t = 0;
  int i;

  zSeqInit( &seq );
  for( i=0; i<10; i++ ){
    v = zVecCreateList( 1, (double)i );
    zSeqEnqueue( &seq, v, DT );
  }
  while( !zListIsEmpty( &seq ) ){
    cp = zSeqDequeue( &seq );
    printf( "%g ", t );
    zVecPrint( cp->data.v );
    t += cp->data.dt;
    zSeqCellFree( cp );
  }
  return 0;
}
