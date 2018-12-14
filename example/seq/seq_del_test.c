#include <zm/zm_seq.h>

#define DT 0.1

int main(void)
{
  zSeq seq;
  zVec v;
  zSeqListCell *cp;
  double t = 0;
  register int i;

  zSeqInit( &seq );
  for( i=0; i<10; i++ ){
    v = zVecCreateList( 1, (double)i );
    zSeqEnqueue( &seq, v, DT );
  }
  while( !zListIsEmpty( &seq ) ){
    cp = zSeqDequeue( &seq );
    printf( "%f ", t );
    zVecWrite( cp->data.v );
    t += cp->data.dt;
    zSeqListCellFree( cp );
  }
  return 0;
}
