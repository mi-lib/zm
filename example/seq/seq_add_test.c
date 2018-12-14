#include <zm/zm_seq.h>

#define DT 0.1

int main(void)
{
  zSeq seq;
  zVec v;
  register int i;

  zSeqInit( &seq );
  for( i=0; i<10; i++ ){
    v = zVecCreateList( 1, (double)i );
    zSeqEnqueue( &seq, v, DT );
  }
  zSeqWriteFile( &seq, "dummy" );
  zSeqFree( &seq );
  return 0;
}
