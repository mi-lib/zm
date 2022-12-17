#include <zm/zm_seq.h>

int main(void)
{
  zSeq seq;
  zSeqListCell *cp;
  register int i=0;

  zSeqScanFile( &seq, "test" );
  zListForEach( &seq, cp ){
    if( !cp ) continue;
    printf( "jump to->%2d : %f ", i, cp->data.dt );
    zVecPrint( cp->data.v );
    i++;
  }
  zSeqFree( &seq );
  return 0;
}
