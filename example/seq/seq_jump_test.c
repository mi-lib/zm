#include <zm/zm_seq.h>

int main(void)
{
  zSeq seq;
  zSeqListCell *cp;
  int i;

  zSeqScanFile( &seq, "test" );
  /* a warning is issued at the end. */
  for( i=0; i<=zListSize(&seq); i++ ){
    if( !( cp = zSeqJump( &seq, i ) ) ) continue;
    printf( "jump to %2d : %g ", i, cp->data.dt );
    zVecPrint( cp->data.v );
  }
  zSeqFree( &seq );
  return 0;
}
