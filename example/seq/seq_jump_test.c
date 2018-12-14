#include <zm/zm_seq.h>

int main(void)
{
  zSeq seq;
  zSeqListCell *cp;
  register int i;

  zSeqReadFile( &seq, "test" );
  for( i=0; i<zListNum(&seq)+1; i++ ){
    if( !( cp = zSeqJump( &seq, i ) ) ) continue;
    printf( "jump to->%2d : %f ", i, cp->data.dt );
    zVecWrite( cp->data.v );
  }
  zSeqFree( &seq );
  return 0;
}
