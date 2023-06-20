#include <zm/zm_complex.h>

void eval_test(ZTK *ztk)
{
  zComplex c;

  ZTKRewind( ztk );
  do{
    if( ZTKKeyCmp( ztk, "complex" ) ){
      eprintf( "original string = %s\t-> ", ZTKVal(ztk) );
      zComplexPrint( zComplexFromZTK( &c, ztk ) );
      zEndl();
    }
  } while( ZTKKeyNext(ztk) );
}

int main(int argc, char *argv[])
{
  ZTK ztk;

  ZTKInit( &ztk );

  eprintf("\nparsing...\n");
  ZTKParse( &ztk, "complex.ztk" );
  eprintf("done.\n\n");
  ZTKPrint( &ztk );
  eval_test( &ztk );

  ZTKDestroy( &ztk );
  return 0;
}
