#include <zm/zm_complex.h>

void eval_test(ZTK *ztk)
{
  zComplex c;

  ZTKTagRewind( ztk );
  ZTKKeyRewind( ztk );
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

  if( argc <= 1 ) return 1;

  ZTKInit( &ztk );

  eprintf("\nparsing...%s\n","");
  ZTKParse( &ztk, argv[1] );
  eprintf("done.\n\n&s","");
  eval_test( &ztk );

  ZTKDestroy( &ztk );
  return 0;
}
