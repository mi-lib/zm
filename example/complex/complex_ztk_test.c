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

void reg_tf_def(ZTK *ztk)
{
  char *key[] = {
    "complex",
  };
  ZTKDefReg( ztk, "val", key );
}

int main(int argc, char *argv[])
{
  ZTK ztk;

  if( argc <= 1 ) return 1;

  ZTKInit( &ztk );
  reg_tf_def( &ztk );

  eprintf("\nparsing...\n");
  ZTKParse( &ztk, argv[1] );
  eprintf("done.\n\n");
  eval_test( &ztk );

  ZTKDestroy( &ztk );
  return 0;
}
