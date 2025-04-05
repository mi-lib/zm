#include <zm/zm.h>

enum{
  HISTOGRAM_HELP = 0,
  HISTOGRAM_NUM,
  HISTOGRAM_OUTPUT,
};
zOption option[] = {
  { "h",    "help",     NULL, "display this messages",   NULL, false },
  { "n",    "num", "<value>", "number of bins",          "10", false },
  { "o", "output", "<name>",  "basename of output files", "h", false },
  { NULL, NULL, NULL, NULL, NULL, false },
};

int zm_histogram_usage(char *argv)
{
  ZECHO( "Usage: %s <filename>\n", argv );
  ZECHO( "<option>" );
  zOptionHelp( option );
  exit( EXIT_SUCCESS );
}

FILE *zm_histogram_read_option(char *argv[])
{
  zStrAddrList arglist;
  FILE *fp;

  if( !zOptionRead( option, argv+1, &arglist ) )
    exit( EXIT_FAILURE );
  if( option[HISTOGRAM_HELP].flag )
    zm_histogram_usage( argv[0] );
  if( zListIsEmpty( &arglist ) )
    return stdin;
  if( !( fp = fopen( zListHead(&arglist)->data, "r" ) ) ){
    ZOPENERROR( zListHead(&arglist)->data );
    exit( EXIT_FAILURE );
  }
  return fp;
}

bool zm_histogram_sample(FILE *fp, zVecList *list)
{
  char linebuf[BUFSIZ], tokenbuf[BUFSIZ];
  zVec v;
  int field_num = 0, i;
  double val;

  zListInit( list );
  while( !feof( fp ) ){
    if( !fgets( linebuf, BUFSIZ, fp ) ) continue;
    strcpy( tokenbuf, linebuf );
    for( field_num=0; zSDouble( tokenbuf, &val ); field_num++ );
    if( field_num == 0 ) continue;
    if( !( v = zVecAlloc( field_num ) ) ){
      ZALLOCERROR();
      return false;
    }
    strcpy( tokenbuf, linebuf );
    for( i=0; i<field_num; i++ ){
      zSDouble( tokenbuf, &val );
      zVecSetElemNC( v, i, val );
    }
    if( !zVecListInsertHead( list, v ) ){
      ZALLOCERROR();
      return false;
    }
  }
  if( zListIsEmpty( list ) ){
    ZRUNERROR( "empty data file" );
    return false;
  }
  return true;
}

bool zm_histogram_output(zVecList *list)
{
  zHistogram histogram;
  zVec data;
  zVecListCell *cp;
  int i, data_id;
  char filename[BUFSIZ];
  FILE *fp;
  bool retval = true;

  if( !( data = zVecAlloc( zListSize(list) ) ) ){
    ZALLOCERROR();
    return false;
  }
  for( data_id=0; data_id<zVecSizeNC(zListTail(list)->data); data_id++ ){
    i = 0;
    zListForEach( list, cp ){
      zVecSetElem( data, i, zVecElem(cp->data,data_id) );
      i++;
    }
    if( !zHistogramCreateAuto( &histogram, zVecBuf(data), zVecSizeNC(data), atoi(option[HISTOGRAM_NUM].arg) ) ){
      retval = false;
      break;
    }
    sprintf( filename, "%s%d", option[HISTOGRAM_OUTPUT].arg, data_id );
    if( !( fp = fopen( filename, "w" ) ) ){
      ZOPENERROR( filename );
      retval = false;
      break;
    }
    zHistogramFPrint( fp, &histogram );
    fclose( fp );
  }
  zVecFree( data );
  return retval;
}

int main(int argc, char *argv[])
{
  FILE *fp;
  zVecList list;

  fp = zm_histogram_read_option( argv );
  if( zm_histogram_sample( fp, &list ) )
    zm_histogram_output( &list );
  zVecListDestroy( &list );
  if( fp != stdin ) fclose( fp );
  return 0;
}
