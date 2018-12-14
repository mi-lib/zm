/* seqmod - sequence file modifier */

#include <zm/zm.h>

zSeq seq;
char *name;
char seqfile[BUFSIZ];

/* ******************************************************* */

enum{
  SEQMOD_HELP = 0,
  SEQMOD_COUNT,
  SEQMOD_IP,
  SEQMOD_REV,
  SEQMOD_APPEND,
};
zOption option[] = {
  { "h", "help",        NULL, "display this messages", NULL, false },
  { "c", "count",       NULL, "count steps",           NULL, false },
  { "i", "interpolate", "<dt>", "interpolate sequence every dt second", NULL, false },
  { "r", "reverse",     NULL, "reverse sequence", NULL, false },
  { "a", "append",      NULL, "append sequences", NULL, false },
  { NULL, NULL, NULL, NULL, NULL, false },
};

int seqmodUsage(void)
{
  ZECHO( "usage: seqmod <option> filename [filename]" );
  ZECHO( "<option>" );
  zOptionHelp( option );
  exit( 0 );
}

bool seqmodLoadSequence(char *name)
{
  char buf[BUFSIZ];

  if( !zSeqReadFile( &seq, name ) )
    return false;
  zGetBasename( name, buf, BUFSIZ );
  zCutSuffix( buf );
  strcpy( seqfile, buf );
  return true;
}

void seqmodCount(void)
{
  printf( "0-%d steps\n", zListNum(&seq)-1 );
}

void seqmodInterpolate(double dt)
{
  register int i;
  int dim, step;
  double tt = 0;
  zIP ip;
  zSeq ipseq;
  zSeqListCell *cp;
  zVec v;
  char filename[BUFSIZ];

  dim = zVecSize(zListHead(&seq)->data.v);
  /* total time */
  zListForEach( &seq, cp )
    tt += cp->data.dt;
  /* interpolation */
  if( !( v = zVecAlloc( dim ) ) ){
    ZALLOCERROR();
    exit( 1 );
  }
  zIPCreateSpline( &ip, &seq, ZSPLINE_FIX_EDGE, v, ZSPLINE_FIX_EDGE, v );
  /* zero-fill sequence */
  step = tt / dt;
  zListInit( &ipseq );
  for( i=0; i<=step; i++ ){
    if( !( v = zVecAlloc( dim ) ) || !zSeqEnqueue( &ipseq, v, dt ) ){
      ZALLOCERROR();
      exit( 1 );
    }
    zIPVec( &ip, dt*i, v );
  }
  /* output */
  sprintf( filename, "%s.i", seqfile );
  zAddSuffix( filename, ZSEQ_SUFFIX, filename, BUFSIZ );
  zSeqWriteFile( &ipseq, filename );
  zSeqFree( &ipseq );
  zIPDestroy( &ip );
}

void seqmodReverse(void)
{
  zSeq revseq;
  zSeqListCell *cp;
  char filename[BUFSIZ];

  zSeqInit( &revseq );
  while( !zListIsEmpty( &seq ) ){
    cp = zSeqDequeue( &seq );
    zStackPush( &revseq, cp );
  }
  sprintf( filename, "%s.r", seqfile );
  zAddSuffix( filename, ZSEQ_SUFFIX, filename, BUFSIZ );
  zSeqWriteFile( &revseq, filename );
  zSeqFree( &revseq );
}

bool seqmodAppend(zStrList *arglist)
{
  zStrListCell *cell;
  zSeq subseq;
  char filename[BUFSIZ];

  if( zListNum(arglist) < 2 ) seqmodUsage();
  cell = zListTail(arglist);
  if( !seqmodLoadSequence( cell->data ) )
    return false;
  cell = zListCellNext(cell);
  zListToHead( arglist, cell ){
    if( !zSeqReadFile( &subseq, cell->data ) ) continue;
    zListAppend( &seq, &subseq );
  }
  sprintf( filename, "%s.a.%s", seqfile, ZSEQ_SUFFIX );
  zSeqWriteFile( &seq, filename );
  return true;
}

bool seqmodOperate(int argc, char *argv[])
{
  zStrList arglist;
  bool ret = true;

  if( !zOptionRead( option, argv, &arglist ) ) return false;
  if( option[SEQMOD_HELP].flag )
    seqmodUsage();
  if( option[SEQMOD_APPEND].flag ){
    ret = seqmodAppend( &arglist );
  } else
  if( zListIsEmpty(&arglist) ||
      !seqmodLoadSequence( zListTail(&arglist)->data ) ){
    ret = false;
  }
  zStrListDestroy( &arglist, false );
  if( !ret ) return false;

  if( option[SEQMOD_COUNT].flag )
    seqmodCount();
  if( option[SEQMOD_IP].flag )  /* interpolate sequence */
    seqmodInterpolate( atof(option[SEQMOD_IP].arg) );
  if( option[SEQMOD_REV].flag ) /* reverse sequence */
    seqmodReverse();
  return true;
}

int main(int argc, char *argv[])
{
  if( argc < 2 ) seqmodUsage();
  if( seqmodOperate( argc, argv+1 ) )
    zSeqFree( &seq );
  return 0;
}
