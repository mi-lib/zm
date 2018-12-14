#include <zm/zm_parse.h>

#if 0
#define STR "  +12.3e-4 * a - ( 5+b )/7  "
#elif 0
#define STR "2*3+1*4/3"
#elif 0
#define STR "sin(a)"
#elif 0
#define STR "   "
#elif 0
#define STR "(3+3)!"
#elif 0
#define STR "1 2 3"
#elif 1
#define STR "sin(0.5*pi)"
#elif 0
#define STR "2*sin(1+2,3*4+5,6)"
#elif 0
#define STR "testfunc(1+2,3*4+5,6,2*3,5)"
#else
#define STR "{3*(1+2)+4}"
#endif

int main(int argc, char *argv[])
{
  char buf[BUFSIZ];
  zmParser parser;

  argc > 1 ? strcpy( buf, argv[1] ) : strcpy( buf, STR );
  zmParserInit( &parser );
  printf( "expr: %s\n", buf );
  zmParse( &parser, buf );
  printf( "eval = %f\n", zmParseEval( &parser ) );
  zmParserDestroy( &parser );
  return 0;
}
