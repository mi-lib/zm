/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_parse - parser of mathematical expression.
 */

#include <zm/zm_parse.h>
#include <zm/zm_stat.h>

static zmParseNum *_zmParseNumAllocName(char *str);
static zmParseNum *_zmParseNumAllocVal(char *str, double val);
static zmParseNum *_zmParseNumAlloc(char *str);
static void _zmParseNumFree(zmParseNum *num);
static double _zmParseNumEval(zmParseNum *num);

static zmParseVar *_zmParseVarListReg(zmParseVarList *list, char *str);
static void _zmParseVarListDestroy(zmParseVarList *list);
static bool _zmParseVarAssoc(zmParseVarList *list, char *str, zmParseCell *cell);
static bool _zmParseVarBuiltin(zmParseVarList *list);
static double _zmParseVarEval(zmParseVar *var);
#ifdef DEBUG
static void _zmParseVarListWrite(zmParseVarList *list);
#endif

static zmParseFunc *_zmParseFuncListReg(zmParseFuncList *list, char *str, zmParseFuncID id, int argc);
static void _zmParseFuncListDestroy(zmParseFuncList *list);
static bool _zmParseFuncBuiltin(zmParseFuncList *list);
static int _zmParseFuncGetArg(zmParseCell *arg, int argc, double *argv);
static double _zmParseFuncEval(zmParseFunc *func, zmParseCell *arg);
#ifdef DEBUG
static void _zmParseFuncListWrite(zmParseFuncList *list);
#endif

static zmParseOp *_zmParseOpFind(char *str);
static double _zmParseOpEval(zmParseOp *op, double arg1, double arg2);

static zmParseCell *_zmParseCellAlloc(zmParseCellType type, char *str, zmParseVarList *varlist, zmParseFuncList *funclist);
static void _zmParseCellFree(zmParseCell *cell);
static double _zmParseCellEval(zmParseCell *cell);
static char *_zmParseCellStr(zmParseCell *cell);
#ifdef DEBUG
static void _zmParseCellWrite(zmParseCell *cell, int indent);
#endif

static char _zmParseParChar(zmParseStatus status);
static void _zmParseTokFree(zmParseTok *tok);
static int _zmParseTokLevel(zmParseTok *tok);
static zmParseStatus _zmParseTokStatus(zmParseTok *tok);
static int _zmParseTokOpSpec(zmParseTok *tok);
static int _zmParseTokPriority(zmParseTok *tok);
static char *_zmParseTokStr(zmParseTok *tok);
#ifdef DEBUG
static void _zmParseTokWrite(zmParseTok *tok);
#endif

static void _zmParseStackInit(zmParseTok *stack);
static void _zmParseStackDestroy(zmParseTok *stack);
#ifdef DEBUG
static void _zmParseStackWrite(zmParseTok *stack);
#endif

static zmParseTok *_zmParseTokAlloc(zmParser *parser, zmParseCellType type, char *str);
static void _zmParseStackPush(zmParser *parser, zmParseTok *tok, int level, zmParseStatus status);
static zmParseTok *_zmParseStackPop(zmParser *parser, zmParseTok *toknext);

static zmParseTok *_zmParseTokNum(zmParser *parser, char *str, char *tkn, size_t size);
static zmParseTok *_zmParseTokVar(zmParser *parser, char *str, char *tkn, size_t size);
static zmParseTok *_zmParseTokOp(zmParser *parser, char *str, char *tkn, size_t size);
static zmParseTok *_zmParseTokFunc(zmParser *parser, char *str, char *tkn, size_t size);
static zmParseTok *_zmParseTokPar(zmParser *parser, char *str, char *tkn, size_t size);

static bool _zmParseStackPrefix(zmParser *parser, zmParseTok *toknext);
static bool _zmParseStackInfix(zmParser *parser, zmParseTok *toknext);
static bool _zmParseStackPreinfix(zmParser *parser, zmParseTok *toknext);
static bool _zmParseStackPostfix(zmParser *parser, zmParseTok *toknext);
static bool _zmParseStack(zmParser *parser, zmParseTok *toknext);

static bool _zmParseIn(zmParser *parser, zmParseStatus status);
static zmParseTok *_zmParsePar(zmParser *parser, zmParseTok *tok);

/* set of delimiters and operators */

static char __zm_parse_delim[] = {
  EOF, '\t', '\v', '\f', '\n', '\r', ' ', '\0'
};

static char __zm_parse_op[] = "+-*/%^!,;:=";

/* number */

zmParseNum *_zmParseNumAllocName(char *str)
{
  zmParseNum *num;

  if( !( num = zAlloc( zmParseNum, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !zNameSet( num, str ) ){
    ZALLOCERROR();
    zFree( num );
    return NULL;
  }
  return num;
}

zmParseNum *_zmParseNumAllocVal(char *str, double val)
{
  zmParseNum *num;

  if( ( num = _zmParseNumAllocName( str ) ) )
    num->val = val;
  return num;
}

zmParseNum *_zmParseNumAlloc(char *str)
{
  return _zmParseNumAllocVal( str, atof( str ) );
}

void _zmParseNumFree(zmParseNum *num)
{
  zNameDestroy( num );
  zFree( num );
}

double _zmParseNumEval(zmParseNum *num){ return num->val; }

/* list of variables */

zmParseVar *_zmParseVarListReg(zmParseVarList *list, char *str)
{
  zmParseVarCell *cp;

  zListForEach( list, cp ){
    if( strcmp( zName(&cp->data), str ) == 0 )
      return &cp->data;
  }
  if( !( cp = zAlloc( zmParseVarCell, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !zNameSet( &cp->data, str ) ){
    ZALLOCERROR();
    zFree( cp );
    return NULL;
  }
  cp->data.cell = NULL; /* dummy */
  zListInsertHead( list, cp );
  return &cp->data;
}

void _zmParseVarListDestroy(zmParseVarList *list)
{
  zmParseVarCell *cp;

  while( !zListIsEmpty( list ) ){
    zListDeleteHead( list, &cp );
    zNameDestroy( &cp->data );
    zFree( cp );
  }
}

bool _zmParseVarAssoc(zmParseVarList *list, char *str, zmParseCell *cell)
{
  zmParseVar *var;

  if( !( var = _zmParseVarListReg( list, str ) ) )
    return false;
  if( var->cell )
    _zmParseCellFree( var->cell );
  var->cell = cell;
  return true;
}

bool _zmParseVarBuiltin(zmParseVarList *list)
{
  zmParseCell *c_pi, *c_e;

  if( !( c_pi = zAlloc( zmParseCell, 1 ) ) ||
      !( c_pi->dat.num = _zmParseNumAllocVal( "pi", zPI ) ) ){
    ZALLOCERROR();
    return false;
  }
  if( !( c_e = zAlloc( zmParseCell, 1 ) ) ||
      !( c_e->dat.num = _zmParseNumAllocVal( "e", exp(1.0) ) ) ){
    ZALLOCERROR();
    return false;
  }
  c_pi->type = c_e->type = ZM_PARSE_CELL_NUM;
  c_pi->arg1 = c_pi->arg2 = NULL;
  c_e->arg1 = c_e->arg2 = NULL;
  if( !_zmParseVarAssoc( list, "pi", c_pi ) ||
      !_zmParseVarAssoc( list, "e",  c_e  ) ){
    return false;
  }
  return true;
}

double _zmParseVarEval(zmParseVar *var)
{
  return var->cell ? _zmParseCellEval( var->cell ) : 0;
}

#ifdef DEBUG
void _zmParseVarListWrite(zmParseVarList *list) /* for debug */
{
  zmParseVarCell *cp;

  zListForEach( list, cp ){
    printf( "%s = %f\n", zName(&cp->data), _zmParseVarEval(&cp->data) );
  }
}
#endif

/* list of functions */

zmParseFunc *_zmParseFuncListReg(zmParseFuncList *list, char *str, zmParseFuncID id, int argc)
{
  zmParseFuncCell *cp;

  zListForEach( list, cp ){
    if( strcmp( zName(&cp->data), str ) == 0 )
      return &cp->data;
  }
  if( !( cp = zAlloc( zmParseFuncCell, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !zNameSet( &cp->data, str ) ){
    ZALLOCERROR();
    zFree( cp );
    return NULL;
  }
  cp->data.id = id;
  cp->data.argc = argc;
  cp->data.cell = NULL;
  zListInsertHead( list, cp );
  return &cp->data;
}

void _zmParseFuncListDestroy(zmParseFuncList *list)
{
  zmParseFuncCell *cp;

  while( !zListIsEmpty( list ) ){
    zListDeleteHead( list, &cp );
    zNameDestroy( &cp->data );
    zFree( cp );
  }
}

bool _zmParseFuncBuiltin(zmParseFuncList *list)
{
  if( !_zmParseFuncListReg( list, "sin", ZM_PARSE_FUNC_SIN, 1 ) ||
      !_zmParseFuncListReg( list, "cos", ZM_PARSE_FUNC_COS, 1 ) ||
      !_zmParseFuncListReg( list, "tan", ZM_PARSE_FUNC_TAN, 1 ) ||
      !_zmParseFuncListReg( list, "cot", ZM_PARSE_FUNC_COT, 1 ) ||
      !_zmParseFuncListReg( list, "asin", ZM_PARSE_FUNC_ASIN, 1 ) ||
      !_zmParseFuncListReg( list, "acos", ZM_PARSE_FUNC_ACOS, 1 ) ||
      !_zmParseFuncListReg( list, "atan", ZM_PARSE_FUNC_ATAN, 1 ) ||
      !_zmParseFuncListReg( list, "acot", ZM_PARSE_FUNC_ACOT, 1 ) ||
      !_zmParseFuncListReg( list, "sinh", ZM_PARSE_FUNC_SINH, 1 ) ||
      !_zmParseFuncListReg( list, "cosh", ZM_PARSE_FUNC_COSH, 1 ) ||
      !_zmParseFuncListReg( list, "tanh", ZM_PARSE_FUNC_TANH, 1 ) ||
      !_zmParseFuncListReg( list, "coth", ZM_PARSE_FUNC_COTH, 1 ) ||
      !_zmParseFuncListReg( list, "exp", ZM_PARSE_FUNC_EXP, 1 ) ||
      !_zmParseFuncListReg( list, "log", ZM_PARSE_FUNC_LOG, 1 ) ||
      !_zmParseFuncListReg( list, "ln", ZM_PARSE_FUNC_LN, 1 ) ||
      !_zmParseFuncListReg( list, "abs", ZM_PARSE_FUNC_ABS, 1 ) ||
      !_zmParseFuncListReg( list, "sqrt", ZM_PARSE_FUNC_SQRT, 1 ) ||
      !_zmParseFuncListReg( list, "deg", ZM_PARSE_FUNC_DEG, 1 ) ||
      !_zmParseFuncListReg( list, "rad", ZM_PARSE_FUNC_RAD, 1 ) ||
      !_zmParseFuncListReg( list, "testfunc", ZM_PARSE_FUNC_USER, 3 ) ){
    return false;
  }
  return true;
}

int _zmParseFuncGetArg(zmParseCell *arg, int argc, double *argv)
{
  int c = 0;

  if( !arg ) return 0;
  if( arg->type == ZM_PARSE_CELL_CHUNK &&
      arg->dat.op->type == ZM_PARSE_OP_COMMA ){
    c = _zmParseFuncGetArg( arg->arg1, argc, argv );
    if( c >= argc ) return c;
    argv[c] = _zmParseCellEval( arg->arg2 );
    return c+1;
  }
  *argv = _zmParseCellEval( arg );
  return 1;
}

double _zmParseFuncEval(zmParseFunc *func, zmParseCell *arg)
{
  double *argv, ret = 0;

  if( !( argv = zAlloc( double, func->argc ) ) ){
    ZALLOCERROR();
    return 0;
  }
  _zmParseFuncGetArg( arg, func->argc, argv );

  switch( func->id ){
  case ZM_PARSE_FUNC_SIN:  ret = sin( argv[0] ); break;
  case ZM_PARSE_FUNC_COS:  ret = cos( argv[0] ); break;
  case ZM_PARSE_FUNC_TAN:  ret = tan( argv[0] ); break;
  case ZM_PARSE_FUNC_COT:  ret = 1.0 / tan( argv[0] ); break;
  case ZM_PARSE_FUNC_ASIN: ret = asin( argv[0] ); break;
  case ZM_PARSE_FUNC_ACOS: ret = acos( argv[0] ); break;
  case ZM_PARSE_FUNC_ATAN: ret = atan( argv[0] ); break;
  case ZM_PARSE_FUNC_ACOT: ret = atan( 1.0 / argv[0] ); break;
  case ZM_PARSE_FUNC_SINH: ret = sinh( argv[0] ); break;
  case ZM_PARSE_FUNC_COSH: ret = cosh( argv[0] ); break;
  case ZM_PARSE_FUNC_TANH: ret = tanh( argv[0] ); break;
  case ZM_PARSE_FUNC_COTH: ret = 1.0 / tanh( argv[0] ); break;
  case ZM_PARSE_FUNC_EXP:  ret = exp( argv[0] ); break;
  case ZM_PARSE_FUNC_LOG:  ret = log10( argv[0] ); break;
  case ZM_PARSE_FUNC_LN:   ret = log( argv[0] ); break;
  case ZM_PARSE_FUNC_ABS:  ret = fabs( argv[0] ); break;
  case ZM_PARSE_FUNC_SQRT: ret = sqrt( argv[0] ); break;
  case ZM_PARSE_FUNC_DEG:  ret = zRad2Deg( argv[0] ); break;
  case ZM_PARSE_FUNC_RAD:  ret = zDeg2Rad( argv[0] ); break;
  case ZM_PARSE_FUNC_USER: ret = 0; break; /* dummy */
  default: ;
  }
  zFree( argv );
  return ret;
}

#ifdef DEBUG
void _zmParseFuncListWrite(zmParseFuncList *list) /* for debug */
{
  zmParseFuncCell *cp;

  zListForEach( list, cp ){
    printf( "[#%d] %s <%d>\n", cp->data.id, zName(&cp->data), cp->data.argc );
    if( cp->data.cell )
      _zmParseCellWrite( cp->data.cell, 0 );
  }
}
#endif

/* operator */

zmParseOp *_zmParseOpFind(char *str)
{
  zmParseOp *op;

  for( op=(zmParseOp *)zm_parse_op; op->str; op++ ){
    if( strcmp( op->str, str ) == 0 )
      return op;
  }
  return NULL;
}

double _zmParseOpEval(zmParseOp *op, double arg1, double arg2)
{
  switch( op->type ){
  case ZM_PARSE_OP_ADD: return arg1 + arg2;
  case ZM_PARSE_OP_SUB: return arg1 - arg2;
  case ZM_PARSE_OP_MUL: return arg1 * arg2;
  case ZM_PARSE_OP_DIV:
  case ZM_PARSE_OP_RAT:
    return zIsTiny( arg2 ) ? ( arg2 >= 0 ? HUGE_VAL : -HUGE_VAL ) : arg1 / arg2;
  case ZM_PARSE_OP_MOD: return (int)arg1 % (int)arg2;
  case ZM_PARSE_OP_POW: return pow( arg1, arg2 );
  case ZM_PARSE_OP_FCT: return zFacto( arg1 );
  case ZM_PARSE_OP_COMMA: return arg1; /* dummy */
  case ZM_PARSE_OP_SEP:   return arg1; /* dummy */
  case ZM_PARSE_OP_LET:   return arg2; /* dummy */
  default: ;
  }
  return 0;
}

/* cell */

zmParseCell *_zmParseCellAlloc(zmParseCellType type, char *str, zmParseVarList *varlist, zmParseFuncList *funclist)
{
  zmParseCell *cell;

  if( !( cell = zAlloc( zmParseCell, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  switch( ( cell->type = type ) ){
  case ZM_PARSE_CELL_NUM:
    if( !( cell->dat.num = _zmParseNumAlloc( str ) ) )
      return NULL;
    break;
  case ZM_PARSE_CELL_VAR:
    cell->dat.var = _zmParseVarListReg( varlist, str );
    break;
  case ZM_PARSE_CELL_OP:
    cell->dat.op = _zmParseOpFind( str );
    break;
  case ZM_PARSE_CELL_FUNC:
    cell->dat.func = _zmParseFuncListReg( funclist, str, ZM_PARSE_FUNC_USER, 0 );
    break;
  case ZM_PARSE_CELL_PAR:
    strcpy( cell->dat.parident, str );
    break;
  default: ;
  }
  cell->arg1 = cell->arg2 = NULL;
  return cell;
}

void _zmParseCellFree(zmParseCell *cell)
{
  if( cell->type == ZM_PARSE_CELL_NUM )
    _zmParseNumFree( cell->dat.num );
  if( cell->arg1 ){
    _zmParseCellFree( cell->arg1 );
    zFree( cell->arg1 );
  }
  if( cell->arg2 ){
    _zmParseCellFree( cell->arg2 );
    zFree( cell->arg2 );
  }
}

double _zmParseCellEval(zmParseCell *cell)
{
  double arg1, arg2;

  switch( cell->type ){
  case ZM_PARSE_CELL_NUM: return _zmParseNumEval( cell->dat.num );
  case ZM_PARSE_CELL_VAR: return _zmParseVarEval( cell->dat.var );
  case ZM_PARSE_CELL_FUNC: return _zmParseFuncEval( cell->dat.func, cell->arg1 );
  case ZM_PARSE_CELL_CHUNK:
    arg1 = cell->arg1 ? _zmParseCellEval( cell->arg1 ) : 0;
    arg2 = cell->arg2 ? _zmParseCellEval( cell->arg2 ) : 0;
    return _zmParseOpEval( cell->dat.op, arg1, arg2 );
  default: ;
  }
  return 0;
}

char *_zmParseCellStr(zmParseCell *cell)
{
  if( !cell ) return zNullStr();
  switch( cell->type ){
  case ZM_PARSE_CELL_NUM:   return zName( cell->dat.num );
  case ZM_PARSE_CELL_VAR:   return zName( cell->dat.var );
  case ZM_PARSE_CELL_FUNC:  return zName( cell->dat.func );
  case ZM_PARSE_CELL_OP:    return (char *)cell->dat.op->str;
  case ZM_PARSE_CELL_PAR:   return (char *)cell->dat.parident;
  case ZM_PARSE_CELL_CHUNK: return _zmParseCellStr( cell->arg1 );
  default: ;
  }
  return NULL;
}

#ifdef DEBUG
void _zmParseCellWrite(zmParseCell *cell, int indent) /* for debug */
{
  zIndent( indent );
  printf( "type: " );
  zIndent( indent );
  switch( cell->type ){
  case ZM_PARSE_CELL_NUM:
    printf( "number = %s\n", zName(cell->dat.num) );    break;
  case ZM_PARSE_CELL_VAR:
    printf( "variable = %s\n", zName(cell->dat.var) );  break;
  case ZM_PARSE_CELL_OP:
    printf( "operator = %s\n", cell->dat.op->str );     break;
  case ZM_PARSE_CELL_FUNC:
    printf( "function = %s\n", zName(cell->dat.func) ); break;
  case ZM_PARSE_CELL_PAR:
    printf( "parenthesis = %s\n", cell->dat.parident ); break;
  case ZM_PARSE_CELL_CHUNK:
    printf( "[chunk] operator = %s\n", cell->dat.op->str ); break;
  default:
    printf( "invalid\n" ); break;
  }
  if( cell->arg1 ){
    zIndent( indent+2 );
    printf( "<arg1>\n" );
    _zmParseCellWrite( cell->arg1, indent+2 );
  }
  if( cell->arg2 ){
    zIndent( indent+2 );
    printf( "<arg2>\n" );
    _zmParseCellWrite( cell->arg2, indent+2 );
  }
}
#endif

/* token stack */

char _zmParseParChar(zmParseStatus status)
{
  switch( status ){
  case ZM_PARSE_IN_PAR:
  case ZM_PARSE_IN_FUNC:     return '(';
  case ZM_PARSE_OUT_PAR:
  case ZM_PARSE_OUT_FUNC:    return ')';
  case ZM_PARSE_IN_BRACE:    return '{';
  case ZM_PARSE_OUT_BRACE:   return '}';
  case ZM_PARSE_IN_BRACKET:  return '[';
  case ZM_PARSE_OUT_BRACKET: return ']';
  default: return ' ';
  }
}

void _zmParseTokFree(zmParseTok *tok)
{
  if( !tok ) return;
  if( tok->cell )
    _zmParseCellFree( tok->cell );
  zFree( tok );
}

int _zmParseTokLevel(zmParseTok *tok)
{
  return tok ? tok->level : 0;
}

zmParseStatus _zmParseTokStatus(zmParseTok *tok)
{
  return tok ? tok->status : ZM_PARSE_FLAT;
}

int _zmParseTokOpSpec(zmParseTok *tok)
{
  if( !tok )
    return ZM_PARSE_OP_INVALID;
  if( tok->cell->type == ZM_PARSE_CELL_OP )
    return tok->cell->dat.op->spec;
  return tok->prev && tok->prev->level != tok->level ?
    ZM_PARSE_OP_INVALID : _zmParseTokOpSpec(tok->prev);
}

int _zmParseTokPriority(zmParseTok *tok)
{
  if( !tok )
    return ZM_PARSE_MAX_PRIORITY;
  if( tok->cell->type == ZM_PARSE_CELL_OP )
    return tok->cell->dat.op->priority;
  return tok->prev && tok->prev->level != tok->level ?
    ZM_PARSE_MAX_PRIORITY : _zmParseTokPriority(tok->prev);
}

char *_zmParseTokStr(zmParseTok *tok)
{
  return tok ? _zmParseCellStr( tok->cell ) : zNullStr();
}

#ifdef DEBUG
void _zmParseTokWrite(zmParseTok *tok) /* for debug */
{
  printf( "status = " );
  switch( tok->status ){
  case ZM_PARSE_FLAT:        printf( "flat\n" );               break;
  case ZM_PARSE_IN_PAR:      printf( "in parenthesis\n" );     break;
  case ZM_PARSE_OUT_PAR:     printf( "out of parenthesis\n" ); break;
  case ZM_PARSE_IN_BRACE:    printf( "in brace\n" );           break;
  case ZM_PARSE_OUT_BRACE:   printf( "out of brace\n" );       break;
  case ZM_PARSE_IN_BRACKET:  printf( "in bracket\n" );         break;
  case ZM_PARSE_OUT_BRACKET: printf( "out of bracket\n" );     break;
  case ZM_PARSE_IN_FUNC:     printf( "in function\n" );        break;
  case ZM_PARSE_OUT_FUNC:    printf( "out of function\n" );    break;
  default:                   printf( "___\n" );                break;
  }
  printf( " level = %d\n", tok->level );
  printf( " op.spec. = " );
  switch( _zmParseTokOpSpec(tok) ){
  case ZM_PARSE_OP_INVALID:  printf( "invalid\n" );   break;
  case ZM_PARSE_OP_PREFIX:   printf( "prefix\n" );    break;
  case ZM_PARSE_OP_INFIX:    printf( "infix\n" );     break;
  case ZM_PARSE_OP_PREINFIX: printf( "pre-infix\n" ); break;
  case ZM_PARSE_OP_POSTFIX:  printf( "postfix\n" );   break;
  }
  printf( " priority = %d\n", _zmParseTokPriority(tok) );
  _zmParseCellWrite( tok->cell, 1 );
}
#endif

void _zmParseStackInit(zmParseTok *stack)
{
  stack->status = ZM_PARSE_FLAT;
  stack->level = 0;
  stack->cell = NULL;
  stack->prev = NULL;
}

void _zmParseStackDestroy(zmParseTok *stack)
{
  zmParseTok *tok;

  while( stack->prev ){
    tok = stack->prev;
    stack->prev = tok->prev;
    _zmParseTokFree( tok );
  }
}

#ifdef DEBUG
void _zmParseStackWrite(zmParseTok *stack) /* for debug */
{
  zmParseTok *tok;

  for( tok=stack->prev; tok; tok=tok->prev )
    _zmParseTokWrite( tok );
}
#endif

/* parser */

void zmParserInit(zmParser *parser)
{
  _zmParseStackInit( &parser->stack );
  zListInit( &parser->varlist );
  _zmParseVarBuiltin( &parser->varlist );
  zListInit( &parser->funclist );
  _zmParseFuncBuiltin( &parser->funclist );
  parser->level = 0;
  parser->status = ZM_PARSE_FLAT;
}

void zmParserDestroy(zmParser *parser)
{
  _zmParseStackDestroy( &parser->stack );
  _zmParseVarListDestroy( &parser->varlist );
  _zmParseFuncListDestroy( &parser->funclist );
}

zmParseTok *_zmParseTokAlloc(zmParser *parser, zmParseCellType type, char *str)
{
  zmParseCell *cell;
  zmParseTok *tok;

  if( !( cell = _zmParseCellAlloc( type, str, &parser->varlist, &parser->funclist ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( tok = zAlloc( zmParseTok, 1 ) ) ){
    ZALLOCERROR();
    _zmParseCellFree( cell );
    return NULL;
  }
  tok->level = parser->level;
  tok->status = parser->status;
  tok->cell = cell;
  tok->prev = NULL;
  return tok;
}

void _zmParseStackPush(zmParser *parser, zmParseTok *tok, int level, zmParseStatus status)
{
  if( !tok ) return;
  tok->level = level;
  tok->status = status;
  tok->prev = parser->stack.prev;
  parser->stack.prev = tok;
}

zmParseTok *_zmParseStackPop(zmParser *parser, zmParseTok *toknext)
{
  zmParseTok *tok;

  if( !parser->stack.prev ||
      ( toknext && parser->stack.prev->level != toknext->level ) )
    return NULL;
  tok = parser->stack.prev;
  parser->stack.prev = tok->prev;
  return tok;
}

zmParseTok *_zmParseTokNum(zmParser *parser, char *str, char *tkn, size_t size)
{
  return *zSNumToken( str, tkn, size ) ?
    _zmParseTokAlloc( parser, ZM_PARSE_CELL_NUM, tkn ) : NULL;
}

zmParseTok *_zmParseTokVar(zmParser *parser, char *str, char *tkn, size_t size)
{
  char *cp, *sp;

  cp = tkn;
  if( isalpha( *str ) ){
    for( sp=str; isalnum(*sp); *cp++=*sp++ );
    zStrCopyNC( str, sp );
  }
  *cp = '\0';
  return *tkn ? _zmParseTokAlloc( parser, ZM_PARSE_CELL_VAR, tkn ) : NULL;
}

zmParseTok *_zmParseTokOp(zmParser *parser, char *str, char *tkn, size_t size)
{
  zmParseOp *op;
  char *cp, *sp;

  for( cp=tkn, sp=str; zIsIncludedChar(*sp,__zm_parse_op); )
    *cp++ = *sp++;
  *cp = '\0';
  zStrCopyNC( str, sp );
  return ( op = _zmParseOpFind( tkn ) ) ?
    _zmParseTokAlloc( parser, ZM_PARSE_CELL_OP, tkn ) : NULL;
}

zmParseTok *_zmParseTokFunc(zmParser *parser, char *str, char *tkn, size_t size)
{
  char *cp, *sp;

  cp = tkn;
  if( isalpha( *str ) ){
    for( sp=str; isalnum(*sp); *cp++=*sp++ );
    if( *sp != '(' ) return NULL;
    parser->status = ZM_PARSE_IN_FUNC;
    parser->level++;
    zStrCopyNC( str, sp+1 );
  }
  *cp = '\0';
  return *tkn ? _zmParseTokAlloc( parser, ZM_PARSE_CELL_FUNC, tkn ) : NULL;
}

zmParseTok *_zmParseTokPar(zmParser *parser, char *str, char *tkn, size_t size)
{
  switch( *str ){
  case '(':
    parser->status = ZM_PARSE_IN_PAR;
    parser->level++;
    zStrCopyNC( str, str+1 );
    return _zmParseTokAlloc( parser, ZM_PARSE_CELL_PAR, "(" );
  case ')':
    parser->status = parser->status == ZM_PARSE_IN_FUNC ? ZM_PARSE_OUT_FUNC : ZM_PARSE_OUT_PAR;
    parser->level--;
    zStrCopyNC( str, str+1 );
    return _zmParseTokAlloc( parser, ZM_PARSE_CELL_PAR, ")" );
  case '{':
    parser->status = ZM_PARSE_IN_BRACE;
    parser->level++;
    zStrCopyNC( str, str+1 );
    return _zmParseTokAlloc( parser, ZM_PARSE_CELL_PAR, "{" );
  case '}':
    parser->status = ZM_PARSE_OUT_BRACE;
    parser->level--;
    zStrCopyNC( str, str+1 );
    return _zmParseTokAlloc( parser, ZM_PARSE_CELL_PAR, "}" );
  case '[':
    parser->status = ZM_PARSE_IN_BRACKET;
    parser->level++;
    zStrCopyNC( str, str+1 );
    return _zmParseTokAlloc( parser, ZM_PARSE_CELL_PAR, "[" );
  case ']':
    parser->status = ZM_PARSE_OUT_BRACKET;
    parser->level--;
    zStrCopyNC( str, str+1 );
    return _zmParseTokAlloc( parser, ZM_PARSE_CELL_PAR, "]" );
  default: ;
  }
  return NULL;
}

bool _zmParseStackPrefix(zmParser *parser, zmParseTok *toknext)
{
  zmParseTok *op, *arg;

  arg = _zmParseStackPop( parser, toknext );
  op = _zmParseStackPop( parser, toknext );
  if( !op || !arg || op->cell->type != ZM_PARSE_CELL_OP ){
    ZRUNERROR( "syntax error around %s %s ", _zmParseTokStr(op), _zmParseTokStr(arg) );
    _zmParseTokFree( op );
    _zmParseTokFree( arg );
    return false;
  }
  op->cell->arg1 = arg->cell;
  op->cell->type = ZM_PARSE_CELL_CHUNK;
  arg->cell = NULL;
  _zmParseTokFree( arg );
  _zmParseStackPush( parser, op, op->level, op->status );
  return true;
}

bool _zmParseStackInfix(zmParser *parser, zmParseTok *toknext)
{
  zmParseTok *op, *arg1, *arg2;

  arg2 = _zmParseStackPop( parser, toknext );
  op = _zmParseStackPop( parser, toknext );
  arg1 = _zmParseStackPop( parser, toknext );
  if( !op || !arg1 || !arg2 || op->cell->type != ZM_PARSE_CELL_OP ){
    ZRUNERROR( "syntax error around %s %s %s ", _zmParseTokStr(arg1), _zmParseTokStr(op), _zmParseTokStr(arg2) );
    _zmParseTokFree( op );
    _zmParseTokFree( arg1 );
    _zmParseTokFree( arg2 );
    return false;
  }
  op->cell->arg1 = arg1->cell;
  op->cell->arg2 = arg2->cell;
  op->cell->type = ZM_PARSE_CELL_CHUNK;
  arg1->cell = NULL;
  _zmParseTokFree( arg1 );
  arg2->cell = NULL;
  _zmParseTokFree( arg2 );
  _zmParseStackPush( parser, op, op->level, op->status );
  return true;
}

bool _zmParseStackPreinfix(zmParser *parser, zmParseTok *toknext)
{
  zmParseTok *op, *arg1, *arg2;

  arg2 = _zmParseStackPop( parser, toknext );
  op = _zmParseStackPop( parser, toknext );
  arg1 = _zmParseStackPop( parser, toknext );
  if( !op || !arg2 || op->cell->type != ZM_PARSE_CELL_OP ){
    ZRUNERROR( "syntax error around %s %s %s ", _zmParseTokStr(arg1), _zmParseTokStr(op), _zmParseTokStr(arg2) );
    _zmParseTokFree( op );
    _zmParseTokFree( arg1 );
    _zmParseTokFree( arg2 );
    return false;
  }
  op->cell->arg1 = arg1 ? arg1->cell : NULL;
  op->cell->arg2 = arg2->cell;
  op->cell->type = ZM_PARSE_CELL_CHUNK;
  if( arg1 ){
    arg1->cell = NULL;
    _zmParseTokFree( arg1 );
  }
  arg2->cell = NULL;
  _zmParseTokFree( arg2 );
  _zmParseStackPush( parser, op, op->level, op->status );
  return true;
}

bool _zmParseStackPostfix(zmParser *parser, zmParseTok *toknext)
{
  zmParseTok *op, *arg;

  op = _zmParseStackPop( parser, toknext );
  arg = _zmParseStackPop( parser, toknext );
  if( !op || !arg || op->cell->type != ZM_PARSE_CELL_OP ){
    ZRUNERROR( "syntax error around %s %s ", _zmParseTokStr(arg), _zmParseTokStr(op) );
    _zmParseTokFree( op );
    _zmParseTokFree( arg );
    return false;
  }
  op->cell->arg1 = arg->cell;
  op->cell->type = ZM_PARSE_CELL_CHUNK;
  arg->cell = NULL;
  _zmParseTokFree( arg );
  _zmParseStackPush( parser, op, op->level, op->status );
  return true;
}

bool _zmParseStack(zmParser *parser, zmParseTok *toknext)
{
  bool ret = true;

  switch( _zmParseTokOpSpec( parser->stack.prev ) ){
  case ZM_PARSE_OP_PREFIX:
    ret = _zmParseStackPrefix( parser, toknext );   break;
  case ZM_PARSE_OP_INFIX:
    ret = _zmParseStackInfix( parser, toknext );    break;
  case ZM_PARSE_OP_PREINFIX:
    ret = _zmParseStackPreinfix( parser, toknext ); break;
  case ZM_PARSE_OP_POSTFIX:
    ret = _zmParseStackPostfix( parser, toknext );  break;
  default: ;
  }
  return ret;
}

bool _zmParseIn(zmParser *parser, zmParseStatus status)
{
  zmParseTok *tok, *tokpar;

  if( _zmParseTokStatus( parser->stack.prev ) != status ){
    ZRUNERROR( "mismatch of %c and %c ", _zmParseParChar(_zmParseTokStatus(parser->stack.prev)), _zmParseParChar(status) );
    return false;
  }
  while( 1 ){
    _zmParseStack( parser, NULL );
    tok = _zmParseStackPop( parser, NULL );
    if( tok->cell->type == ZM_PARSE_CELL_PAR ){
      _zmParseTokFree( tok );
      return true;
    }
    if( !parser->stack.prev ) break;
    if( parser->stack.prev->cell->type == ZM_PARSE_CELL_PAR ){
      tokpar = _zmParseStackPop( parser, tok );
      _zmParseTokFree( tokpar );
      break;
    } else
    if( parser->stack.prev->cell->type == ZM_PARSE_CELL_FUNC ){
      tokpar = _zmParseStackPop( parser, tok );
      tokpar->cell->arg1 = tok->cell;
      tok->cell = NULL;
      _zmParseTokFree( tok );
      tok = tokpar;
      tok->level = _zmParseTokLevel( parser->stack.prev );
      tok->status = _zmParseTokStatus( parser->stack.prev );
    }
    _zmParseStackPush( parser, tok, tok->level, tok->status );
  }
  _zmParseStackPush( parser, tok, parser->level, ( parser->status = _zmParseTokStatus( parser->stack.prev ) ) );
  return true;
}

zmParseTok *_zmParsePar(zmParser *parser, zmParseTok *tok)
{
  static zmParseTok tok_par;

  switch( parser->status ){
  case ZM_PARSE_OUT_PAR:
    _zmParseTokFree( tok );
    tok = _zmParseIn( parser, ZM_PARSE_IN_PAR ) ? &tok_par : NULL;
    break;
  case ZM_PARSE_OUT_BRACE:
    _zmParseTokFree( tok );
    tok = _zmParseIn( parser, ZM_PARSE_IN_BRACE ) ? &tok_par : NULL;
    break;
  case ZM_PARSE_OUT_BRACKET:
    _zmParseTokFree( tok );
    tok = _zmParseIn( parser, ZM_PARSE_IN_BRACKET ) ? &tok_par : NULL;
    break;
  case ZM_PARSE_OUT_FUNC:
    _zmParseTokFree( tok );
    tok = _zmParseIn( parser, ZM_PARSE_IN_FUNC ) ? & tok_par : NULL;
    break;
  default:
    _zmParseStackPush( parser, tok, parser->level, parser->status );
  }
  return tok;
}

void zmParse(zmParser *parser, char *str)
{
  zmParseTok *tok;
  char tkn[BUFSIZ], *cp;

  do{
    if( ( cp = zSSkipIncludedChar( str, __zm_parse_delim ) ) )
      zStrCopyNC( str, cp );
    if( ( tok = _zmParseTokOp( parser, str, tkn, BUFSIZ ) ) ){
      while( _zmParseTokPriority( tok ) >= _zmParseTokPriority( parser->stack.prev ) )
        _zmParseStack( parser, tok );
      _zmParseStackPush( parser, tok, parser->level, parser->status );
    } else
    if( ( tok = _zmParseTokFunc( parser, str, tkn, BUFSIZ ) ) ||
        ( tok = _zmParseTokNum( parser, str, tkn, BUFSIZ ) ) ||
        ( tok = _zmParseTokVar( parser, str, tkn, BUFSIZ ) ) ){
      _zmParseStackPush( parser, tok, parser->level, parser->status );
    } else
    if( ( tok = _zmParseTokPar( parser, str, tkn, BUFSIZ ) ) ){
      tok = _zmParsePar( parser, tok );
    } else
      tok = NULL;
  } while( tok );
  while( _zmParseTokOpSpec( parser->stack.prev ) != ZM_PARSE_OP_INVALID )
    _zmParseStack( parser, NULL );
}

double zmParseEval(zmParser *parser)
{
  zmParseTok *tok;
  double ret;

  if( !( tok = _zmParseStackPop( parser, NULL ) ) )
    return 0;
  if( parser->stack.prev ){
    ZRUNERROR( "syntax error around %s %s ", _zmParseTokStr(parser->stack.prev), _zmParseTokStr(tok) );
    return 0;
  }
  ret = _zmParseCellEval( tok->cell );
  _zmParseTokFree( tok );
  return ret;
}
