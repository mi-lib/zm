/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_parse.h
 * \brief parser of mathematical expression.
 * \author Zhidao
 */

#ifndef __ZM_PARSE_H__
#define __ZM_PARSE_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/* pre-declaration of zmParseCell class */

struct _zmParseCell;

/* number */

typedef struct{
  Z_NAMED_CLASS
  double val;
} zmParseNum;

/* list of variables */

typedef struct{
  Z_NAMED_CLASS
  struct _zmParseCell *cell;
} zmParseVar;
zListClass( zmParseVarList, zmParseVarCell, zmParseVar );

/* list of functions */

typedef enum{
  ZM_PARSE_FUNC_SIN,
  ZM_PARSE_FUNC_COS,
  ZM_PARSE_FUNC_TAN,
  ZM_PARSE_FUNC_COT,
  ZM_PARSE_FUNC_ASIN,
  ZM_PARSE_FUNC_ACOS,
  ZM_PARSE_FUNC_ATAN,
  ZM_PARSE_FUNC_ACOT,
  ZM_PARSE_FUNC_SINH,
  ZM_PARSE_FUNC_COSH,
  ZM_PARSE_FUNC_TANH,
  ZM_PARSE_FUNC_COTH,
  ZM_PARSE_FUNC_EXP,
  ZM_PARSE_FUNC_LOG,
  ZM_PARSE_FUNC_LN,
  ZM_PARSE_FUNC_ABS,
  ZM_PARSE_FUNC_SQRT,
  ZM_PARSE_FUNC_DEG,
  ZM_PARSE_FUNC_RAD,
  ZM_PARSE_FUNC_USER,
} zmParseFuncID;

typedef struct{
  Z_NAMED_CLASS
  zmParseFuncID id;
  int argc;
  struct _zmParseCell *cell;
} zmParseFunc;
zListClass( zmParseFuncList, zmParseFuncCell, zmParseFunc );

/* operator */

typedef enum{
  ZM_PARSE_OP_VOID=0,
  ZM_PARSE_OP_ADD,
  ZM_PARSE_OP_SUB,
  ZM_PARSE_OP_MUL,
  ZM_PARSE_OP_DIV,
  ZM_PARSE_OP_MOD,
  ZM_PARSE_OP_POW,
  ZM_PARSE_OP_FCT,
  ZM_PARSE_OP_COMMA,
  ZM_PARSE_OP_SEP,
  ZM_PARSE_OP_RAT,
  ZM_PARSE_OP_LET,
} zmParseOpType;

typedef enum{
  ZM_PARSE_OP_INVALID=-1,
  ZM_PARSE_OP_PREFIX,
  ZM_PARSE_OP_INFIX,
  ZM_PARSE_OP_PREINFIX,
  ZM_PARSE_OP_POSTFIX,
} zmParseOpSpec;

#define ZM_PARSE_MAX_PRIORITY 999

typedef struct{
  const char *str;
  zmParseOpType type;
  zmParseOpSpec spec;
  int priority;
} zmParseOp;
static const zmParseOp zm_parse_op[] = {
  { "+", ZM_PARSE_OP_ADD,   ZM_PARSE_OP_PREINFIX, 40 },
  { "-", ZM_PARSE_OP_SUB,   ZM_PARSE_OP_PREINFIX, 40 },
  { "*", ZM_PARSE_OP_MUL,   ZM_PARSE_OP_INFIX,    30 },
  { "/", ZM_PARSE_OP_DIV,   ZM_PARSE_OP_INFIX,    30 },
  { "%", ZM_PARSE_OP_MOD,   ZM_PARSE_OP_INFIX,    30 },
  { "^", ZM_PARSE_OP_POW,   ZM_PARSE_OP_INFIX,    20 },
  { "!", ZM_PARSE_OP_FCT,   ZM_PARSE_OP_POSTFIX,  10 },
  { ",", ZM_PARSE_OP_COMMA, ZM_PARSE_OP_INFIX,    60 },
  { ";", ZM_PARSE_OP_SEP,   ZM_PARSE_OP_INFIX,    80 },
  { ":", ZM_PARSE_OP_RAT,   ZM_PARSE_OP_INFIX,    50 },
  { "=", ZM_PARSE_OP_LET,   ZM_PARSE_OP_INFIX,    70 },
  { NULL, ZM_PARSE_OP_VOID, ZM_PARSE_OP_INVALID,  ZM_PARSE_MAX_PRIORITY },
};

/* cell */

typedef enum{
  ZM_PARSE_CELL_INVALID = -1,
  ZM_PARSE_CELL_NUM,
  ZM_PARSE_CELL_VAR,
  ZM_PARSE_CELL_OP,
  ZM_PARSE_CELL_FUNC,
  ZM_PARSE_CELL_PAR,
  ZM_PARSE_CELL_CHUNK, /* partially parsed cell */
} zmParseCellType;

typedef struct _zmParseCell{
  zmParseCellType type;
  union{
    zmParseNum *num;
    zmParseVar *var;
    zmParseOp *op;
    zmParseFunc *func;
    char parident[2];
  } dat;
  struct _zmParseCell *arg1;
  struct _zmParseCell *arg2;
} zmParseCell;

/* token stack */

typedef enum{
  ZM_PARSE_NOCHANGE = -1,
  ZM_PARSE_FLAT,
  ZM_PARSE_IN_PAR,
  ZM_PARSE_OUT_PAR,
  ZM_PARSE_IN_BRACE,
  ZM_PARSE_OUT_BRACE,
  ZM_PARSE_IN_BRACKET,
  ZM_PARSE_OUT_BRACKET,
  ZM_PARSE_IN_FUNC,
  ZM_PARSE_OUT_FUNC,
} zmParseStatus;

typedef struct _zmParseTok{
  zmParseStatus status;
  int level;
  zmParseCell *cell;
  struct _zmParseTok *prev;
} zmParseTok;

/* parser */

typedef struct{
  zmParseTok stack;
  zmParseVarList varlist;
  zmParseFuncList funclist;
  int level;
  zmParseStatus status;
} zmParser;

void zmParserInit(zmParser *parser);
void zmParserDestroy(zmParser *parser);
void zmParse(zmParser *parser, char *str);
double zmParseEval(zmParser *parser);

__END_DECLS

#endif /* __ZM_PARSE_H__ */
