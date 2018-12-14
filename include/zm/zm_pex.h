/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_pex - polynomial expression class.
 */

#ifndef __ZM_PEX_H__
#define __ZM_PEX_H__

#include <zm/zm_stat.h>
#include <zm/zm_le.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zPex
 * polynomial expression class
 * ********************************************************** */

typedef zVec zPex;

/* n-dim polynomial has n+1 coeficients including the constant term. */
#define zPexDim(p)               (int)( zVecSizeNC(p) - 1 )
#define zPexCoeff(p,i)           zVecElem(p,i)

#define zPexCoeffHigh(p,i)       zPexCoeff(p,zPexDim(p)-(i))
#define zPexCoeffLow(p,i)        zPexCoeff(p,i)

#define zPexSetDim(p,d)          zVecSetSize(p,(d)+1)
#define zPexSetCoeff(p,i,c)      zVecSetElem(p,i,c)

#define zPexSetCoeffHigh(p,i,c)  zPexSetCoeff(p,zPexDim(p)-(i),c)
#define zPexSetCoeffLow(p,i,c)   zPexSetCoeff(p,i,c)

#define zPexIsEqual(p1,p2)       zVecIsEqual(p1,p2)

#define zPexTouchup(p)           zVecTouchup(p)

#define zPexAlloc(n)             zVecAlloc((n)+1)
#define zPexFree(p)              zVecFree(p)
#define zPexCopy(s,d)            zVecCopy(s,d)
#define zPexCopyArray(a,n,p)     zVecCopyArray(a,(n)+1,p)
#define zPexClone(p)             zVecClone(p)
#define zPexCloneArray(a,n)      zVecCloneArray(a,(n)+1)

#define zPexSetCoeffList(p,a...) zVecSetElemList(p,##a)
#define zPexCreateList(s,a...)   zVecCreateList((s)+1,##a)

/*! \brief regulate polynomial.
 *
 * zPexRgl() regulates the polynomial expression \a p.
 * The dimension of \a p is reduced so that the coefficient
 * of the greatest dimension is non-zero.
 * \return
 * zPexRgl() returns a pointer \a p.
 */
__EXPORT zPex zPexRgl(zPex *p);

/*! \brief direct addision and substraction of polynomial expression.
 *
 * zPexAddDRC() adds the polynomial expression \a p2 directly
 * to \a p1.
 *
 * zPexSubDRC() substracts the polynomial expression \a p2
 * directly from \a p1.
 * \return
 * Both zPexAddDRC() and zPexSubDRC() return a pointer \a p1.
 * \notes
 * The dimension of \a p1 must not be less than that of \a p2.
 */
__EXPORT zPex zPexAddDRC(zPex p1, zPex p2);
__EXPORT zPex zPexSubDRC(zPex p1, zPex p2);

/* METHOD:
 * zPexAdd, zPexSub, zPexMul, zPexDiv
 * - arithmatic of polynomial expression.
 * [SYNOPSIS]
 * zPex zPexAdd(zPex p1, zPex p2);
 * zPex zPexSub(zPex p1, zPex p2);
 * zPex zPexMul(zPex p1, zPex p2);
 * bool zPexDiv(zPex p, zPex f, zPex *q, zPex *r);
 * [DESCRIPTION]
 * 'zPexAdd()' adds the two polynomial expression 'p1'
 * and 'p2', and returns the pointer to the result, which is
 * newly allocated.
 *
 * 'zPexSub()' substracts the polynomial expression 'p2'
 * from 'p1', and returns the pointer to the result, which
 * is newly allocated.
 *
 * 'zPexMul()' multiplies the polynomial expression 'p1'
 * by 'p2', and returns the pointer to the result, which is
 * newly allocated.
 *
 * 'zPexDiv()' divides the polynomial expression 'p'
 * by 'f'. The quotient is set where is pointed by 'q', and
 * the surplus is set whare is pointed by 'r'.
 * [RETURN VALUE]
 * 'zPexAdd()', 'zPexSub()' and 'zPexMul()' return
 * pointers to the results.
 *
 * 'zPexDiv()' returns a boolean value. If it fails to
 * allocate working space, the false value is returned.
 * Otherwise, the true value is returned.
 */
__EXPORT zPex zPexAdd(zPex p1, zPex p2);
__EXPORT zPex zPexSub(zPex p1, zPex p2);
__EXPORT zPex zPexMul(zPex p1, zPex p2);
__EXPORT bool zPexDiv(zPex p, zPex f, zPex *q, zPex *r);

/* METHOD:
 * zPexExp
 * - expand factors into polynomial expression.
 * [SYNOPSIS]
 * zPex zPexExp(zVec factor);
 * [DESCRIPTION]
 * 'zPexExp()' expands the factorial style into a
 * polynomial expression.
 * When 'factor' is [a1 a2 ... aN], for example,
 * it means that the factorial style (x-a1)(x-a2)...(x-aN)
 * is given.
 * The result newly allocated is returned.
 * [RETURN VALUE]
 * 'zPexExp()' returns the pointer to the result.
 */
__EXPORT zPex zPexExp(zVec factor);

/* METHOD:
 * zPexModulo - modulo of a primary expression.
 * [SYNOPSIS]
 * zPex zPexModulo(zPex p1, double a, zPex p2);
 * [DESCRIPTION]
 * 'zPexModulo()' rearranges a given polynomial
 * expression 'p1' modulo a primary expression x+'a'.
 * Namely:
 *     'p1': c_0 + c_1 x^1 + ... + c_n x^n
 *  -> 'p2': b_0 + b_1(x+'a')^1 + ... + b_n (x+'a')^n
 * where 'p1' and 'p2' are equivalent.
 * The resultant series of coefficients are stored
 * in 'p2'. Note that 'p2' is that modulo x+'a',
 * so that it should not be manipulated by ordinary
 * PEX methods.
 * [NOTES]
 * 'p1' and 'p2' must be in the same dimension.
 * [RETURN VALUE]
 * 'zPexModulo()' returns a pointer 'p2' if 'p2'
 * has the same dimension with 'p1'. Otherwise,
 * the null pointer is returned.
 */
__EXPORT zPex zPexModulo(zPex p1, double a, zPex p2);

/* METHOD:
 * zPexDif, zPexIntg
 * - differentiation and integration of polynomial expression.
 * [SYNOPSIS]
 * zPex zPexDif(zPex p);
 * zPex zPexIntg(zPex p);
 * [DESCRIPTION]
 * 'zPexDif()' differentiates the polynoimal expression
 * 'p', and returns the pointer to the result, which
 * is newly allocated.
 *
 * 'zPexIntg()' integrates the polynoimal
 * expression 'p'. The integration constant will be zero
 * temporary. And it returns the pointer to the result,
 * which is newly allocated.
 * [RETURN VALUE]
 * All these functions return pointers to the results.
 */
__EXPORT zPex zPexDif(zPex p);
__EXPORT zPex zPexIntg(zPex p);

/* METHOD:
 * zPexVal, zPexCVal, zPexDifVal
 * - evaluate polynomial expression for one argument.
 * [SYNOPSIS]
 * double zPexVal(zPex p, double arg);
 * zComplex *zPexCVal(zPex p, zComplex *arg, zComplex *c);
 * double zPexDifVal(zPex p, int dim, double arg);
 * [DESCRIPTION]
 * 'zPexVal()' calculates the value of polynomial
 * expression 'p', giving an argument 'arg'.
 *
 * 'zPexCVal()' calculates the value of polynomial
 * expression 'p', giving a complex number argument 'arg'.
 * The result is stored in 'c'.
 *
 * 'zPexDifVal()' calculates the 'dim'th differential
 * value of polynomial expression 'p', giving the argument
 * 'arg'.
 * When 'dim' == 0, 'zPexDivVal( p, 0, x )' returns
 * the same value with 'zPexVal( p, x )'.
 * [RETURN VALUE]
 * 'zPexVal()' and 'zPexDifVal()' returns the value
 * calculated.
 * 'zPexCVal()' returns a pointer 'c'.
 */
__EXPORT double zPexVal(zPex p, double arg);
__EXPORT zComplex *zPexCVal(zPex p, zComplex *arg, zComplex *c);
__EXPORT double zPexDifVal(zPex p, int dim, double arg);

/* zPexFWrite, zPexWrite, zPexFExpr, zPexExpr
 * - expression of polynomial expression.
 *
 * 'zPexFWrite()' writes the coefficients of the polynomial
 * expression 'p' in ascending absolute value order to the
 * current position of the file pointed by 'fp'.
 * And 'zPexWrite()' writes it simply to the standerd output.
 *
 * 'zPexFExpr()' writes the expression 'p' in decending
 * absolute value order in the expression form to the current
 * position of file 'fp'. 'zPexExpr()' writes the expression
 * 'p' to the standerd output.
 * [RETURN VALUE]
 * 'zPexFWrite()', 'zPexWrite()', 'zPexFExpr()'
 * and 'zPexExpr()' return no values.
 */
__EXPORT zPex zPexFRead(FILE *fp);
#define zPexRead(p)     zPexFRead( stdin, p )
__EXPORT void zPexFWrite(FILE *fp, zPex p);
#define zPexWrite(p)    zPexFWrite( stdout, p )
__EXPORT void zPexFExpr(FILE *fp, zPex p, char c);
#define zPexExpr(p,c)   zPexFExpr( stdout, (p), (c) )
#define zPexFExprX(f,p) zPexFExpr( (f), (p), 'x' )
#define zPexExprX(p)    zPexExpr( p, 'x' )

__END_DECLS

#include <zm/zm_pex_eq.h> /* polynomial equation solver */

#endif /* __ZM_PEX_H__ */
