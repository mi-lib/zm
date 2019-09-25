/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_pex - polynomial expression class.
 */

#ifndef __ZM_PEX_H__
#define __ZM_PEX_H__

#include <zm/zm_stat.h>
#include <zm/zm_cvec.h>
#include <zm/zm_le.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zPex
 * polynomial expression class
 * ********************************************************** */

typedef zVec zPex;

/* n-dim polynomial has n+1 coeficients including the constant term. */
#define zPexDim(p)               (int)( zVecSizeNC(p) - 1 )
#define zPexCoeff(p,i)           zVecElemNC(p,i)

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

/*! \brief arithmatic of polynomial expression.
 *
 * zPexAdd() adds two polynomial expressions \a p1 and \a p2,
 * and returns a pointer to the result, which is newly
 * allocated.
 *
 * zPexSub() substracts a polynomial expression \a p2 from
 * \a  p1, and returns a pointer to the result, which is
 * newly allocated.
 *
 * zPexMul() multiplies a polynomial expression \a p1 by
 * \a p2, and returns a pointer to the result, which is
 * newly allocated.
 *
 * zPexDiv() divides a polynomial expression \a p by \a f.
 * The quotient is set where is pointed by \a q, and the
 * surplus is set whare is pointed by \a r.
 * \return
 * zPexAdd(), zPexSub() and zPexMul() return a pointer to
 * the result.
 *
 * zPexDiv() returns a boolean value. If it fails to allocate
 * memory, the false value is returned. Otherwise, the true
 * value is returned.
 */
__EXPORT zPex zPexAdd(zPex p1, zPex p2);
__EXPORT zPex zPexSub(zPex p1, zPex p2);
__EXPORT zPex zPexMul(zPex p1, zPex p2);
__EXPORT bool zPexDiv(zPex p, zPex f, zPex *q, zPex *r);

/*! \brief expand factors into a polynomial expression.
 *
 * zPexExp() and zPexCExp() expand a factorial style into a polynomial
 * expression. Suppose \a factor is [a1 a2 ... aN], the corresponding
 * factorial style is (x-a1)(x-a2)...(x-aN).
 *
 * zPexExp() accepts real-number factors given by a real vector, while
 * zPexCExp() does complex-number factors by a complex vector.
 *
 * The result is newly allocated for both functions.
 * \note
 * The factors given to zPexCExp() have to include only real numbers
 * and paired co-conjugate imaginary numbers.
 * \return
 * zPexExp() and zPexCExp() return a pointer to the result.
 */
__EXPORT zPex zPexExp(zVec factor);
__EXPORT zPex zPexCExp(zCVec factor);

/*! \brief compute modulo of a primary expression.
 *
 * zPexModulo() rearranges a given polynomial expression
 * \a p1 modulo a primary expression x+\a a.
 * Namely:
 *     \a p1: c_0 + c_1 x^1 + ... + c_n x^n
 *  -> \a p2: b_0 + b_1(x+\a a)^1 + ... + b_n (x+\a a)^n,
 * where \a p1 and \a p2 are equivalent.
 *
 * The resultant series of coefficients are stored in
 * \a p2. Note that \a p2 is that modulo x+\a a, so that
 * it should not be manipulated by ordinary PEX methods.
 * \notes
 * \a p1 and \a p2 must be in the same dimension.
 * \return
 * zPexModulo() returns a pointer \a p2 if \a p2 has the
 * same dimension with \a p1. Otherwise, the null pointer
 * is returned.
 */
__EXPORT zPex zPexModulo(zPex p1, double a, zPex p2);

/*! \brief differentiate and integrate a polynomial expression.
 *
 * zPexDif() differentiates a polynoimal expression \a p,
 * and returns a pointer to the result, which is newly
 * allocated.
 *
 * zPexIntg() integrates a polynoimal expression \a p.
 * The integration constant will be zero temporary. It
 * returns a pointer to the result, which is newly allocated.
 * \return
 * These functions return a pointer to the result.
 */
__EXPORT zPex zPexDif(zPex p);
__EXPORT zPex zPexIntg(zPex p);

/*! \brief evaluate a polynomial expression for one argument.
 *
 * zPexVal() calculates the value of a polynomial
 * expression \a p, given an argument \a arg.
 *
 * zPexCVal() calculates the value of a polynomial
 * expression \a p, given a complex number argument
 * \a arg. The result is stored in \a c.
 *
 * zPexDifVal() calculates the \a dim th differential
 * value of a polynomial expression \a p, given the
 * argument \a arg.
 * If \a dim is 0, zPexDivVal( \a p, 0, \a x ) returns
 * the same value with zPexVal( \a p, \a x ).
 * \return
 * zPexVal() and zPexDifVal() returns the value
 * calculated.
 * zPexCVal() returns a pointer \a c.
 */
__EXPORT double zPexVal(zPex p, double arg);
__EXPORT zComplex *zPexCVal(zPex p, zComplex *arg, zComplex *c);
__EXPORT double zPexDifVal(zPex p, int dim, double arg);

/*! \brief scan a polynomial expression from a ZTK processor. */
__EXPORT zPex zPexFromZTK(ZTK *ztk);

/*! \brief expression of a polynomial expression.
 *
 * zPexFScan() scans coefficients of a polynomial
 * expression \a p in ascending absolute value order
 * from the current position of a file \a fp.
 * zPexScan() scans them from the standerd input.
 *
 * zPexFPrint() prints coefficients of a polynomial
 * expression \a p in ascending absolute value order
 * to the current position of a file \a fp.
 * zPexPrint() prints them to the standerd output.
 *
 * zPexFExpr() prints the expression \a p in decending
 * absolute value order in the expression form to the
 * current position of a file a fp. zPexExpr() prints
 * the expression \a p to the standerd output.
 * \return
 * zPexFPrint(), zPexPrint(), zPexFExpr() and zPexExpr()
 * return no values.
 */
__EXPORT zPex zPexFScan(FILE *fp);
#define zPexScan(p)     zPexFScan( stdin, p )
__EXPORT void zPexFPrint(FILE *fp, zPex p);
#define zPexPrint(p)    zPexFPrint( stdout, p )
__EXPORT void zPexFExpr(FILE *fp, zPex p, char c);
#define zPexExpr(p,c)   zPexFExpr( stdout, (p), (c) )
#define zPexFExprX(f,p) zPexFExpr( (f), (p), 'x' )
#define zPexExprX(p)    zPexExpr( p, 'x' )

__END_DECLS

#include <zm/zm_pex_eq.h> /* polynomial equation solver */

#endif /* __ZM_PEX_H__ */
