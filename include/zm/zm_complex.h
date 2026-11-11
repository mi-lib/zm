/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_complex.h
 * \brief complex number class.
 * \author Zhidao
 */

#ifndef __ZM_COMPLEX_H__
#define __ZM_COMPLEX_H__

#include <zm/zm_misc.h>

/* ********************************************************** */
/*! \defgroup complex complex number.
 * \{ *//* ************************************************** */

/*! \struct zComplex
 * \brief complex number class.
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zComplex ){
  double re; /*!< \brief real part */
  double im; /*!< \brief imaginary part */
#ifdef __cplusplus
 public:
  zComplex() : re{0}, im{0} {}
  zComplex(double _re, double _im) : re{_re}, im{_im} {}
  ~zComplex(){}
  zComplex &create(double _re, double _im);
  zComplex &createPolar(double r, double angle);
  zComplex &operator=(double src);
  zComplex &zero();
  zComplex &touchup(double tol);
  bool isTol(double tol);
  bool isReal(double tol);
  bool isCoconj(const zComplex &anotherval);
  bool operator==(const zComplex &anotherval);
  bool operator!=(const zComplex &anotherval);
  double sqrAbs();
  double abs();
  double arg();
  zComplex conj();
  zComplex operator-();
  zComplex &operator+=(const zComplex &c);
  zComplex &operator-=(const zComplex &c);
  zComplex &operator*=(double k);
  zComplex &operator*=(const zComplex &c);
  zComplex &operator/=(double k);
  zComplex &operator/=(const zComplex &c);
  zComplex pow(double z);
  zComplex pow(const zComplex &z);
  zComplex pow(double z, const zComplex &ref);
  zComplex log(double base);
  zComplex log(const zComplex &base);
  zComplex normalize();
  zComplex &fromString(const char *str);
  zComplex &fromZTK(ZTK *ztk);
  friend zComplex operator+(const zComplex &c1, const zComplex &c2);
  friend zComplex operator-(const zComplex &c1, const zComplex rhs);
  friend zComplex operator*(const zComplex &c, double k);
  friend zComplex operator*(double k, const zComplex &c);
  friend zComplex operator*(const zComplex &c1, const zComplex &c2);
  friend zComplex operator/(const zComplex &c, double k);
  friend zComplex operator/(double k, const zComplex &c);
  friend zComplex operator/(const zComplex &c1, const zComplex &c2);
  friend std::ostream &operator<<(std::ostream &stream, const zComplex &complex);
#endif /* __cplusplus */
};

__BEGIN_DECLS

/*! \brief create a complex number.
 *
 * zComplexCreate() creates a complex number whose real part
 * is \a r and imaginary part is \a i. The result is put
 * where \a c points.
 *
 * _zComplexCreate() evaluates both \a r and \a i before
 * assigning them, so that the macros based on it work even if
 * \a c points the same with an operand.
 */
#define _zComplexCreate(c,r,i) do{\
  double __z_complex_re = (r), __z_complex_im = (i);\
  (c)->re = __z_complex_re;\
  (c)->im = __z_complex_im;\
} while(0)
__ZM_EXPORT zComplex *zComplexCreate(zComplex *c, double r, double i);

/*! \brief create a complex number from the polar expression.
 *
 * zComplexCreatePolar() creates a complex number from the polar
 * expression with radius \a r and argument angle \a t
 * in Gaussian plane, where \a t is given in radians.
 * The result is put where \a c points.
 */
#define _zComplexCreatePolar(c,r,t) _zComplexCreate( c, (r)*cos(t), (r)*sin(t) )
__ZM_EXPORT zComplex *zComplexCreatePolar(zComplex *c, double r, double t);

/*! \brief copy a complex number to another.
 *
 * zComplexCopy() copies a complex number \a src to another \a dest. */
#define _zComplexCopy(src,dest) ( *(dest) = *(src) )
__ZM_EXPORT zComplex *zComplexCopy(const zComplex *src, zComplex *dest);

/*! \brief create a zero complex number.
 *
 * zComplexZero() creates a zero complex number \a c by
 * setting both real and imaginary parts for zeros. */
#define _zComplexZero(c) _zComplexCreate(c,0,0)
#define zComplexZero(c) zComplexCreate(c,0,0)

/*! \brief touchup a complex number.
 *
 * zComplexTouchup() replaces real part or imaginary part of a complex number \a c for zero if either
 * value relative to the other part is smaller than \a tol.
 * \return
 * zComplexTouchup() returns a pointer \a c.
 */
__ZM_EXPORT zComplex *zComplexTouchup(zComplex *c, double tol);

/*! \brief test if a complex number is under the tolerance.
 *
 * zComplexIsTol() tests if the absolute values of
 * both real and imaginary parts of a given complex
 * number \a c is less than the tolerance \a tol,
 * where \a tol must be positive.
 * The result is returned as a boolean value.
 */
#define zComplexIsTol(c,tol) ( zIsTol((c)->re,tol) && zIsTol((c)->im,tol) )

/*! \brief test if a complex nubmer is under the default tolerance. */
#define zComplexIsTiny(c)    zComplexIsTol( c, zTOL )

/*! \brief check if a complex number is a real number. */
#define zComplexIsReal(c,tol) zIsTol( (c)->im, tol )

/*! \brief check if two complex numbers are equal. */
#define zComplexEqual(c1,c2,tol) ( zEqual( (c1)->re, (c2)->re, tol ) && zEqual( (c1)->im, (c2)->im, tol ) )

/*! \brief check if two complex numbers are exactly the same. */
#define zComplexMatch(c1,c2) ( (c1)->re == (c2)->re && (c1)->im == (c2)->im )

/*! \brief check if two complex numbers are co-conjugate. */
#define zComplexCoconj(c1,c2,tol) ( zEqual( (c1)->re, (c2)->re, tol ) && zEqual( (c1)->im, -(c2)->im, tol ) )

/* arithmatics */

/*! \brief squared absolute value of a complex number
 * \retval the squared absolute value of a complex number \a c.
 */
#define _zComplexSqrAbs(c) ( zSqr((c)->re) + zSqr((c)->im) )
__ZM_EXPORT double zComplexSqrAbs(const zComplex *c);

/*! \brief absolute value of a complex number
 * \retval the absolute value of a complex number \a c.
 */
#define _zComplexAbs(c) sqrt( _zComplexSqrAbs(c) )
__ZM_EXPORT double zComplexAbs(const zComplex *c);

/*! \brief argument angle of a complex number
 * \retval the argument angle of a complex number \a c,
 * which is in radians between -PI and PI.
 */
#define _zComplexArg(c) atan2( (c)->im, (c)->re )
__ZM_EXPORT double zComplexArg(const zComplex *c);

/*! \brief the complex conjugate of a complex number.
 * zComplexConj() computes the complex conjugate of a complex
 * number \a c. The result is stored where is pointed by \a cc.
 * \retval a pointer \a cc.
 */
#define _zComplexConj(c,cc) _zComplexCreate( cc, (c)->re, -(c)->im )
__ZM_EXPORT zComplex *zComplexConj(const zComplex *c, zComplex *cc);

/*! \brief arithmatics of complex numbers.
 *
 * zComplexAdd() adds two complex numbers \a c1 and \a c2.
 * The result is put into \a c.
 *
 * zComplexSub() subtracts a complex number \a c2 from another
 * \a c1. The result is put into \a c.
 *
 * zComplexRev() reverses a complex number \a c.
 * The result is put into \a rc.
 *
 * zComplexMul() multiplies a complex number \a c by a real
 * number \a k. The result is put into \a ec.
 *
 * zComplexDiv() divides a complex number \a c by a real number
 * \a k. The result is put into \a rc.
 *
 * zComplexCMul() multiplies a complex number \a c1 by another
 * \a c2. The result is put into \a c.
 *
 * zComplexCMulConj() multiplies a complex number \a c1 by the
 * conjugate of another \a c2. The result is put into \a c.
 *
 * zComplexCDiv() divides a complex number \a c1 by \another
 * \a c2. The result is put into \a c.
 *
 * _zComplexAdd(), _zComplexSub(), _zComplexRev(), _zComplexMul(),
 * _zComplexDiv(), _zComplexCMul(), and _zComplexCMulConj() are
 * macro versions of the above functions, which do not return
 * values. _zComplexDiv() does not check if \a k is zero.
 * \return
 * These functions returns a pointer to the result.
 *
 * zComplexDiv() returns the null pointer if \a k is zero.
 */
#define _zComplexAdd(c1,c2,c)      _zComplexCreate( c, (c1)->re+(c2)->re, (c1)->im+(c2)->im )
#define _zComplexSub(c1,c2,c)      _zComplexCreate( c, (c1)->re-(c2)->re, (c1)->im-(c2)->im )
#define _zComplexRev(c,rc)         _zComplexCreate( rc, -(c)->re, -(c)->im )
#define _zComplexMul(c,k,ec)       _zComplexCreate( ec, (c)->re*(k), (c)->im*(k) )
#define _zComplexDiv(c,k,ec)       _zComplexCreate( ec, (c)->re/(k), (c)->im/(k) )
#define _zComplexCMul(c1,c2,c)     _zComplexCreate( c, (c1)->re*(c2)->re - (c1)->im*(c2)->im, (c1)->re*(c2)->im + (c1)->im*(c2)->re )
#define _zComplexCMulConj(c1,c2,c) _zComplexCreate( c, (c1)->re*(c2)->re + (c1)->im*(c2)->im, -(c1)->re*(c2)->im + (c1)->im*(c2)->re )

__ZM_EXPORT zComplex *zComplexAdd(const zComplex *c1, const zComplex *c2, zComplex *c);
__ZM_EXPORT zComplex *zComplexSub(const zComplex *c1, const zComplex *c2, zComplex *c);
__ZM_EXPORT zComplex *zComplexRev(const zComplex *c, zComplex *rc);
__ZM_EXPORT zComplex *zComplexMul(const zComplex *c, double k, zComplex *ec);
__ZM_EXPORT zComplex *zComplexDiv(const zComplex *c, double k, zComplex *rc);
__ZM_EXPORT zComplex *zComplexCMul(const zComplex *c1, const zComplex *c2, zComplex *c);
__ZM_EXPORT zComplex *zComplexCMulConj(const zComplex *c1, const zComplex *c2, zComplex *c);
/*! \brief invert a complex number. */
__ZM_EXPORT zComplex *zComplexInv(const zComplex *c, zComplex *ic);
__ZM_EXPORT zComplex *zComplexCDiv(const zComplex *c1, const zComplex *c2, zComplex *c);

#define _zComplexAddDRC(c1,c2) _zComplexAdd( c1, c2, c1 )
#define _zComplexSubDRC(c1,c2) _zComplexSub( c1, c2, c1 )
#define _zComplexRevDRC(c)     _zComplexRev( c, c )
#define _zComplexMulDRC(c,k)   _zComplexMul( c, k, c )
#define _zComplexDivDRC(c,k)   _zComplexDiv( c, k, c )

#define zComplexAddDRC(c1,c2)  zComplexAdd( c1, c2, c1 )
#define zComplexSubDRC(c1,c2)  zComplexSub( c1, c2, c1 )
#define zComplexRevDRC(c)      zComplexRev( c, c )
#define zComplexMulDRC(c,k)    zComplexMul( c, k, c )
#define zComplexDivDRC(c,k)    zComplexDiv( c, k, c )

/*! \brief multiply two complex numbers. */
__ZM_EXPORT zComplex *zComplexCMulDRC(zComplex *c1, const zComplex *c2);
/*! \brief multiply a complex numbers by the conjugate of another complex number. */
__ZM_EXPORT zComplex *zComplexCMulConjDRC(zComplex *c1, const zComplex *c2);
/*! \brief divide a complex number by another. */
__ZM_EXPORT zComplex *zComplexCDivDRC(zComplex *c1, const zComplex *c2);

/*! \brief power and logarithm of complex numbers.
 *
 * zComplexPow() calculates a complex number \a c raised to the
 * power of a real number \a z. The result is put into \a pc.
 *
 * zComplexCPow() calculates a complex number \a c raised to the
 * power of another complex number \a z. The result is put into
 * \a pc.
 *
 * zComplexLog() calculates the base \a base (a real number)
 * logarithm of a complex number \a c. The result is put into
 * \a lc.
 *
 * zComplexCLog() calculates the base \a base (a complex number)
 * logarithm of another complex number \a c. The result is put
 * into \a lc.
 * \return
 * These functions returns a pointer to the result.
 * \notes
 * A number raised to the power of a complex number and a logarithm
 * of a complex number are indeed not unique because one complex
 * has non-unique argument angles.
 * zComplexPow(), zComplexCPow(), zComplexLog() and zComplexCLog()
 * calculate only one solution with the argument angle of \a c
 * between -pi and pi.
 */
__ZM_EXPORT zComplex *zComplexPow(const zComplex *c, double z, zComplex *pc);
__ZM_EXPORT zComplex *zComplexPowRef(const zComplex *c, double z, const zComplex *ref, zComplex *pc);
__ZM_EXPORT zComplex *zComplexCPow(const zComplex *c, const zComplex *z, zComplex *pc);
__ZM_EXPORT zComplex *zComplexLog(const zComplex *c, double base, zComplex *lc);
__ZM_EXPORT zComplex *zComplexCLog(const zComplex *c, const zComplex *base, zComplex *lc);

/*! \brief normalize a complex number.
 *
 * zComplexNormalize() normalizes a complex number \a c so as
 * to magnify its absolute value to be 1.
 * The result is stored where is pointed by \a nc.
 * \retval a pointer \a nc
 * \retval the null pointer if the absolute value of \a c is 0
 */
__ZM_EXPORT zComplex *zComplexNormalize(const zComplex *c, zComplex *nc);

/* I/O */

/*! \brief read a complex number from a string. */
__ZM_EXPORT zComplex *zComplexFromStr(zComplex *c, const char *str);

/*! \brief read a complex number from a ZTK format processor. */
__ZM_EXPORT zComplex *zComplexFromZTK(zComplex *c, ZTK *ztk);

/*! \brief primt a complex number.
 *
 * zComplexFPrint() prints a complex number \a c to the
 * current position of a file \a fp in the following style:
 *  x + y i
 */
__ZM_EXPORT void zComplexFPrint(FILE *fp, const zComplex *c);
/*! \brief print a complex number to the standard output. */
#define zComplexPrint(c) zComplexFPrint( stdout, (c) )

/*! \brief print the coordinates of a complex number.
 *
 * zComplexCoordFPrint() prints the coordinates of a
 * complex number \a c on Gaussian plane to the current
 * position of \a file \a fp in the following style.
 *  x y
 */
__ZM_EXPORT void zComplexCoordFPrint(FILE *fp, const zComplex *c);
/*! \brief prints the coordinates of a complex number
 * to the standard output. */
#define zComplexCoordPrint(c) zComplexCoordFPrint( stdout, (c) )

/*! \brief check if a complex number is a member of an array.
 *
 * zComplexValIsIncluded() checks if a complex number \a c is
 * included in an array \a array.
 * zComplexValConjIsIncluded() checks if conjugate of \a c is
 * included in an array \a array.
 * For both functions, \a size is the size of the array, and
 * \a tol is the tolerance to regard two values are the same.
 * \return
 * zComplexValIsIncluded() returns the true value if \a c is
 * included in \a \array. Otherwise, the false value is returned.
 * zComplexValConjIsIncluded() returns the true value if conjugate
 * of \a c is included in \a \array. Otherwise, the false value
 * is returned.
 */
__ZM_EXPORT bool zComplexValIsIncluded(const zComplex *array, int size, const zComplex *c, double tol);
__ZM_EXPORT bool zComplexValConjIsIncluded(const zComplex *array, int size, const zComplex *c, double tol);

/*! \} */

__END_DECLS

#ifdef __cplusplus
inline zComplex &zComplex::create(double _re, double _im){ _zComplexCreate( this, _re, _im ); return *this; }
inline zComplex &zComplex::createPolar(double r, double angle){ _zComplexCreatePolar( this, r, angle ); return *this; }
inline zComplex &zComplex::operator=(double src){ _zComplexCreate( this, src, 0 ); return *this; }
inline zComplex &zComplex::zero(){ _zComplexZero( this ); return *this; }
inline zComplex &zComplex::touchup(double tol = zTOL){ return *zComplexTouchup( this, tol ); }
inline bool zComplex::isTol(double tol = zTOL){ return zComplexIsTol( this, tol ); }
inline bool zComplex::isReal(double tol = zTOL){ return zComplexIsReal( this, tol ); }
inline bool zComplex::isCoconj(const zComplex &anotherval){ return zComplexCoconj( this, const_cast<zComplex*>(&anotherval), zTOL ); }
inline bool zComplex::operator==(const zComplex &anotherval){ return zComplexEqual( this, const_cast<zComplex*>(&anotherval), zTOL ); }
inline bool zComplex::operator!=(const zComplex &anotherval){ return !zComplexEqual( this, const_cast<zComplex*>(&anotherval), zTOL ); }

inline double zComplex::sqrAbs(){ return _zComplexSqrAbs( this ); }
inline double zComplex::abs(){ return _zComplexAbs( this ); }
inline double zComplex::arg(){ return _zComplexArg( this ); }
inline zComplex zComplex::conj(){ zComplex conjugate; _zComplexConj( this, &conjugate ); return conjugate; }
inline zComplex zComplex::operator-(){ zComplex rev; _zComplexRev( this, &rev ); return rev; }
inline zComplex &zComplex::operator+=(const zComplex &c){ _zComplexAddDRC( this, &c ); return *this; }
inline zComplex &zComplex::operator-=(const zComplex &c){ _zComplexSubDRC( this, &c ); return *this; }
inline zComplex &zComplex::operator*=(double k){ _zComplexMulDRC( this, k ); return *this; }
inline zComplex &zComplex::operator*=(const zComplex &c){ return *zComplexCMulDRC( this, &c ); }
inline zComplex &zComplex::operator/=(double k){ zComplexDivDRC( this, k ); return *this; }
inline zComplex &zComplex::operator/=(const zComplex &c){ return *zComplexCDivDRC( this, &c ); }
inline zComplex zComplex::pow(double z){ zComplex ret; zComplexPow( this, z, &ret ); return ret; }
inline zComplex zComplex::pow(const zComplex &z){ zComplex ret; zComplexCPow( this, &z, &ret ); return ret; }
inline zComplex zComplex::pow(double z, const zComplex &ref){ zComplex ret; zComplexPowRef( this, z, &ref, &ret ); return ret; }
inline zComplex zComplex::log(double base){ zComplex ret; zComplexLog( this, base, &ret ); return ret; }
inline zComplex zComplex::log(const zComplex &base){ zComplex ret; zComplexCLog( this, &base, &ret ); return ret; }
inline zComplex zComplex::normalize(){ zComplex ret; zComplexNormalize( this, &ret ); return ret; }
inline zComplex &zComplex::fromString(const char *str){ return *zComplexFromStr( this, str ); }
inline zComplex &zComplex::fromZTK(ZTK *ztk){ return *zComplexFromZTK( this, ztk ); }

inline zComplex operator+(const zComplex &c1, const zComplex &c2){ zComplex ret; _zComplexAdd( &c1, &c2, &ret ); return ret; }
inline zComplex operator-(const zComplex &c1, const zComplex &c2){ zComplex ret; _zComplexSub( &c1, &c2, &ret ); return ret; }
inline zComplex operator*(const zComplex &c, double k){ zComplex ret; _zComplexMul( &c, k, &ret ); return ret; }
inline zComplex operator*(double k, const zComplex &c){ zComplex ret; _zComplexMul( &c, k, &ret ); return ret; }
inline zComplex operator*(const zComplex &c1, const zComplex &c2){ zComplex ret; _zComplexCMul( &c1, &c2, &ret ); return ret; }
inline zComplex operator/(const zComplex &c, double k){ zComplex ret; zComplexDiv( &c, k, &ret ); return ret; }
inline zComplex operator/(double k, const zComplex &c){ zComplex ret; zComplexInv( &c, &ret ); _zComplexMulDRC( &ret, k ); return ret; }
inline zComplex operator/(const zComplex &c1, const zComplex &c2){ zComplex ret; zComplexCDiv( &c1, &c2, &ret ); return ret; }

inline std::ostream &operator<<(std::ostream &stream, const zComplex &complex){
  if( complex.re == 0 ){
    if( complex.im == 0 )
      stream << 0;
    else{
      if( complex.im < 0 ) stream << '-';
      if( fabs( complex.im ) != 1 ) stream << complex.im;
      stream << 'i';
    }
  } else{
    stream << complex.re;
    if( complex.im != 0 ){
      stream << ( complex.im > 0 ? '+' : '-' );
      if( fabs( complex.im ) != 1 ) stream << complex.im;
      stream << 'i';
    }
  }
  return stream;
}
#endif /* __cplusplus */

#include <zm/zm_complex_pe.h> /* polynomial equation solver */

#endif /* __ZM_COMPLEX_H__ */
