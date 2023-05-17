/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_tridiag - linear equation: tridiagonal equation.
 */

#ifndef __ZM_LE_TRIDIAG_H__
#define __ZM_LE_TRIDIAG_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief tridiag equation solver.
 *
 * zTridiagSolve() solves the tridiagonal equation 'm x = d', where,
 * 'm' = | b0 c0  0  .  .  0 |
 *       | a1 b1 c1  .  .  0 |
 *       |  0 a2 b2 c2  .  0 |
 *       |  0  0 a3  .     . |
 *       |     .        .  . |
 *       |     .           . |
 *       |  0  0  .  . an bn |
 * (Thus, the first component of \a a and the last component of \a c
 *  does not have any meanings.)
 * The answer will be set on \a ans.
 *
 * zTridiagSolveDST() destructively modifies \a a, \a b, \a c and \a d
 * while the calculation.
 * \return
 * zTridiagSolve() and zTridiagSolveDST() return a pointer to \a ans
 * in successful cases, or the null pointer otherwise.
 */
__ZM_EXPORT zVec zTridiagSolveDST(zVec a, zVec b, zVec c, zVec d, zVec ans);
__ZM_EXPORT zVec zTridiagSolve(zVec a, zVec b, zVec c, zVec d, zVec ans);

__END_DECLS

#endif /* __ZM_LE_TRIDIAG_H__ */
