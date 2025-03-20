/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_tridiag - linear equation: tridiagonal equation.
 */

#ifndef __ZM_LE_TRIDIAG_H__
#define __ZM_LE_TRIDIAG_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief tridiagonal equation solver.
 *
 * zLETridiagSolve() solves a linear tridiagonal equation m x = d, where
 *  m = | b0 c0  0  .    .      0   |
 *      | a1 b1 c1  .    .      0   |
 *      |  0 a2 b2 c2    .      0   |
 *      |  0  0 a3  .    .      .   |
 *      |     .          .      .   |
 *      |     .          .      .   |
 *      |  0  0  .  . a(n-1) b(n-1) |
 *  d = [ d0 d1 ... d(n-1) ]^T
 * Vectors \a a, \a b, \a c, and \a d store the above (a0, ..., a(n-1)), (b0, ..., b(n-1)), (c0, ..., c(n-1)),
 * and (d0, ..., d(n-1)), respectively. Note that the first component of \a a and the last component of
 * \a c does not have any meanings.
 * The answer will be set in \a ans.
 *
 * zLETridiagSolveDST() destructively modifies \a a, \a b, \a c and \a d during solving the equation.
 * \return
 * zLETridiagSolve() and zLETridiagSolveDST() return a pointer \a ans if they succeed. Otherwise, they
 * return the null pointer.
 */
__ZM_EXPORT zVec zLETridiagSolveDST(zVec a, zVec b, zVec c, zVec d, zVec ans);
__ZM_EXPORT zVec zLETridiagSolve(const zVec a, const zVec b, const zVec c, const zVec d, zVec ans);

__END_DECLS

#endif /* __ZM_LE_TRIDIAG_H__ */
