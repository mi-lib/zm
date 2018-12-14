/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_tridiag - linear equation: tridiagonal equation.
 */

#ifndef __ZM_LE_TRIDIAG_H__
#define __ZM_LE_TRIDIAG_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zTridiagSolveDST, zTridiagSolve
 * - tridiag equation solver.
 *
 * 'zTridiagSolve()' solves the tridiagonal equation 'm x = d'.
 * where, 'm' = | b0 c0  0  .  .  0 |
 *              | a1 b1 c1  .  .  0 |
 *              |  0 a2 b2 c2  .  0 |
 *              |  0  0 a3  .     . |
 *              |     .        .  . |
 *              |     .           . |
 *              |  0  0  .  . an bn |
 * (Thus, the first component of 'a' and the last component of 'c'
 *  does not have any meanings.)
 * The answer will be set on 'ans'.
 * #
 * 'zTridiagSolveDST()' destructively modifies 'a',
 * 'b', 'c' and 'd' while calculating.
 * [RETURN VALUE]
 * Each of 'zTridiagSolve()' and 'zTridiagSolveDST()'
 * returns a pointer to 'ans' if succeeding to getting the
 * answer, or the null pointer if failing.
 */
__EXPORT zVec zTridiagSolveDST(zVec a, zVec b, zVec c, zVec d, zVec ans);
__EXPORT zVec zTridiagSolve(zVec a, zVec b, zVec c, zVec d, zVec ans);

__END_DECLS

#endif /* __ZM_LE_TRIDIAG_H__ */
