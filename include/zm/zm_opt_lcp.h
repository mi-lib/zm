/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lcp - optimization tools: linear complementarity problem.
 */

#ifndef __ZM_OPT_LCP_H__
#define __ZM_OPT_LCP_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zLCPSolveLemke, zLCPSolveIP
 * - linear complementarity problem solver.
 *
 * 'zLCPSolveLemke()' solves a linear complementarity problem
 * (LCP) by Lemke s method (1965 C. E. Lemke) which finds a
 * combination of two vectors 'w' and 'z' which satisfies:
 *   'w - m z = q', 'w >= 0' and 'z >= 0'.
 * where 'm' is a square matrix.
 * #
 * 'zLCPSolveIP()' also solves LCP by Potra s predictor-
 * corrector algorithm for infeasible-interior-point method
 * (1994 F. Potra and R Sheng).
 * [RETURN VALUE]
 * 'zLCPSolveLemke()' and 'zLCPSolveIP()' return the true
 * value if succeeding to get the solution, or the false
 * value if there are vector/matrix size mismatch, bad memory
 * allocation, non-feasible solution and infinite solution
 * are found.
 * [NOTES]
 * Since 'w' is often not required to be answered, it is
 * allowed to give the null pointer for 'w' to skip the
 * operation.
 */
__EXPORT bool zLCPSolveLemke(zMat m, zVec q, zVec w, zVec z);
__EXPORT bool zLCPSolveIP(zMat m, zVec q, zVec w, zVec z);

__END_DECLS

#endif /* __ZM_OPT_LCP_H__ */
