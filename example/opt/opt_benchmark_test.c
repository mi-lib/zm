#include <zm/zm_opt.h>

#define TEST 13
/*
 0: Ackley
 1: Rosenbrock
 2: Beale
 3: Goldstein & Price
 4: Booth
 5: Bukin No. 6
 6: Matyas
 7: Levi No. 13
 8: three-hump camel
 9: Easom
10: Eggholder
11: McCormick
12: Styblinski-Tang
 *: mixed Gaussian
 */

#if TEST == 0
/* Ackley function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 20-20*exp(-2.0*sqrt((zSqr(x1)+zSqr(x2))/2))-exp((cos(zPIx2*x1)+cos(zPIx2*x2))/2);
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-33.0,-33.0 );
  *max = zVecCreateList( 2, 33.0, 33.0 );
  *ans = zVecCreateList( 2,  0.0,  0.0 );
  *val_ans = -zE;
}
#elif TEST == 1
/* Rosenbrock function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 100 * zSqr(x2-x1*x1) + zSqr(x1-1);
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-5.0,-5.0 );
  *max = zVecCreateList( 2, 5.0, 5.0 );
  *ans = zVecCreateList( 2, 1.0, 1.0 );
  *val_ans = 0;
}
#elif TEST == 2
/* Beale function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return zSqr(1.5-x1+x1*x2) + zSqr(2.25-x1+x1*x2*x2) + zSqr(2.625-x1+x1*x2*x2*x2);
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-4.5,-4.5 );
  *max = zVecCreateList( 2, 4.5, 4.5 );
  *ans = zVecCreateList( 2, 3.0, 0.5 );
  *val_ans = 0;
}
#elif TEST == 3
/* Goldstein & Price function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return ( 1 + zSqr(x1+x2+1) * (19-14*x1+3*x1*x1-14*x2+6*x1*x2+3*x2*x2) )
       * ( 30 + zSqr(2*x1-3*x2) * (18-32*x1+12*x1*x1+48*x2-36*x1*x2+27*x2*x2) );
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-2.0,-2.0 );
  *max = zVecCreateList( 2, 2.0, 2.0 );
  *ans = zVecCreateList( 2, 0.0,-1.0 );
  *val_ans = 3.0;
}
#elif TEST == 4
/* Booth function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return zSqr(x1+2*x2-7) + zSqr(2*x1+x2-5);
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-10.0,-10.0 );
  *max = zVecCreateList( 2, 10.0, 10.0 );
  *ans = zVecCreateList( 2, 1.0, 3.0 );
  *val_ans = 0;
}
#elif TEST == 5
/* Bukin function No. 6 */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 100 * sqrt(fabs(x2-0.01*x1*x1)) + 0.01 * fabs(x1+10);
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-15.0,-3.0 );
  *max = zVecCreateList( 2, -5.0, 3.0 );
  *ans = zVecCreateList( 2,-10.0, 1.0 );
  *val_ans = 0;
}
#elif TEST == 6
/* Matyas function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 0.26*(x1*x1+x2*x2) -0.48*x1*x2;
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-10.0,-10.0 );
  *max = zVecCreateList( 2, 10.0, 10.0 );
  *ans = zVecCreateList( 2, 0.0, 0.0 );
  *val_ans = 0;
}
#elif TEST == 7
/* Levi function No. 13 */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return zSqr(sin(zPI*3*x1)) + zSqr(x1-1)*(1+zSqr(sin(zPI*3*x2))) + zSqr(x2-1)*(1+zSqr(sin(zPIx2*x2)));
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-10.0,-10.0 );
  *max = zVecCreateList( 2, 10.0, 10.0 );
  *ans = zVecCreateList( 2, 1.0, 1.0 );
  *val_ans = 0;
}
#elif TEST == 8
/* three-hump camel function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 2*x1*x1 - 1.05*pow(x1,4) + pow(x1,6)/6 + x1*x2 + x2*x2;
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-5.0,-5.0 );
  *max = zVecCreateList( 2, 5.0, 5.0 );
  *ans = zVecCreateList( 2, 0.0, 0.0 );
  *val_ans = 0;
}
#elif TEST == 9
/* Easom function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return -cos(x1)*cos(x2)*exp(-(zSqr(x1-zPI)+zSqr(x2-zPI)));
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-100.0,-100.0 );
  *max = zVecCreateList( 2, 100.0, 100.0 );
  *ans = zVecCreateList( 2, zPI, zPI );
  *val_ans = -1;
}
#elif TEST == 10
/* eggholder function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return -(x2+47)*sin(sqrt(fabs(x2+0.5*x1+47))) - x1*sin(sqrt(fabs(x1-x2-47)));
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-512.0,-512.0 );
  *max = zVecCreateList( 2, 512.0, 512.0 );
  *ans = zVecCreateList( 2, 512.0, 404.2319 );
  *val_ans = -959.6407;
}
#elif TEST == 11
/* McCormick function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return sin(x1+x2) + zSqr(x1-x2) - 1.5*x1 + 2.5*x2;
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-1.5,-3.0 );
  *max = zVecCreateList( 2, 4.0, 4.0 );
  *ans = zVecCreateList( 2,-0.54719,-1.54719 );
  *val_ans = -2.9133;
}
#elif TEST == 12
/* Styblinski-Tang function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 0.5 * ( pow(x1,4)-16*x1*x1+5*x1 + pow(x2,4)-16*x2*x2+5*x2 );
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-5.0,-3.0 );
  *max = zVecCreateList( 2, 4.0, 4.0 );
  *ans = zVecCreateList( 2,-2.903534,-2.903534 );
  *val_ans = -39.166165 * 2;
}
#else
/* mixed Gaussian */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return   -exp(-0.3*(zSqr(x1-5)+zSqr(x2-7)))
         -5*exp(-2.0*(zSqr(x1+4)+zSqr(x2+5)))
         -2*exp(-0.5*(zSqr(x1-5)+zSqr(x2+8)))
         -3*exp(-0.3*(zSqr(x1+6)+zSqr(x2-6)));
}
void create_boundary_answer(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-10.0,-10.0 );
  *max = zVecCreateList( 2, 10.0, 10.0 );
  *ans = zVecCreateList( 2, -4.0, -5.0 );
  *val_ans = -5.0;
}
#endif


void output(char *method, zVec var, double val)
{
  printf( "[%s]\n\t", method ); zVecPrint( var );
  printf( " %.10g\n", val );
}

int main(int argc, char *argv[])
{
  zVec min, max, var, ans;
  double val, val_ans;

  create_boundary_answer( &min, &max, &ans, &val_ans );
  output( "*answer*", ans, val_ans );
  var = zVecAlloc( zVecSizeNC(min) );

  zOptSolveDIRECT( testfunc, NULL, min, max, 0, zTOL, var, &val );
  output( "DIRECT", var, val );

  zOptSolveNM( testfunc, NULL, min, max, 0, zTOL, var, &val );
  output( "Nelder-Mead", var, val );

  zOptSolveGADefault( testfunc, NULL, min, max, 0, zTOL, var, &val );
  output( "GA", var, val );

  zOptSolvePSODefault( testfunc, NULL, min, max, 0, zTOL, var, &val );
  output( "PSO", var, val );

  zVecFreeAO( 3, min, max, ans );
  return 0;
}
