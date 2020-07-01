#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <zm/zm_opt.h>
#include <time.h>

typedef struct{
  char *name;
  double (* testfunc)(zVec,void*);
  void (* create_boundary_answer)(zVec*,zVec*,zVec*,double*);
} sample_t;

#define TESTFUNC( name )         testfunc_##name
#define BOUNDARY_ANSWER( name )  create_boundary_answer_##name
#define SAMPLE( name )           sample_##name
#define SAMPLE_DEF( name ) sample_t SAMPLE( name ) = { #name, TESTFUNC( name ), BOUNDARY_ANSWER( name ) }

/* Ackley function */
double TESTFUNC( ackley )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 20-20*exp(-2.0*sqrt((zSqr(x1)+zSqr(x2))/2))-exp((cos(zPIx2*x1)+cos(zPIx2*x2))/2);
}
void BOUNDARY_ANSWER( ackley )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-33.0,-33.0 );
  *max = zVecCreateList( 2, 33.0, 33.0 );
  *ans = zVecCreateList( 2,  0.0,  0.0 );
  *val_ans = -zE;
}
SAMPLE_DEF( ackley );

/* Rosenbrock function */
double TESTFUNC( rosenbrock )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 100 * zSqr(x2-x1*x1) + zSqr(x1-1);
}
void BOUNDARY_ANSWER( rosenbrock )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-5.0,-5.0 );
  *max = zVecCreateList( 2, 5.0, 5.0 );
  *ans = zVecCreateList( 2, 1.0, 1.0 );
  *val_ans = 0;
}
SAMPLE_DEF( rosenbrock );

/* Beale function */
double TESTFUNC( beale )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return zSqr(1.5-x1+x1*x2) + zSqr(2.25-x1+x1*x2*x2) + zSqr(2.625-x1+x1*x2*x2*x2);
}
void BOUNDARY_ANSWER( beale )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-4.5,-4.5 );
  *max = zVecCreateList( 2, 4.5, 4.5 );
  *ans = zVecCreateList( 2, 3.0, 0.5 );
  *val_ans = 0;
}
SAMPLE_DEF( beale );

/* Goldstein & Price function */
double TESTFUNC( gp )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return ( 1 + zSqr(x1+x2+1) * (19-14*x1+3*x1*x1-14*x2+6*x1*x2+3*x2*x2) )
       * ( 30 + zSqr(2*x1-3*x2) * (18-32*x1+12*x1*x1+48*x2-36*x1*x2+27*x2*x2) );
}
void BOUNDARY_ANSWER( gp )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-2.0,-2.0 );
  *max = zVecCreateList( 2, 2.0, 2.0 );
  *ans = zVecCreateList( 2, 0.0,-1.0 );
  *val_ans = 3.0;
}
SAMPLE_DEF( gp );

/* Booth function */
double TESTFUNC( booth )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return zSqr(x1+2*x2-7) + zSqr(2*x1+x2-5);
}
void BOUNDARY_ANSWER( booth )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-10.0,-10.0 );
  *max = zVecCreateList( 2, 10.0, 10.0 );
  *ans = zVecCreateList( 2, 1.0, 3.0 );
  *val_ans = 0;
}
SAMPLE_DEF( booth );

/* Bukin function No. 6 */
double TESTFUNC( bukin6 )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 100 * sqrt(fabs(x2-0.01*x1*x1)) + 0.01 * fabs(x1+10);
}
void BOUNDARY_ANSWER( bukin6 )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-15.0,-3.0 );
  *max = zVecCreateList( 2, -5.0, 3.0 );
  *ans = zVecCreateList( 2,-10.0, 1.0 );
  *val_ans = 0;
}
SAMPLE_DEF( bukin6 );

/* Matyas function */
double TESTFUNC( matyas )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 0.26*(x1*x1+x2*x2) -0.48*x1*x2;
}
void BOUNDARY_ANSWER( matyas )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-10.0,-10.0 );
  *max = zVecCreateList( 2, 10.0, 10.0 );
  *ans = zVecCreateList( 2, 0.0, 0.0 );
  *val_ans = 0;
}
SAMPLE_DEF( matyas );

/* Levi function No. 13 */
double TESTFUNC( levi13 )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return zSqr(sin(zPI*3*x1)) + zSqr(x1-1)*(1+zSqr(sin(zPI*3*x2))) + zSqr(x2-1)*(1+zSqr(sin(zPIx2*x2)));
}
void BOUNDARY_ANSWER( levi13 )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-10.0,-10.0 );
  *max = zVecCreateList( 2, 10.0, 10.0 );
  *ans = zVecCreateList( 2, 1.0, 1.0 );
  *val_ans = 0;
}
SAMPLE_DEF( levi13 );

/* three-hump camel function */
double TESTFUNC( thcamel )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 2*x1*x1 - 1.05*pow(x1,4) + pow(x1,6)/6 + x1*x2 + x2*x2;
}
void BOUNDARY_ANSWER( thcamel )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-5.0,-5.0 );
  *max = zVecCreateList( 2, 5.0, 5.0 );
  *ans = zVecCreateList( 2, 0.0, 0.0 );
  *val_ans = 0;
}
SAMPLE_DEF( thcamel );

/* Easom function */
double TESTFUNC( easom )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return -cos(x1)*cos(x2)*exp(-(zSqr(x1-zPI)+zSqr(x2-zPI)));
}
void BOUNDARY_ANSWER( easom )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-100.0,-100.0 );
  *max = zVecCreateList( 2, 100.0, 100.0 );
  *ans = zVecCreateList( 2, zPI, zPI );
  *val_ans = -1;
}
SAMPLE_DEF( easom );

/* eggholder function */
double TESTFUNC( eggholder )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return -(x2+47)*sin(sqrt(fabs(x2+0.5*x1+47))) - x1*sin(sqrt(fabs(x1-x2-47)));
}
void BOUNDARY_ANSWER( eggholder )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-512.0,-512.0 );
  *max = zVecCreateList( 2, 512.0, 512.0 );
  *ans = zVecCreateList( 2, 512.0, 404.2319 );
  *val_ans = -959.6407;
}
SAMPLE_DEF( eggholder );

/* McCormick function */
double TESTFUNC( mccormick )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return sin(x1+x2) + zSqr(x1-x2) - 1.5*x1 + 2.5*x2;
}
void BOUNDARY_ANSWER( mccormick )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-1.5,-3.0 );
  *max = zVecCreateList( 2, 4.0, 4.0 );
  *ans = zVecCreateList( 2,-0.54719,-1.54719 );
  *val_ans = -2.9133;
}
SAMPLE_DEF( mccormick );

/* Styblinski-Tang function */
double TESTFUNC( st )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return 0.5 * ( pow(x1,4)-16*x1*x1+5*x1 + pow(x2,4)-16*x2*x2+5*x2 );
}
void BOUNDARY_ANSWER( st )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-5.0,-3.0 );
  *max = zVecCreateList( 2, 4.0, 4.0 );
  *ans = zVecCreateList( 2,-2.903534,-2.903534 );
  *val_ans = -39.166165 * 2;
}
SAMPLE_DEF( st );

/* Gaussian mixture model */
double TESTFUNC( gmm )(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return   -exp(-0.3*(zSqr(x1-5)+zSqr(x2-7)))
         -5*exp(-2.0*(zSqr(x1+4)+zSqr(x2+5)))
         -2*exp(-0.5*(zSqr(x1-5)+zSqr(x2+8)))
         -3*exp(-0.3*(zSqr(x1+6)+zSqr(x2-6)));
}
void BOUNDARY_ANSWER( gmm )(zVec *min, zVec *max, zVec *ans, double *val_ans)
{
  *min = zVecCreateList( 2,-10.0,-10.0 );
  *max = zVecCreateList( 2, 10.0, 10.0 );
  *ans = zVecCreateList( 2, -4.0, -5.0 );
  *val_ans = -5.0;
}
SAMPLE_DEF( gmm );

sample_t *sample[] = {
  &SAMPLE( ackley ),     /* Ackley */
  &SAMPLE( rosenbrock ), /* Rosenbrock */
  &SAMPLE( beale ),      /* Beale */
  &SAMPLE( gp ),         /* Goldstein & Price */
  &SAMPLE( booth ),      /* Booth */
  &SAMPLE( bukin6 ),     /* Bukin No. 6 */
  &SAMPLE( matyas ),     /* Matyas */
  &SAMPLE( levi13 ),     /* Levi No. 13 */
  &SAMPLE( thcamel ),    /* three-hump camel */
  &SAMPLE( easom ),      /* Easom */
  &SAMPLE( eggholder ),  /* Eggholder */
  &SAMPLE( mccormick ),  /* McCormick */
  &SAMPLE( st ),         /* Styblinski-Tang */
  &SAMPLE( gmm ),        /* mixed Gaussian */
  NULL,
};

long timespect2nanosec(struct timespec *tp1, struct timespec *tp2)
{
  return ( tp2->tv_sec - tp1->tv_sec ) * 1000000000 + ( tp2->tv_nsec - tp2->tv_nsec );
}

void output(char *method, struct timespec *tp1, struct timespec *tp2, zVec var, double val)
{
  printf( " (%s)", method );
  strlen( method ) > 5 ? printf( "\t" ) : printf( "\t\t" );
  printf( "%ld %.10g\t", timespect2nanosec(tp1,tp2), val );
  zVecPrint( var );
}

int main(int argc, char *argv[])
{
  zVec min, max, var, ans;
  double val, val_ans;
  sample_t **sp;
  struct timespec tp1, tp2;

  for( sp=sample; *sp; sp++ ){
    printf( "[%s]\n", (*sp)->name );

    (*sp)->create_boundary_answer( &min, &max, &ans, &val_ans );
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tp1 );
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tp2 );
    output( "*answer*", &tp1, &tp2, ans, val_ans );
    var = zVecAlloc( zVecSizeNC(min) );

    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tp1 );
    zOptSolveDIRECT( (*sp)->testfunc, NULL, min, max, 0, zTOL, var, &val );
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tp2 );
    output( "DIRECT", &tp1, &tp2, var, val );

    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tp1 );
    zOptSolveNM( (*sp)->testfunc, NULL, min, max, 0, zTOL, var, &val );
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tp2 );
    output( "Nelder-Mead", &tp1, &tp2, var, val );

    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tp1 );
    zOptSolveGADefault( (*sp)->testfunc, NULL, min, max, 0, zTOL, var, &val );
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tp2 );
    output( "GA", &tp1, &tp2, var, val );

    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tp1 );
    zOptSolvePSODefault( (*sp)->testfunc, NULL, min, max, 0, zTOL, var, &val );
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &tp2 );
    output( "PSO", &tp1, &tp2, var, val );

    zVecFreeAO( 3, min, max, ans );
  }
  return 0;
}
