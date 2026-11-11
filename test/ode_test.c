#include <zm/zm_ode.h>

/* ODE quadrature methods to be tested */
typedef struct{
  const char *name;
  zODE *(* create)(zODE*,int,int,zVec (*)(double,zVec,void*,zVec));
  void (* destroy)(zODE*);
  zVec (* update)(zODE*,double,zVec,double,void*);
  void (* init_history)(zODE*,double,zVec,void*);
  zVec (* create_history)(zODE*,double*,zVec,double,void*);
  int step;  /* step number of history (only for multistep methods) */
  int order; /* theoretical order of accuracy */
} method_t;

#define METHOD(type,step,order) { #type, zODECreate##type, zODEDestroy##type, zODEUpdate##type, zODEInitHistory##type, zODECreateHistory##type, step, order }

static const method_t method_onestep[] = {
  METHOD( Euler,  0, 1 ),
  METHOD( Heun,   0, 3 ), /* named after the third-order Runge-Kutta method by Heun */
  METHOD( RK4,    0, 4 ),
  METHOD( RKG,    0, 4 ),
  METHOD( RKF45,  0, 5 ), /* the fifth-order solution is adopted */
  METHOD( CK45,   0, 5 ),
  METHOD( DP45,   0, 5 ),
  METHOD( BEuler, 0, 1 ),
  METHOD( Gauss,  0, 4 ),
  METHOD( Radau,  0, 3 ),
  METHOD( TR,     0, 2 ),
};
#define METHOD_ONESTEP_NUM ( sizeof(method_onestep)/sizeof(method_t) )

static void method_assign(zODE *ode, const method_t *method)
{
  ode->create = method->create;
  ode->destroy = method->destroy;
  ode->update = method->update;
  ode->init_history = method->init_history;
  ode->create_history = method->create_history;
  ode->order = method->order;
  zODEAssignFunc( ode, NULL, NULL );
}

static zODE *method_create(zODE *ode, const method_t *method, int dim, int step, zVec (* f)(double,zVec,void*,zVec))
{
  method_assign( ode, method );
  return zODECreate( ode, dim, step, f );
}

/* x' = -y, y' = x (unit circle) */
zVec circle(double t, zVec x, void *dummy, zVec dx)
{
  zVecSetElemNC( dx, 0,-zVecElemNC(x,1) );
  zVecSetElemNC( dx, 1, zVecElemNC(x,0) );
  return dx;
}

/* error of the solution of the unit circle from the exact one at the time t */
double circle_error(zVec x, double t, double r)
{
  return sqrt( zSqr( zVecElemNC(x,0) - r*cos(t) ) + zSqr( zVecElemNC(x,1) - r*sin(t) ) );
}

/* integrate the unit circle from t=0 to t=T with a one-step method, and return the error at T. */
double circle_integrate(const method_t *method, double dt, double T)
{
  zODE ode;
  zVec x;
  int i, n;
  double err;

  n = (int)( T / dt + 0.5 );
  method_create( &ode, method, 2, method->step, circle );
  x = zVecCreateList( 2, 1.0, 0.0 );
  for( i=0; i<n; i++ )
    zODEUpdate( &ode, i*dt, x, dt, NULL );
  err = circle_error( x, n*dt, 1.0 );
  zVecFree( x );
  zODEDestroy( &ode );
  return err;
}

/* observed order of accuracy from errors with time steps dt and dt/2 */
double observed_order(double err1, double err2)
{
  return log( err1 / err2 ) / log( 2 );
}

#define ORDER_TOL 0.2

/* 1. order of accuracy of one-step methods */
void assert_ode_order_onestep(void)
{
  const double dt = 0.05, T = 2.0;
  int i;
  double order;
  bool result = true;

  for( i=0; i<METHOD_ONESTEP_NUM; i++ ){
    order = observed_order( circle_integrate( &method_onestep[i], dt, T ), circle_integrate( &method_onestep[i], 0.5*dt, T ) );
    if( fabs( order - method_onestep[i].order ) > ORDER_TOL ){
      eprintf( "%s: observed order = %g (expected %d)\n", method_onestep[i].name, order, method_onestep[i].order );
      result = false;
    }
  }
  zAssert( zODEUpdate (order of accuracy of one-step methods), result );
}

/* x' = t^6 for t>=0 and x' = 0 for t<0, whose solution is x = t^7/7 for t>=0. Since x is constant and
 * x' is zero for t<0, the initial histories of Adams's and Gear's methods are consistent with it. */
zVec poly6(double t, zVec x, void *dummy, zVec dx)
{
  zVecSetElemNC( dx, 0, t > 0 ? pow( t, 6 ) : 0 );
  return dx;
}

/* integrate poly6 from t=0 to t=1 with a multistep method, and return the error at t=1. */
double poly6_integrate(const method_t *method, int step, double dt)
{
  zODE ode;
  zVec x;
  int i, n;
  double err;

  n = (int)( 1.0 / dt + 0.5 );
  method_create( &ode, method, 1, step, poly6 );
  x = zVecCreateList( 1, 0.0 );
  zODEInitHistory( &ode, 0, x, NULL );
  for( i=0; i<n; i++ )
    zODEUpdate( &ode, i*dt, x, dt, NULL );
  err = fabs( zVecElemNC(x,0) - 1.0/7 );
  zVecFree( x );
  zODEDestroy( &ode );
  return err;
}

/* 2. order of accuracy of multistep methods */
void assert_ode_order_multistep(void)
{
  const method_t adams = METHOD( Adams, 0, 0 ), gear = METHOD( Gear, 0, 0 );
  const double dt = 0.00625;
  int step;
  double order;
  bool result = true;

  for( step=1; step<=6; step++ ){ /* the order of Adams's method with k steps is k. */
    order = observed_order( poly6_integrate( &adams, step, dt ), poly6_integrate( &adams, step, 0.5*dt ) );
    if( fabs( order - step ) > ORDER_TOL ){
      eprintf( "Adams (%d steps): observed order = %g (expected %d)\n", step, order, step );
      result = false;
    }
  }
  /* the order of Gear's method with k steps is k. Since the nonlinear equation is solved with the
   * tolerance zTOL, the error saturates for more steps. */
  for( step=1; step<=4; step++ ){
    order = observed_order( poly6_integrate( &gear, step, dt ), poly6_integrate( &gear, step, 0.5*dt ) );
    if( fabs( order - step ) > ORDER_TOL ){
      eprintf( "Gear (%d steps): observed order = %g (expected %d)\n", step, order, step );
      result = false;
    }
  }
  zAssert( zODEUpdate (order of accuracy of multistep methods), result );
}

/* integrate the unit circle with a large time step, and check the transition of the norm of the state.
 * sgn = 1: the norm monotonically increases, sgn = -1: monotonically decreases, sgn = 0: is preserved. */
bool circle_norm_check(const method_t *method, int sgn)
{
  const double dt = 0.5, T = 100.0, tol = 1.0e-9;
  zODE ode;
  zVec x;
  int i, n;
  double norm, norm_prev = 1.0;
  bool result = true;

  n = (int)( T / dt + 0.5 );
  method_create( &ode, method, 2, method->step, circle );
  x = zVecCreateList( 2, 1.0, 0.0 );
  for( i=0; i<n; i++ ){
    zODEUpdate( &ode, i*dt, x, dt, NULL );
    norm = zVecNorm( x );
    if( ( sgn == 0 && fabs( norm - 1 ) > tol ) ||
        ( sgn > 0 && norm <= norm_prev ) || ( sgn < 0 && norm >= norm_prev ) ){
      eprintf( "%s: norm = %.15g -> %.15g at step %d\n", method->name, norm_prev, norm, i );
      result = false;
      break;
    }
    norm_prev = norm;
  }
  zVecFree( x );
  zODEDestroy( &ode );
  return result;
}

/* 3. qualitative behaviors for a harmonic oscillation */
void assert_ode_circle_norm(void)
{
  const method_t euler = METHOD( Euler, 0, 1 ), beuler = METHOD( BEuler, 0, 1 ), radau = METHOD( Radau, 0, 3 );
  const method_t tr = METHOD( TR, 0, 2 ), gauss = METHOD( Gauss, 0, 4 );
  bool result = true;

  if( !circle_norm_check( &euler, 1 ) ) result = false;  /* explicit Euler method amplifies */
  if( !circle_norm_check( &beuler, -1 ) ) result = false; /* backward Euler method damps */
  if( !circle_norm_check( &radau, -1 ) ) result = false;  /* Radau method damps */
  if( !circle_norm_check( &tr, 0 ) ) result = false;     /* trapezoidal method preserves norm */
  if( !circle_norm_check( &gauss, 0 ) ) result = false;  /* Gauss method preserves norm */
  zAssert( zODEUpdate (norm of harmonic oscillation), result );
}

#define STIFF_A 100.0

/* x' = -A x */
zVec stiff(double t, zVec x, void *dummy, zVec dx)
{
  zVecSetElemNC( dx, 0, -STIFF_A * zVecElemNC(x,0) );
  return dx;
}

/* stability functions R(z) of A-stable one-step methods, with which x <- R(z) x for x' = -A x and z = -A dt */
double stiff_r_beuler(double z){ return 1 / ( 1 - z ); }
double stiff_r_tr(double z){ return ( 1 + z/2 ) / ( 1 - z/2 ); }
double stiff_r_gauss(double z){ return ( 1 + z/2 + z*z/12 ) / ( 1 - z/2 + z*z/12 ); }
double stiff_r_radau(double z){ return ( 1 + z/3 ) / ( 1 - 2*z/3 + z*z/6 ); }

/* 4. stability for a stiff equation */
void assert_ode_stiff(void)
{
  const method_t method[] = { /* A-stable one-step methods */
    METHOD( BEuler, 0, 1 ),
    METHOD( TR,     0, 2 ),
    METHOD( Gauss,  0, 4 ),
    METHOD( Radau,  0, 3 ),
  };
  double (* r[])(double) = { stiff_r_beuler, stiff_r_tr, stiff_r_gauss, stiff_r_radau };
  const method_t gear = METHOD( Gear, 3, 3 );
  const double dt = 0.1; /* A dt = 10, for which explicit methods diverge */
  const double tol = 1.0e-8;
  const int n = 20;
  zODE ode;
  zVec x;
  int i, k;
  double xabs, xabs_prev, x_exact;
  bool result = true;

  /* one-step methods: |x| monotonically decreases, and x coincides with R(z)^i, even after it becomes
   * much smaller than zTOL, since the nonlinear equation is solved with a tolerance relative to x. */
  for( k=0; k<sizeof(method)/sizeof(method_t); k++ ){
    method_create( &ode, &method[k], 1, method[k].step, stiff );
    x = zVecCreateList( 1, 1.0 );
    for( x_exact=xabs_prev=1.0, i=0; i<5*n; i++ ){
      zODEUpdate( &ode, i*dt, x, dt, NULL );
      x_exact *= r[k]( -STIFF_A * dt );
      if( ( xabs = fabs( zVecElemNC(x,0) ) ) >= xabs_prev || fabs( zVecElemNC(x,0) - x_exact ) > tol * fabs( x_exact ) ){
        eprintf( "%s: |x| = %g -> %g at step %d (exact: %g)\n", method[k].name, xabs_prev, xabs, i, fabs( x_exact ) );
        result = false;
        break;
      }
      xabs_prev = xabs;
    }
    zVecFree( x );
    zODEDestroy( &ode );
  }
  /* multistep method: |x| does not monotonically decrease, but is bounded and decays. */
  method_create( &ode, &gear, 1, gear.step, stiff );
  x = zVecCreateList( 1, 1.0 );
  zODEInitHistory( &ode, 0, x, NULL );
  for( i=0; i<n; i++ ){
    zODEUpdate( &ode, i*dt, x, dt, NULL );
    if( ( xabs = fabs( zVecElemNC(x,0) ) ) > 1.0 ){
      eprintf( "%s: |x| = %g at step %d\n", gear.name, xabs, i );
      result = false;
      break;
    }
  }
  if( fabs( zVecElemNC(x,0) ) > 1.0e-3 ){
    eprintf( "%s: |x| = %g after %d steps\n", gear.name, fabs( zVecElemNC(x,0) ), n );
    result = false;
  }
  zVecFree( x );
  zODEDestroy( &ode );
  zAssert( zODEUpdate (stability for a stiff equation), result );
}

#define FORCED_C 1.0e-13

/* x' = c - x */
zVec forced(double t, zVec x, void *c, zVec dx)
{
  zVecSetElemNC( dx, 0, *(double *)c - zVecElemNC(x,0) );
  return dx;
}

/* 5. implicit methods starting from the zero state */
void assert_ode_implicit_zero(void)
{
  const method_t method[] = {
    METHOD( BEuler, 0, 1 ),
    METHOD( TR,     0, 2 ),
    METHOD( Gauss,  0, 4 ),
    METHOD( Radau,  0, 3 ),
    METHOD( Gear,   1, 1 ), /* Gear's method with one step is identical with backward Euler method */
  };
  double (* r[])(double) = { stiff_r_beuler, stiff_r_tr, stiff_r_gauss, stiff_r_radau, stiff_r_beuler };
  const double cs[] = { 0, 1.0, FORCED_C }, dt = 0.1, tol = 1.0e-8;
  const int n = 10;
  zODE ode;
  zVec x;
  int i, j, k;
  double c, x_exact;
  bool result = true;

  /* since y = x - c satisfies y' = -y, the solution is x = c ( 1 - R(-dt)^i ). The tolerance of the residual
   * of the nonlinear equation is zero at the zero state, which has to be updated properly even if c dt is
   * much smaller than zTOL, and has to remain zero if c = 0. */
  for( j=0; j<sizeof(cs)/sizeof(double); j++ ){
    c = cs[j];
    for( k=0; k<sizeof(method)/sizeof(method_t); k++ ){
      method_create( &ode, &method[k], 1, method[k].step, forced );
      x = zVecCreateList( 1, 0.0 );
      zODEInitHistory( &ode, 0, x, &c );
      for( x_exact=1, i=0; i<n; i++ ){
        zODEUpdate( &ode, i*dt, x, dt, &c );
        x_exact *= r[k]( -dt );
        if( fabs( zVecElemNC(x,0) - c * ( 1 - x_exact ) ) > tol * fabs( c * ( 1 - x_exact ) ) ){
          eprintf( "%s: x = %.17g at step %d (exact: %.17g)\n", method[k].name, zVecElemNC(x,0), i, c * ( 1 - x_exact ) );
          result = false;
          break;
        }
      }
      zVecFree( x );
      zODEDestroy( &ode );
    }
  }
  zAssert( zODEUpdate (implicit methods from the zero state), result );
}

/* 6. rejection of the zero time step */
void assert_ode_zero_dt(void)
{
  zODE ode;
  zVec x;
  int k;
  bool result = true;

  for( k=0; k<METHOD_ONESTEP_NUM; k++ ){
    x = zVecCreateList( 2, 1.0, 0.0 );
    method_create( &ode, &method_onestep[k], 2, 0, circle );
    if( zODEUpdate( &ode, 0, x, 0, NULL ) ) result = false;
    zODEDestroy( &ode );
    method_assign( &ode, &method_onestep[k] );
    zODECreateDC( &ode, 2, 0, circle );
    if( zODEUpdateDC( &ode, 0, x, 0, NULL ) ) result = false;
    zODEDestroyDC( &ode );
    if( zVecElemNC(x,0) != 1.0 || zVecElemNC(x,1) != 0.0 ){
      eprintf( "%s: x = (%g, %g) after a zero time step\n", method_onestep[k].name, zVecElemNC(x,0), zVecElemNC(x,1) );
      result = false;
    }
    zVecFree( x );
  }
  zAssert( zODEUpdate (rejection of the zero time step), result );
}

/* x'' = -x */
zVec oscillator(double t, zVec x, zVec dx, void *dummy, zVec ddx)
{
  zVecSetElemNC( ddx, 0, -zVecElemNC(x,0) );
  return ddx;
}

/* energy error of a harmonic oscillation x'' = -x with x(0) = 1 and x'(0) = 0 */
double oscillator_energy_error(zVec x, zVec dx)
{
  return 0.5 * ( zSqr( zVecElemNC(x,0) ) + zSqr( zVecElemNC(dx,0) ) ) - 0.5;
}

/* integrate a harmonic oscillation for a long period, and return the maximum absolute energy errors in
 * the first and the second halves of the period. If drift is not null, it is set for true if the energy
 * monotonically decreases. */
void oscillator_integrate(zODE2 *ode, double *err1, double *err2, bool *drift)
{
  const double dt = 0.1;
  const int n = 10000;
  zVec x, dx;
  int i;
  double e, e_prev = 0;

  zODE2Create( ode, 1, 0, oscillator );
  x = zVecCreateList( 1, 1.0 );
  dx = zVecCreateList( 1, 0.0 );
  zODE2InitHistory( ode, 0, x, dx, NULL );
  *err1 = *err2 = 0;
  if( drift ) *drift = true;
  for( i=0; i<n; i++ ){
    zODE2Update( ode, i*dt, x, dx, dt, NULL );
    e = oscillator_energy_error( x, dx );
    if( i < n/2 ){
      if( fabs( e ) > *err1 ) *err1 = fabs( e );
    } else{
      if( fabs( e ) > *err2 ) *err2 = fabs( e );
    }
    if( drift && e >= e_prev ) *drift = false;
    e_prev = e;
  }
  zVecFree( x );
  zVecFree( dx );
  zODE2Destroy( ode );
}

/* 7. symplecticity of second-order ODE solvers */
void assert_ode2_symplectic(void)
{
  zODE2 ode;
  double err1, err2;
  bool drift, result = true;

  /* symplectic methods: the energy error is bounded. */
  zODE2Assign( &ode, Leapfrog, NULL, NULL, NULL, NULL );
  oscillator_integrate( &ode, &err1, &err2, NULL );
  if( err2 > 1.1 * err1 ){
    eprintf( "Leapfrog: energy error %g in the first half, %g in the second half\n", err1, err2 );
    result = false;
  }
  zODE2Assign( &ode, Symplectic, NULL, NULL, NULL, NULL );
  oscillator_integrate( &ode, &err1, &err2, NULL );
  if( err2 > 1.1 * err1 ){
    eprintf( "Symplectic: energy error %g in the first half, %g in the second half\n", err1, err2 );
    result = false;
  }
  /* a non-symplectic method: the energy monotonically drifts. */
  zODE2Assign( &ode, Regular, NULL, NULL, NULL, NULL );
  zODE2AssignRegular( &ode, RK4 );
  oscillator_integrate( &ode, &err1, &err2, &drift );
  if( !drift || err2 < 1.5 * err1 ){
    eprintf( "RK4: energy error %g in the first half, %g in the second half (%s)\n", err1, err2, drift ? "monotonic" : "not monotonic" );
    result = false;
  }
  zAssert( zODE2Update (symplecticity), result );
}

/* rejection of the zero time step by a second-order ODE solver */
bool ode2_zero_dt_check(zODE2 *ode, const char *name)
{
  zVec x, v;
  bool result = true;

  zODE2Create( ode, 1, 0, oscillator );
  x = zVecCreateList( 1, 1.0 );
  v = zVecCreateList( 1, 0.0 );
  if( zODE2Update( ode, 0, x, v, 0, NULL ) || zVecElemNC(x,0) != 1.0 || zVecElemNC(v,0) != 0.0 ){
    eprintf( "%s: x = %g, v = %g after a zero time step\n", name, zVecElemNC(x,0), zVecElemNC(v,0) );
    result = false;
  }
  if( zODE2Update( ode, 0, x, v, 0.1, NULL ) != x ){
    eprintf( "%s: zODE2Update does not return x\n", name );
    result = false;
  }
  zVecFree( x );
  zVecFree( v );
  zODE2Destroy( ode );
  return result;
}

/* 8. rejection of the zero time step by second-order ODE solvers */
void assert_ode2_zero_dt(void)
{
  zODE2 ode;
  bool result = true;

  zODE2Assign( &ode, Regular, NULL, NULL, NULL, NULL );
  zODE2AssignRegular( &ode, RK4 );
  if( !ode2_zero_dt_check( &ode, "Regular" ) ) result = false;
  zODE2Assign( &ode, Symplectic, NULL, NULL, NULL, NULL );
  if( !ode2_zero_dt_check( &ode, "Symplectic" ) ) result = false;
  zODE2Assign( &ode, Leapfrog, NULL, NULL, NULL, NULL );
  if( !ode2_zero_dt_check( &ode, "Leapfrog" ) ) result = false;
  zAssert( zODE2Update (rejection of the zero time step), result );
}

/* 9. accuracy of deferred correction */
void assert_ode_dc_accuracy(void)
{
  const double dt = 0.5, T = 10.0;
  zODE ode;
  zVec x;
  int i, k, n;
  double err, err_plain, tol;
  bool is_adaptive, result = true;

  n = (int)( T / dt + 0.5 );
  /* the error in each time step is controlled to be in the order of ZM_ODE_TOL. */
  tol = 10 * n * ZM_ODE_TOL;
  /* multistep methods are not available for deferred correction, since their histories are
   * disturbed by trial updates. */
  for( k=0; k<METHOD_ONESTEP_NUM; k++ ){
    err_plain = circle_integrate( &method_onestep[k], dt, T );
    /* embedded Runge-Kutta methods control the step size by themselves, so that the deferred correction
     * does not necessarily reduce the error. */
    is_adaptive = method_onestep[k].create == zODECreateRKF45 ||
                  method_onestep[k].create == zODECreateCK45 ||
                  method_onestep[k].create == zODECreateDP45;
    method_assign( &ode, &method_onestep[k] );
    zODECreateDC( &ode, 2, 0, circle );
    x = zVecCreateList( 2, 1.0, 0.0 );
    for( i=0; i<n; i++ )
      zODEUpdateDC( &ode, i*dt, x, dt, NULL );
    if( ( err = circle_error( x, n*dt, 1.0 ) ) > tol || ( !is_adaptive && err > err_plain ) ){
      eprintf( "%s: error = %g with deferred correction (tolerance %g, error without it %g)\n", method_onestep[k].name, err, tol, err_plain );
      result = false;
    }
    zVecFree( x );
    zODEDestroyDC( &ode );
  }
  zAssert( zODEUpdateDC (accuracy), result );
}

/* integrate the unit circle from t=0 to t=T with deferred correction, and return the error at T. */
double circle_integrate_dc(const method_t *method, double dt, double T)
{
  zODE ode;
  zVec x;
  int i, n;
  double err;

  n = (int)( T / dt + 0.5 );
  method_assign( &ode, method );
  zODECreateDC( &ode, 2, 0, circle );
  x = zVecCreateList( 2, 1.0, 0.0 );
  for( i=0; i<n; i++ )
    zODEUpdateDC( &ode, i*dt, x, dt, NULL );
  err = circle_error( x, n*dt, 1.0 );
  zVecFree( x );
  zODEDestroyDC( &ode );
  return err;
}

/* 10. order of accuracy of deferred correction */
void assert_ode_dc_order(void)
{
  /* Richardson's extrapolation raises the order of accuracy p by one, or by two for symmetric methods
   * (trapezoidal method). The time steps are chosen so that they are not divided by the step size control
   * and the errors are not dominated by the round-off error or the tolerance of the nonlinear equation
   * solver. */
  const struct{
    method_t method;
    double dt;
    int order;
  } test[] = {
    { METHOD( Euler,  0, 1 ), 0.001, 2 },
    { METHOD( BEuler, 0, 1 ), 0.001, 2 },
    { METHOD( Heun,   0, 3 ), 0.04,  4 },
    { METHOD( RK4,    0, 4 ), 0.04,  5 },
    { METHOD( RKG,    0, 4 ), 0.04,  5 },
    { METHOD( TR,     0, 2 ), 0.02,  4 },
    { METHOD( Radau,  0, 3 ), 0.04,  4 },
  };
  const double T = 2.0;
  int k;
  double order;
  bool result = true;

  for( k=0; k<sizeof(test)/sizeof(test[0]); k++ ){
    order = observed_order( circle_integrate_dc( &test[k].method, test[k].dt, T ), circle_integrate_dc( &test[k].method, 0.5*test[k].dt, T ) );
    if( fabs( order - test[k].order ) > ORDER_TOL ){
      eprintf( "%s: observed order = %g with deferred correction (expected %d)\n", test[k].method.name, order, test[k].order );
      result = false;
    }
  }
  zAssert( zODEUpdateDC (order of accuracy), result );
}

/* x' = -y, y' = x (circle) with counting the number of calls */
zVec circle_count(double t, zVec x, void *count, zVec dx)
{
  (*(int *)count)++;
  return circle( t, x, NULL, dx );
}

/* integrate a circle with the radius r from t=0 to t=1 with an adaptive step size control, and return the
 * relative error at t=1. The number of calls of the differential function is stored where count points. */
double circle_count_integrate(zODE *ode, bool is_dc, double r, int *count)
{
  const double dt = 0.1;
  const int n = 10;
  zVec x;
  int i;
  double err;

  if( is_dc ) zODECreateDC( ode, 2, 0, circle_count );
  else        zODECreate( ode, 2, 0, circle_count );
  x = zVecCreateList( 2, r, 0.0 );
  for( *count=0, i=0; i<n; i++ )
    is_dc ? zODEUpdateDC( ode, i*dt, x, dt, count ) : zODEUpdate( ode, i*dt, x, dt, count );
  err = circle_error( x, n*dt, r ) / r;
  zVecFree( x );
  if( is_dc ) zODEDestroyDC( ode ); else zODEDestroy( ode );
  return err;
}

/* 11. scale invariance of adaptive step size control */
void assert_ode_adaptive_scale(void)
{
  const double r[] = { 1.0e6, 1.0e12 };
  zODE ode;
  int i, count1, count;
  double err1, err;
  bool is_dc, result = true;

  /* the tolerance is relative for large states, so that the time step is divided in the same way as
   * for the unit circle. */
  for( is_dc=false; ; is_dc=true ){
    if( is_dc ) zODEAssign( &ode, RK4, NULL, NULL ); else zODEAssign( &ode, DP45, NULL, NULL );
    err1 = circle_count_integrate( &ode, is_dc, 1.0, &count1 );
    for( i=0; i<sizeof(r)/sizeof(double); i++ ){
      err = circle_count_integrate( &ode, is_dc, r[i], &count );
      if( count > count1 || err > 10 * err1 ){
        eprintf( "%s: radius %g, relative error = %g, %d calls (radius 1: %g, %d calls)\n", is_dc ? "RK4 with DC" : "DP45", r[i], err, count, err1, count1 );
        result = false;
      }
    }
    if( is_dc ) break;
  }
  zAssert( zODEUpdate (scale invariance of adaptive step size control), result );
}

/* integrate the unit circle from t=0 to t=n dt with an adaptive step size control, and return the error at
 * the end. The number of calls of the differential function is stored where count points. */
double circle_integrate_adaptive(const method_t *method, bool is_dc, double dt, int n, int *count)
{
  zODE ode;
  zVec x;
  int i;
  double err;

  method_assign( &ode, method );
  if( is_dc ) zODECreateDC( &ode, 2, 0, circle_count );
  else        zODECreate( &ode, 2, 0, circle_count );
  x = zVecCreateList( 2, 1.0, 0.0 );
  for( *count=0, i=0; i<n; i++ )
    is_dc ? zODEUpdateDC( &ode, i*dt, x, dt, count ) : zODEUpdate( &ode, i*dt, x, dt, count );
  err = circle_error( x, n*dt, 1.0 );
  zVecFree( x );
  if( is_dc ) zODEDestroyDC( &ode ); else zODEDestroy( &ode );
  return err;
}

/* 12. adaptive step size control for backward integration */
void assert_ode_adaptive_backward(void)
{
  const method_t erk[] = { METHOD( RKF45, 0, 5 ), METHOD( CK45, 0, 5 ), METHOD( DP45, 0, 5 ) };
  const double dt = 0.5;
  const int n = 20;
  int k, count_f, count_b;
  double err_f, err_b;
  bool result = true;

  /* since the circle is symmetric with respect to the time reversal, the time step has to be divided in
   * the same way for the backward integration as for the forward integration. The numbers of calls and
   * the errors might slightly differ for implicit methods due to the round-off errors in the iterations
   * to solve the nonlinear equations. */
  for( k=0; k<METHOD_ONESTEP_NUM+sizeof(erk)/sizeof(method_t); k++ ){
    const method_t *method = k < METHOD_ONESTEP_NUM ? &method_onestep[k] : &erk[k-METHOD_ONESTEP_NUM];
    bool is_dc = k < METHOD_ONESTEP_NUM;
    err_f = circle_integrate_adaptive( method, is_dc, dt, n, &count_f );
    err_b = circle_integrate_adaptive( method, is_dc, -dt, n, &count_b );
    if( abs( count_b - count_f ) > 0.01 * count_f || fabs( err_b - err_f ) > 1.0e-3 * err_f ){
      eprintf( "%s%s: backward error = %g, %d calls (forward error = %g, %d calls)\n", method->name, is_dc ? " with DC" : "", err_b, count_b, err_f, count_f );
      result = false;
    }
  }
  zAssert( zODEUpdate (adaptive step size control for backward integration), result );
}

/* methods with adaptive step size control: embedded Runge-Kutta methods, and one-step methods with
 * deferred correction */
#define METHOD_ADAPTIVE_NUM ( 3 + METHOD_ONESTEP_NUM )
const method_t *method_adaptive(int k, bool *is_dc)
{
  static const method_t erk[] = { METHOD( RKF45, 0, 5 ), METHOD( CK45, 0, 5 ), METHOD( DP45, 0, 5 ) };

  if( ( *is_dc = ( k >= 3 ) ) ) return &method_onestep[k-3];
  return &erk[k];
}

/* 13. prediction of the time step for adaptive step size control */
void assert_ode_step_prediction(void)
{
  const double T = 20.0, tol = 1.0e-4;
  const method_t *method;
  int k, count_s, count_l;
  double err_s, err_l;
  bool is_dc, result = true;

  /* since the time step is predicted from the last step, it does not depend on the interval of updates,
   * so that a long interval does not cost more evaluations than a short interval does. The numbers of
   * calls might slightly differ for implicit methods due to the iterations to solve the nonlinear
   * equations. */
  for( k=0; k<METHOD_ADAPTIVE_NUM; k++ ){
    method = method_adaptive( k, &is_dc );
    err_s = circle_integrate_adaptive( method, is_dc, 1.0, (int)T, &count_s );
    err_l = circle_integrate_adaptive( method, is_dc, 10.0, (int)T/10, &count_l );
    if( count_l > 1.02 * count_s || err_s > tol || err_l > tol ){
      eprintf( "%s%s: interval 10: error = %g, %d calls, interval 1: error = %g, %d calls\n", method->name, is_dc ? " with DC" : "", err_l, count_l, err_s, count_s );
      result = false;
    }
  }
  zAssert( zODEUpdate (prediction of the time step for adaptive step size control), result );
}

/* x' = 1 for t<=0.5, and not a number for t>0.5 */
zVec nan_after_half(double t, zVec x, void *count, zVec dx)
{
  (*(int *)count)++;
  zVecSetElemNC( dx, 0, t > 0.5 ? NAN : 1.0 );
  return dx;
}

/* 14. termination of adaptive step size control for a differential function which is not a number */
void assert_ode_adaptive_nan(void)
{
  const int count_max = 100000;
  const method_t *method;
  zODE ode;
  zVec x;
  int k, count;
  bool is_dc, result = true;

  /* the time step is shortened down to ZM_ODE_DT_TOL around t=0.5, and the rest of the interval is
   * abandoned after the state becomes not a number. Euler method never evaluates the function for t>0.5. */
  for( k=0; k<METHOD_ADAPTIVE_NUM; k++ ){
    method = method_adaptive( k, &is_dc );
    method_assign( &ode, method );
    if( is_dc ) zODECreateDC( &ode, 1, 0, nan_after_half );
    else        zODECreate( &ode, 1, 0, nan_after_half );
    x = zVecCreateList( 1, 0.0 );
    count = 0;
    is_dc ? zODEUpdateDC( &ode, 0, x, 1.0, &count ) : zODEUpdate( &ode, 0, x, 1.0, &count );
    if( count > count_max ){
      eprintf( "%s%s: x = %g, %d calls\n", method->name, is_dc ? " with DC" : "", zVecElemNC(x,0), count );
      result = false;
    }
    zVecFree( x );
    if( is_dc ) zODEDestroyDC( &ode ); else zODEDestroy( &ode );
  }
  zAssert( zODEUpdate (termination for a differential function which is not a number), result );
}

/* integrate the unit circle from t=0 to t=T with a multistep method whose history is created by
 * zODECreateHistory(), and return the error at T. the inconsistent initial history filled with zero
 * degrades the order of accuracy to zero for Gear's method and one for Adams's method. */
double circle_integrate_history(const method_t *method, int step, double dt, double T, bool *result)
{
  zODE ode;
  zVec x;
  double t = 0, err;
  int n;

  method_create( &ode, method, 2, step, circle );
  x = zVecCreateList( 2, 1.0, 0.0 );
  if( zODECreateHistory( &ode, &t, x, dt, NULL ) != x || fabs( t - ( step - 1 ) * dt ) > zTOL ){
    eprintf( "%s (%d steps): t = %g after the history is created (expected %g)\n", method->name, step, t, ( step - 1 ) * dt );
    *result = false;
  }
  for( n=(int)( ( T - t ) / dt + 0.5 ); n>0; n--, t+=dt )
    zODEUpdate( &ode, t, x, dt, NULL );
  err = circle_error( x, t, 1.0 );
  zVecFree( x );
  zODEDestroy( &ode );
  return err;
}

/* 15. order of accuracy of multistep methods with the history created by a starter */
void assert_ode_history_order(void)
{
  const method_t adams = METHOD( Adams, 0, 0 ), gear = METHOD( Gear, 0, 0 );
  const double dt = 0.025, T = 1.0;
  int step;
  double order;
  bool result1 = true, result2 = true;

  for( step=1; step<=6; step++ ){
    order = observed_order( circle_integrate_history( &adams, step, dt, T, &result2 ), circle_integrate_history( &adams, step, 0.5*dt, T, &result2 ) );
    if( fabs( order - step ) > ORDER_TOL ){
      eprintf( "Adams (%d steps): observed order = %g (expected %d)\n", step, order, step );
      result1 = false;
    }
  }
  /* the error saturates for more steps of Gear's method (see assert_ode_order_multistep()). */
  for( step=1; step<=4; step++ ){
    order = observed_order( circle_integrate_history( &gear, step, dt, T, &result2 ), circle_integrate_history( &gear, step, 0.5*dt, T, &result2 ) );
    if( fabs( order - step ) > ORDER_TOL ){
      eprintf( "Gear (%d steps): observed order = %g (expected %d)\n", step, order, step );
      result1 = false;
    }
  }
  zAssert( zODECreateHistory (order of accuracy of multistep methods), result1 );
  zAssert( zODECreateHistory (time after the history is created), result2 );
}

/* 16. history of one-step methods and stability of the history of Gear's method for a stiff equation */
void assert_ode_history_misc(void)
{
  const method_t gear = METHOD( Gear, 0, 0 );
  const double dt = 0.1; /* A dt = 10, with which explicit methods diverge */
  const int n = 50;
  zODE ode;
  zVec x;
  double t;
  int i, k;
  bool result1 = true, result2 = true;

  /* nothing has to be done for one-step methods. */
  for( k=0; k<METHOD_ONESTEP_NUM; k++ ){
    method_create( &ode, &method_onestep[k], 2, 0, circle );
    x = zVecCreateList( 2, 1.0, 0.0 );
    t = 0;
    zODEInitHistory( &ode, t, x, NULL );
    if( zODECreateHistory( &ode, &t, x, dt, NULL ) != x || t != 0 ||
        zVecElemNC(x,0) != 1.0 || zVecElemNC(x,1) != 0.0 ){
      eprintf( "%s: t = %g, x = ( %g, %g ) after the history is created\n", method_onestep[k].name, t, zVecElemNC(x,0), zVecElemNC(x,1) );
      result1 = false;
    }
    zVecFree( x );
    zODEDestroy( &ode );
  }
  /* the history of Gear's method has to be created without spoiling the stability for a stiff equation. */
  for( k=1; k<=6; k++ ){
    method_create( &ode, &gear, 1, k, stiff );
    x = zVecCreateList( 1, 1.0 );
    t = 0;
    zODECreateHistory( &ode, &t, x, dt, NULL );
    for( i=0; i<n; i++, t+=dt )
      zODEUpdate( &ode, t, x, dt, NULL );
    if( !( fabs( zVecElemNC(x,0) ) <= 1.0e-3 ) ){
      eprintf( "Gear (%d steps): |x| = %g after %d steps\n", k, fabs( zVecElemNC(x,0) ), n );
      result2 = false;
    }
    zVecFree( x );
    zODEDestroy( &ode );
  }
  zAssert( zODEInitHistory and zODECreateHistory (one-step methods), result1 );
  zAssert( zODECreateHistory (stability of Gear method for a stiff equation), result2 );
}

/* 17. Gear's method for a state much smaller than the history */
void assert_ode_gear_small_state(void)
{
  const method_t gear = METHOD( Gear, 0, 0 );
  /* coefficients of Gear's method, with which x_{n+1} + sum_i a_i x_{n-i} = b dt f(x_{n+1}) */
  const double a[][6] = {
    { -1, 0, 0, 0, 0, 0 },
    { -4.0/3, 1.0/3, 0, 0, 0, 0 },
    { -18.0/11, 9.0/11, -2.0/11, 0, 0, 0 },
    { -48.0/25, 36.0/25, -16.0/25, 3.0/25, 0, 0 },
    { -300.0/137, 300.0/137, -200.0/137, 75.0/137, -12.0/137, 0 },
    { -120.0/49, 150.0/49, -400.0/147, 75.0/49, -24.0/49, 10.0/147 },
  };
  const double b[] = { 1, 2.0/3, 6.0/11, 12.0/25, 60.0/137, 20.0/49 };
  const double dt = 0.1, x0 = 1.0e-15, tol = 1.0e-8;
  const int n = 10;
  double h[16], s;
  zODE ode;
  zVec x;
  int i, j, k;
  bool result = true;

  /* the history is filled with one, and the state starts from a much smaller value, as it happens to
   * a stiff equation. since f is linear, the solution has to follow the recurrence
   *  x_{n+1} = - sum_i a_i x_{n-i} / ( 1 + b dt A ). */
  for( k=2; k<=6; k++ ){
    method_create( &ode, &gear, 1, k, stiff );
    x = zVecCreateList( 1, 1.0 );
    zODEInitHistory( &ode, 0, x, NULL );
    for( i=0; i<k-1; i++ ) h[i] = 1.0;
    zVecSetElemNC( x, 0, ( h[k-1] = x0 ) );
    for( i=k-1; i<k-1+n; i++ ){
      for( s=0, j=0; j<k; j++ ) s += a[k-1][j] * h[i-j];
      h[i+1] = -s / ( 1 + b[k-1] * dt * STIFF_A );
      zODEUpdate( &ode, ( i - k + 1 ) * dt, x, dt, NULL );
      if( fabs( zVecElemNC(x,0) - h[i+1] ) > tol * fabs( h[i+1] ) ){
        eprintf( "Gear (%d steps): x = %g at step %d (expected %g)\n", k, zVecElemNC(x,0), i-k+2, h[i+1] );
        result = false;
        break;
      }
    }
    zVecFree( x );
    zODEDestroy( &ode );
  }
  zAssert( zODEUpdate (Gear method for a state much smaller than the history), result );
}

/* integrate the oscillator for n steps by the regular update of zODE2 assigned in advance, and return x. */
double ode2_regular_integrate(zODE2 *ode, int num, double dt, int n)
{
  zVec x, v;
  double result;
  int i;

  zODE2Create( ode, 1, num, oscillator );
  x = zVecCreateList( 1, 1.0 );
  v = zVecCreateList( 1, 0.0 );
  for( i=0; i<n; i++ )
    zODE2Update( ode, i*dt, x, v, dt, NULL );
  result = zVecElemNC(x,0);
  zVecFreeAtOnce( 2, x, v );
  zODE2Destroy( ode );
  return result;
}

/* 18. regular update of zODE2 with Gear method, which is identical with backward Euler method for one step.
 * the concatenation of the state in the nonlinear equation of Gear method has to receive the utility
 * given to zODE2Update(). */
void assert_ode2_gear(void)
{
  const double dt = 0.01;
  const int n = 100;
  zODE2 ode;
  double x_gear, x_beuler;
  bool result = true;

  zODE2Assign( &ode, Regular, NULL, NULL, NULL, NULL );
  zODE2AssignRegular( &ode, Gear );
  x_gear = ode2_regular_integrate( &ode, 1, dt, n );
  zODE2Assign( &ode, Regular, NULL, NULL, NULL, NULL );
  zODE2AssignRegular( &ode, BEuler );
  x_beuler = ode2_regular_integrate( &ode, 0, dt, n );
  if( !zIsTol( x_gear - x_beuler, 1.0e-8 ) ){
    eprintf( "Gear: x = %.17g (backward Euler: %.17g)\n", x_gear, x_beuler );
    result = false;
  }
  zAssert( zODE2Update (regular update with Gear method), result );
}

/* error of a harmonic oscillation x'' = -x with x(0) = 1 and x'(0) = 0 at the time 1, integrated with
 * the history initialized (is_create = false) or created (is_create = true). t0 is set for the time
 * after the history is created. */
double ode2_history_error(zODE2 *ode, int num, bool is_create, double dt, double *t0)
{
  zVec x, v;
  double t = 0, result;

  zODE2Create( ode, 1, num, oscillator );
  x = zVecCreateList( 1, 1.0 );
  v = zVecCreateList( 1, 0.0 );
  if( is_create )
    zODE2CreateHistory( ode, &t, x, v, dt, NULL );
  else
    zODE2InitHistory( ode, t, x, v, NULL );
  *t0 = t;
  for( ; t<1-0.5*dt; t+=dt )
    zODE2Update( ode, t, x, v, dt, NULL );
  result = sqrt( zSqr( zVecElemNC(x,0) - cos(t) ) + zSqr( zVecElemNC(v,0) + sin(t) ) );
  zVecFreeAtOnce( 2, x, v );
  zODE2Destroy( ode );
  return result;
}

/* observed order of accuracy of a second-order ODE solver assigned in advance with the history
 * initialized or created. */
double ode2_history_order(zODE2 *ode, int num, bool is_create, double *t0)
{
  double e1, e2;

  e1 = ode2_history_error( ode, num, is_create, 0.025, t0 );
  e2 = ode2_history_error( ode, num, is_create, 0.0125, t0 );
  return log2( e1 / e2 );
}

/* 19. order of accuracy of second-order ODE solvers with the history created */
void assert_ode2_history(void)
{
  zODE2 ode;
  double p, t0;
  int k;
  bool result = true;

  for( k=2; k<=4; k++ ){
    zODE2Assign( &ode, Regular, NULL, NULL, NULL, NULL );
    zODE2AssignRegular( &ode, Adams );
    p = ode2_history_order( &ode, k, true, &t0 );
    if( p < k - 0.2 || !zIsTiny( t0 - (k-1)*0.0125 ) ){
      eprintf( "Regular Adams (%d steps): order %g, t = %g after the history is created\n", k, p, t0 );
      result = false;
    }
    zODE2AssignRegular( &ode, Gear );
    p = ode2_history_order( &ode, k, true, &t0 );
    if( p < k - 0.2 || !zIsTiny( t0 - (k-1)*0.0125 ) ){
      eprintf( "Regular Gear (%d steps): order %g, t = %g after the history is created\n", k, p, t0 );
      result = false;
    }
  }
  /* the velocity at a half step before has to be created with the second-order accuracy. */
  zODE2Assign( &ode, Leapfrog, NULL, NULL, NULL, NULL );
  p = ode2_history_order( &ode, 0, false, &t0 );
  if( p > 1.2 ){
    eprintf( "Leapfrog: order %g with the history initialized\n", p );
    result = false;
  }
  p = ode2_history_order( &ode, 0, true, &t0 );
  if( p < 1.8 || t0 != 0 ){
    eprintf( "Leapfrog: order %g, t = %g with the history created\n", p, t0 );
    result = false;
  }
  /* no history */
  zODE2Assign( &ode, Symplectic, NULL, NULL, NULL, NULL );
  ode2_history_order( &ode, 0, true, &t0 );
  if( t0 != 0 ){
    eprintf( "Symplectic: t = %g with the history created\n", t0 );
    result = false;
  }
  zAssert( zODE2CreateHistory (order of accuracy of second-order ODE solvers), result );
}

/* acceleration linear in time */
zVec linear_acceleration(double t, zVec x, zVec dx, void *dummy, zVec ddx)
{
  zVecSetElemNC( ddx, 0, t );
  return ddx;
}

/* 20. velocity of leapfrog method at a full step for an acceleration varying in time. since the velocity
 * is the trapezoidal integration of the acceleration, it is exact for an acceleration linear in time,
 * if the acceleration at the end of the step is evaluated at the correct time. */
void assert_ode2_leapfrog_time(void)
{
  const double dt = 0.01;
  const int n = 100;
  zODE2 ode;
  zVec x, v;
  double t = 0;
  int i;
  bool result = true;

  zODE2Assign( &ode, Leapfrog, NULL, NULL, NULL, NULL );
  zODE2Create( &ode, 1, 0, linear_acceleration );
  x = zVecCreateList( 1, 0.0 );
  v = zVecCreateList( 1, 0.0 );
  zODE2CreateHistory( &ode, &t, x, v, dt, NULL );
  for( i=0; i<n; i++, t+=dt ){
    zODE2Update( &ode, t, x, v, dt, NULL );
    if( !zIsTol( zVecElemNC(v,0) - 0.5*zSqr(t+dt), 1.0e-12 ) ){
      eprintf( "Leapfrog: v = %.17g at t = %g (exact: %.17g)\n", zVecElemNC(v,0), t+dt, 0.5*zSqr(t+dt) );
      result = false;
      break;
    }
  }
  zVecFreeAtOnce( 2, x, v );
  zODE2Destroy( &ode );
  zAssert( zODE2Update (time of the acceleration at the end of a step of leapfrog method), result );
}

/* x0' = 0 and x1' = -x1 + sqrt(1-x0^2), which is NaN for x0 > 1. */
zVec boundary_ode(double t, zVec x, void *dummy, zVec dx)
{
  zVecSetElemNC( dx, 0, 0 );
  zVecSetElemNC( dx, 1, -zVecElemNC(x,1) + sqrt( 1 - zSqr( zVecElemNC(x,0) ) ) );
  return dx;
}

/* 21. Gear's method for a state close to the boundary of the domain of the differential function.
 * the numerical Jacobian matrix of the nonlinear equation has to be computed by the one-sided difference,
 * since the differential function is not finite on one side. the exact solution is x0 = const. and
 * x1 = c + ( x1(0) - c ) exp(-t) for c = sqrt(1-x0^2). */
void assert_ode_gear_boundary(void)
{
  const method_t gear = METHOD( Gear, 0, 0 );
  const double dt = 0.01, x00 = 1 - 1.0e-9, tol = 1.0e-2;
  const int n = 100;
  double c, x1;
  zODE ode;
  zVec x;
  int i, k;
  bool result = true;

  c = sqrt( 1 - x00*x00 );
  x1 = c + ( 1 - c ) * exp( -n*dt );
  for( k=1; k<=6; k++ ){
    method_create( &ode, &gear, 2, k, boundary_ode );
    x = zVecCreateList( 2, x00, 1.0 );
    zODEInitHistory( &ode, 0, x, NULL );
    for( i=0; i<n; i++ )
      zODEUpdate( &ode, i*dt, x, dt, NULL );
    if( !zVecIsFinite( x ) || !( fabs( zVecElemNC(x,0) - x00 ) < 1.0e-10 ) || !( fabs( zVecElemNC(x,1) - x1 ) < tol ) ){
      eprintf( "Gear (%d steps): x = (%.17g, %g) (expected (%.17g, %g))\n", k, zVecElemNC(x,0), zVecElemNC(x,1), x00, x1 );
      result = false;
    }
    zVecFree( x );
    zODEDestroy( &ode );
  }
  zAssert( zODEUpdate (Gear method close to the boundary of the domain), result );
}

int main(void)
{
  assert_ode_order_onestep();
  assert_ode_order_multistep();
  assert_ode_circle_norm();
  assert_ode_stiff();
  assert_ode_implicit_zero();
  assert_ode_zero_dt();
  assert_ode2_symplectic();
  assert_ode2_zero_dt();
  assert_ode_dc_accuracy();
  assert_ode_dc_order();
  assert_ode_adaptive_scale();
  assert_ode_adaptive_backward();
  assert_ode_step_prediction();
  assert_ode_adaptive_nan();
  assert_ode_history_order();
  assert_ode_history_misc();
  assert_ode_gear_small_state();
  assert_ode2_gear();
  assert_ode2_history();
  assert_ode2_leapfrog_time();
  assert_ode_gear_boundary();
  return EXIT_SUCCESS;
}
