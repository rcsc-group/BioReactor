/**
#Lagrangian Particles
 */
#include "run.h"

coord * loc;
int n_part;
bool part_linear = true;
bool ssi = false;
/**
Tracer particles are advected using a simple forward Euler scheme or,
in 2D, the x/y split semi-implicit method can be used by setting the
`global bool`: `ssi = true`. Learn more on the [test of some schemes
page](timetest.c).
 */
event init (t = 0);

event set_dtmax (i++);

event advance_particles(i++) {
  coord mind = {X0, Y0, Z0}; 
  double V[dimension*n_part];
  interpolate_array ((scalar*){u}, loc, n_part, V, part_linear);
  bool even = true;
  double * ps = (double*)malloc(sizeof(double));
  if (ssi) {    // We make a scratch to store the "old" position 
    if (i % 2 != 0)
      even = false;
    ps = realloc(ps, n_part*sizeof(double));
    for (int j = 0; j < n_part; j++){
      if (even){
        ps[j] = loc[j].x;
      }else{
        ps[j] = loc[j].y;
      }
    }
    
  }
  if (pid() == 0) { // Works in parallel, but is not parallelized
    for (int j = 0; j < n_part; j++) {
      int k = 0;
      foreach_dimension() {//Advance
        if (V[dimension*j + k] < HUGE)
          loc[j].x += dt*V[dimension*j + k];
	if (loc[j].x < mind.x)
	  loc[j].x += L0;
        else if (loc[j].x > (mind.x + L0))
          loc[j].x -= L0;
        k++;
      }
    }
  }
#if _MPI // All theads should know the new coordinates
  MPI_Bcast (&loc[0], (sizeof(coord)*n_part)/sizeof(double),
             MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  /**
     For 2D we can choose to do x/y-split semi-implicit time
     integration. The forward method is used as an initial guess.
   */
  if (ssi){ 
    double tolerance = 1e-5;
    double error = -1., errorj = 0;
    int it = 0;
    int itmax = 50;
    do{ //Newton-Rapshon
      interpolate_array ((scalar*){u}, loc, n_part, V, part_linear);
      if (pid() == 0) {
	for (int j = 0; j < n_part; j++) {
	  if (even) { // x backward
	    double pst = loc[j].x; //Store current guess 
	    loc[j].x = ps[j] + dt*V[dimension*j]; //Update old x position
	    errorj = fabs(pst - loc[j].x); //Compare against the stored guess
	  }else{ // y backward, same deal
	    double pst = loc[j].y; 
	    loc[j].y = ps[j] + dt*V[dimension*j + 1]; 
	    errorj = fabs(pst - loc[j].y);
	  }
	  if (errorj > error) // Find largest error
	    error = errorj; 
        }
      }
#if _MPI // All theads should know the new coordinates
      MPI_Bcast (&loc[0], (sizeof(coord)*n_part)/sizeof(double),
                 MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
      it++;
    }while (error < tolerance && it < itmax);
   // if (it >= itmax)
   //   fprintf(stderr, "Particles: N-R procedure did not converge i = %d, t = %g\n", i ,t);
    for (int j = 0; j < n_part; j++){ //Boundaries
      foreach_dimension(){
	if (loc[j].x < mind.x)
	  loc[j].x += L0;
	else if (loc[j].x > (mind.x + L0))
	  loc[j].x -= L0;
      } //Dimensions
    } // particle iterator
  } //if ssi
  free (ps);
} //event

event free_particles (t = end, last)
  free (loc);

/**
## User function for particle seeding

Initialize particles in a 2D `l`$\times$`l` grid centered at {`xo, yo`} with `n`$\times$`n` particles. 
*/
void init_particles_2D_square_grid (int n, double xo, double yo, double l){
  n_part = sq(n);
  loc = malloc (n_part*sizeof(coord));
  int i = 0;
  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      loc[i].x = xo - l/2. + (double)j*(l/((double)n - 1.));
      loc[i].y = yo - l/2. + (double)k*(l/((double)n - 1.));
      i++;
    }
  }
}

/**
Initialize 1 single particle at xo, yo 
*/
void init_singleparticle (double xo, double yo){
  n_part = 1;
  loc = malloc (n_part*sizeof(coord));
  
  loc[0].x=xo;
  loc[0].y=yo;

}

void init_threeparticles (double xo1, double yo1, double xo2, double yo2, double xo3, double yo3){
  n_part = 3;
  loc = malloc (n_part*sizeof(coord));
  
  loc[0].x=xo1;
  loc[0].y=yo1;

  loc[1].x=xo2;
  loc[1].y=yo2;

  loc[2].x=xo3;
  loc[2].y=yo3;

}

void init_varnumparticles(int numparts,...){

  va_list valist;
  va_start(valist, numparts*2);

  n_part=numparts;
  loc = malloc (n_part*sizeof(coord));

  for (int i = 0; i < numparts; i++)
  {
    loc[i].x=va_arg(valist, double);
    loc[i].y=va_arg(valist, double);
  }

  va_end(valist);
  
}

/**
Initialize particles in a 2D `lx`$\times$`ly` grid centered at {`xo, yo`} with `nx`$\times$`ny` particles. 
*/
void init_particles_2D_rectangular_grid (int nx, int ny, double xo, double yo, double lx, double ly){
  n_part = nx*ny;
  loc = malloc (n_part*sizeof(coord));
  int i = 0;
  for (int j = 0; j < nx; j++) {
    for (int k = 0; k < ny; k++) {
      loc[i].x = xo - lx/2. + (double)j*(lx/((double)nx - 1.));
      loc[i].y = yo - ly/2. + (double)k*(ly/((double)ny - 1.));
      i++;
    }
  }
}


void init_particles_2D_semicircular_arcs (int num_arc, int num_per_arc, double xo, double yo, double D){
  n_part = num_arc*num_per_arc + 1;
  loc = malloc (n_part*sizeof(coord));
  int i = 0;
  for (int j = 0; j < num_arc; j++){
    for (int k = 0; k < num_per_arc; k++){
      loc[i].x = xo + (D/2.)*((double)(j+1)/(double)num_arc)*cos(k*pi/(num_per_arc-1));
      loc[i].y = (D/2.)*((double)(j+1)/(double)num_arc)*sin(k*pi/(num_per_arc-1));
      i++;
    }
  }
  loc[i].x = xo;
  loc[i].y = 0;
}

/**
Initialize particles in a 2D `lx`$\times$`ly` grid centered at {`xo, yo`} with `nx`$\times$`ny` particles. 
*/
void init_particles_2D_rectangular_grid_and_semicirculararc (int nx, int ny, double xo, double yo, double lx, double ly,
 int num_arc, int num_per_arc, double xoo, double yoo, double D){
  
  n_part = nx*ny + num_arc*num_per_arc + 1;
  loc = malloc (n_part*sizeof(coord));
  
  int i = 0;
  
  for (int j = 0; j < nx; j++) {
    for (int k = 0; k < ny; k++) {
      loc[i].x = xo - lx/2. + (double)j*(lx/((double)nx - 1.));
      loc[i].y = yo - ly/2. + (double)k*(ly/((double)ny - 1.));
      i++;
    }
  }

  for (int jj = 0; jj < num_arc; jj++){
    for (int kk = 0; kk < num_per_arc; kk++){
      loc[i].x = xoo + (D/2.)*((double)(jj+1)/(double)num_arc)*cos(kk*pi/(num_per_arc-1));
      loc[i].y = (D/2.)*((double)(jj+1)/(double)num_arc)*sin(kk*pi/(num_per_arc-1));
      i++;
    }
  }
  loc[i].x = xoo;
  loc[i].y = 0;
}

/**
The following function places `n` particles randomly in a circle with 
radius `R` at location {`xo, yo`}. 
*/
void init_particles_random_circle (int n, double xo, double yo, double R) {
  n_part = n;
  loc = malloc (n_part*sizeof(coord));
  int j = 0;
  if (pid() == 0) {
    while (j < n_part) {
      double a = noise();
      double b = noise();
      if (sq(a) + sq(b) < R) {
	loc[j].x = a + xo;
        loc[j].y = b + yo;
	j++;
      }
    }
  }
#if _MPI // All theads should know the randomized coordinates
  MPI_Bcast (&loc[0], (sizeof(coord)*n_part)/sizeof(double),
             MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}
/**
The following function initializes a particle at the centre of each grid cell. 
*/
void init_particles_in_cells(){
  int n = 0;
#if _MPI
#if TREE
  scalar index[];
  z_indexing (index, true);
  stats ff = statsf (index);
  n_part = (int)(ff.max + 1.5);
  loc = malloc (n_part*sizeof(coord));
  while (n < n_part) {
    foreach_dimension()
      loc[n].x = HUGE;
    n++;
  }
  foreach(){
    coord cc = {x, y, z};
    foreach_dimension()
      loc[(int)(index[] + 0.5)].x = cc.x;
  }
#else // MG _MPI
  foreach()
    n++;
  n_part = npe()*n;
  loc = malloc (n_part*sizeof(coord));
  int j = 0;
  while (j < n_part) {
    foreach_dimension()
      loc[j].x = HUGE;
    j++;
  }
  j = 0;
  foreach() {
    coord cc = {x, y, z};
    foreach_dimension()
      loc[j + pid()*n].x = cc.x;
    j++;
  }
#endif // _MPI && (TREE || MG)
  MPI_Allreduce (MPI_IN_PLACE, &loc[0], (n_part*sizeof(coord))/sizeof(double),
                 MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else // !_MPI, it is more concise:
  n_part = 0;
  foreach()
    n_part++;
  loc = malloc (n_part*sizeof(coord));
  foreach() {
    coord cc = {x, y, z};
    foreach_dimension()
      loc[n].x = cc.x;
    n++;
  }
#endif
}
/**
For visualization with bview, the `scatter()` function can 
be used when `BVIEW` $\neq 0$.  
*/ 
#if BVIEW
#include "myscatter.h"
#endif

/**
## Test
* [A test for the x/y-split semi-implicit scheme](parttest.c)

## Usage

* [Tag a portion of a fluid](splash.c)
* [LES of isotropic turbulence](isotropicLES.c)
* [Laminar mixing in 2D](laminarmixing.c)
* [Axisymmetric mixing](coffee.c)
*/
