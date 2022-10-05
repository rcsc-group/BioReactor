#include "navier-stokes/centered.h"
#define mu(f) ( 1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)  )
//#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
//#include "tracer.h"
#include "utils.h"
//#include "henry.h"
#include "henry_oxy2.h"
//#include "henry_oxy.h"
//#include "view.h"
#include "view.h"                    // need to make the animations
#include "draw.h"                    // visualisation helper
#include "tag.h"                     // helps track droplet properties

#include "myparticles.h"
#include "myscatter.h"

const double NN      = 256;      // resolution
const double t_change= 11.2278;  // change in rocking motion
//const double th_cont = 90;     // contact angle
const double t_mix   = 9.5;    // time for releasing tracers
const double t_end_file = 9.5; // t<t_end_file; saving freq=dt_file
const double t_end   = 50.0;    // final time
const double dt_file = 0.1123;   // saving frq before t_end_file
const double dt_file_refine = 0.5; // saving frq after t_end_file
const double dt_video= 0.01;

///*
//&&&&&& Bioreactor parameters &&&&&//
const double L_bio  = 0.4;     // m - reference length scale
const double Th_max = 0.122;   // rad - maximum angle of rotation - 7 degree
const double T_per  = 1.5;     // sec - 40rpm-1.5; 20rpm-3;
//*/


// &&&&& Mesh refinement &&&&&//
int    MINLEVEL = 7;
int    MAXLEVEL = 9;
double F_MAX    = 1e-4;
double U_MAX    = 1e-2;


double LL = 1.0;    // width
double Ly = 0.44;   // height
double L_piv = 0.22; // distance to pivot point
int    i_fig = 500;

scalar ft[], c[], oxy[];   // for tracer and oxygen transfer
scalar * tracers = {ft}, * stracers = {c,oxy};
char buf1[50], buf2[50], buf3[20], buf4[20], buf5[20], buf6[20];
double (* gradient) (double, double, double) = minmod2;
double U0, Ub, Re_w, Re_a, We_w, Fr, rhor, mur, Pe_tracer_1, Pe_tracer_2, Pe_oxy_1, Pe_oxy_2, Th, Th_d, Th_2d, U_bio, w_bio, w_bio_st, T_bio, Th_max2;

//&&&&& Material properties &&&&&//
const double rho_w = 1.0e3;   // kg/m^3
const double rho_a = 1.225; // 0.01e3;   // kg/m^3
const double mu_w  = 1.0e-3;  // Pa-s @ 20C
const double mu_a  = 1.0e-4; // Pa-s - why multiplied 10?
const double grav  = 9.8;     // m/s^2
const double sigma = 0.0728;  // N/m
const double D_tracer_1 = 0.44e-9; // concentration = 1; diffusion in water
const double D_tracer_2 = 1.0e-30; // concentration = 0; diffusion in air
const double D_oxy_1 = 1.98e-5; // m^2/s; concentration = 1; diffusion in air
const double D_oxy_2 = 1.90e-9; // m^2/s; concentration = 0; diffusion in water
const double c_tracer_alpha = 1.0e30; // concentration 1 to concentration 0
const double c_oxy_alpha    = 30;     // c_l=alpha*c_g; c1=alpha*c2; solubility; concentration 1 to concentration 0; air to water = 1/30;

FILE * fp_stats;

//&&&&& Contact angle &&&&&//
/*
vector h[];
h.t[left]  = contact_angle(th_cont*pi/180);
h.t[right] = contact_angle(th_cont*pi/180);
*/

int main(){
  
  //***** Initial conditions *****// 
  L0 = LL;                 // domain size
  //mask (y > 0.9 ? top:none);
  
  origin (-L0/2., -L0/2.); // origin
  init_grid(NN);           // resolution
  
  double H_bio,V_bio;
  H_bio  = L_bio*Ly;
  V_bio  = L_bio/4*(H_bio + 0.5*L_bio*tan(Th_max));
  
  //U_bio  = L_bio*Th_max/T_per;  // reference velocity m/s
  U_bio  = V_bio/(H_bio*0.5)/T_per;
  T_bio  = L_bio/U_bio;
  w_bio  = 2*pi/T_per;          // angular velocity rad/s
  w_bio_st = w_bio*T_bio;       // dimensionless w_bio = 2*pi
  U0     = w_bio_st*Th_max;     // initial velocity (clockwise)
  Ub     = 0.;                  // Boundary velocity

  /*
  dt_file  = (T_per/T_bio)/5;   // saving frq before t_end_file
  dt_file_refine = (T_per/T_bio)/10; // saving frq after t_end_file
  dt_video = (T_per/T_bio)/20;
  */

  //***** Dimensionless number *****//
  Re_w = rho_w*U_bio*L_bio/mu_w;
  Re_a = rho_a*U_bio*L_bio/mu_a;
  We_w = rho_w*U_bio*U_bio*L_bio/sigma;
  Fr   = U_bio/sqrt(grav*L_bio);
  rhor = rho_a/rho_w;
  mur  = mu_a/mu_w;
  Pe_tracer_1 = U_bio*L_bio/D_tracer_1;
  Pe_tracer_2 = U_bio*L_bio/D_tracer_2;
  Pe_oxy_1    = U_bio*L_bio/D_oxy_1;
  Pe_oxy_2    = U_bio*L_bio/D_oxy_2;
    
  rho1 = 1.0;
  rho2 = rho1*rhor;
  mu1  = 1./Re_w;
  mu2  = mur*mu1;
  f.sigma = 1./We_w;

  // Contact angle
  //f.height = h;
  
  // D1 for c=1; D2 for c=0; c1=alpha*c2
  c.D1 = 1./Pe_tracer_1;
  c.D2 = 1./Pe_tracer_2;
  c.alpha = c_tracer_alpha;  // c1 = alpha*c2
  oxy.D1 = 1./Pe_oxy_1;
  oxy.D2 = 1./Pe_oxy_2;
  oxy.alpha = c_oxy_alpha;

  //&&&&& Boundary conditions &&&&&//
  u.n[left]  = dirichlet(0.);
  u.t[left]  = dirichlet(0.);
  u.n[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

  //***** Tracers: slope limiter *****//
  //ft.gradient = minmod2;
  c.gradient = minmod2;
  oxy.gradient = minmod2;

  // Iterations
  NITERMAX = 500;
  TOLERANCE = 1.0e-5;
  
  CFL = 0.02;

  // Pointer of the file to save stats
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");
  }

  run();

  fclose(fp_stats);
}

scalar omega[], viewingfield[], viewingfieldmain[], mylevel[], velnorm[];

//***** Initial conditions  *****//
event init (t = 0)
{
  // Initial phases: y>0 -> volume fraction=0 (air); y<0 -> volume fraction=1
  fraction (f, -y);

  foreach()
  {
	omega[] = 0.0;
  }

  ssi=true;
  //init_particles_2D_semicircular_arcs(1, 30, (double)S0, 0, 1.0*(double)D0);
  init_particles_2D_rectangular_grid (16, 16, 0.0, 0.0, 0.9, 0.9);

  // solid
  // solid (cs,fs, intersection( -(y-0.5*Ly), -(-y-0.5*Ly) ) );

  // Tracer: positive -> ft=1; negative -> ft=0
  //fraction (ft, -(sq(x-0) + sq(y+0.25*LL) - sq(LL*0.1)) );

  // Oxygen: positive -> c=1; negative -> c=0
  // c.phi1 -> f(volume fraction); c.phi2 -> 1-f(1-volume fraction)
  // c.D1 -> volume fraction=1 -> water in bioreactor
  // c.D2 -> volume fraction=0 -> air in bioreactor
  // fraction (c, -(sq(x-0) + sq(y+0.25*LL) - sq(0.1*LL)) );  
  
  // initial velocity field
  // solid vol fraction = 1 (no solid) -> Ux0, solid vol fraction = 0 (w/solid) -> u.x = 0
  //foreach(){
  //u.x[] =  U0*(y+L_piv*cos(Th_max));
  // u.y[] = -U0*(x+L_piv*sin(Th_max));
    //u.x[] =  U0*y;
    //u.y[] = -U0*x;
  //}
}

event tracer(t = t_mix){
  fraction (c, -(sq(x-0) + sq(y+Ly*0.5*0.5) - sq(0.084*Ly)) );

  foreach(){
    if (f[] < 1.)
      oxy[] = 1-f[];
  }
}

event replenish_oxy(t = t_mix; t += 0.1){

  foreach(){
    if (f[] < 1.)
      oxy[] = 1-f[];
  }
}

//&&&&& acceleration &&&&&&//
event acceleration(i++)
{
  ///*
  // Rotation
  if (t >= t_change){
    Th    = Th_max*sin(w_bio_st*t);
    Th_d  = w_bio_st*Th_max*cos(w_bio_st*t);
    Th_2d = -w_bio_st*w_bio_st*Th_max*sin(w_bio_st*t);
  }
  else if (t < t_change){
    Th_max2 = (Th_max-0)/(t_change-0)*t;
    Th    = Th_max2*sin(w_bio_st*t);
    Th_d  = w_bio_st*Th_max2*cos(w_bio_st*t);
    Th_2d = -w_bio_st*w_bio_st*Th_max2*sin(w_bio_st*t);
  }
  //*/
  ///* 
  face vector av = a;
  // 1st: gravitational force, 2nd: Coriolis force
  // 3rd: centrifugal force,   4th: azimuthal force, 5th: no traslational force
  foreach_face(x)
    av.x[] = -sin(Th)/(Fr*Fr) + 2*Th_d*face_value(u.y,0)	\
    + Th_d*Th_d*(x+L_piv*sin(Th)) + Th_2d*(y+L_piv*cos(Th));
  foreach_face(y)
    av.y[] = -cos(Th)/(Fr*Fr) - 2*Th_d*face_value(u.x,0)	\
    + Th_d*Th_d*(y+L_piv*cos(Th)) - Th_2d*(x+L_piv*sin(Th));
  a = av;
  //*/
}


//&&&&& Mesh refinement  &&&&&//
/*event adapt(i++){
  adapt_wavelet ({f,u}, (double []){F_MAX,U_MAX,U_MAX},maxlevel=MAXLEVEL,minlevel=MINLEVEL);
}*/

//***** Log files  *****//
event pro_logfile(i+=100){
  fprintf(stderr, "%d %d %g %d %d \n", i, N, t, mgp.i, mgu.i);    
}

//***** Animations
/*event movies(t += dt_video; t<=t_end)
{
  scalar omega[], m[];
  vorticity (u, omega);
  
  output_ppm(omega, file = "vorticity.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
	     min = -50, max = 50, linear = true, map = cool_warm);
  //output_ppm(u.x, file = "x-vel.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //     linear = true, map = cool_warm);
  //output_ppm(u.y, file = "y-vel.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     linear = true, map = cool_warm);  
  output_ppm(f, file = "vol_fraction.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
	     min = -0.1, max = 1.1, linear = true, map = cool_warm);
  //output_ppm(ft, file = "tracer.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     linear = true, map = cool_warm);
  output_ppm(c, file = "tracer.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  	     min = 0, max = 1, linear = true, map = cool_warm);  
  output_ppm(oxy,file = "oxygen.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
	     min = 0, max = 1, linear = true, map = cool_warm);
}*/

event movies (t += 0.05){

  char timestring[100];
  
  foreach(){
	omega[] = (u.y[1,0] - u.y[-1,0])/(2.*Delta) - (u.x[0,1] - u.x[0,-1])/(2.*Delta);
        velnorm[] = sqrt(sq(u.x[]) + sq(u.y[]));
  	viewingfield[] = 1.0 - oxy[];
        viewingfieldmain[] = 1.0 - f[];
  	mylevel[] = level;
  }

  view(width=1200, height=1200, fov=24.0, ty = 0.0);
  clear();
  draw_vof("f", lw=2);
  cells(lw=0.5);
  scatter(loc, s = 30, pc = {1, 0, 0} );
  //squares("mylevel", map = cool_warm, min = MINLEVEL, max = MAXLEVEL);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Level.mp4");
  
  ///
  
  view(width=1200, height=1200, fov=24.0, ty = 0.0);
  clear();
  draw_vof("f", lw=2);
  scatter(loc, s = 30, pc = {1, 0, 0} );
  squares("omega", map = cool_warm, min = -10.0, max = 10.0);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Vorticity.mp4");

  ///
  
  /*view(width=1200, height=1200, fov=24.0, ty = 0.0);
  clear();
  draw_vof("f", lw=2);
  scatter(loc, s = 30, pc = {1, 0, 0} );
  squares("velnorm", map = cool_warm, min = 0.0, max = 1.0);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("VelField.mp4");*/

  ///
  
  view(width=1200, height=1200, fov=24.0, ty = 0.0);
  clear();
  draw_vof("f", lw=2);
  squares("viewingfield", map = cool_warm, min = 0.0, max = 2.0);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("OxygenC.mp4");

  ///

  view(width=1200, height=1200, fov=24.0, ty = 0.0);
  clear();
  draw_vof("f", lw=2);
  squares("viewingfieldmain", map = cool_warm, min = 0.0, max = 2.0);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Fluids.mp4");
}

/*
event out_files(t+=dt_file; t<=t_end_file)
{
  snprintf(buf1, sizeof(buf1), "Data_all/Data_all_%d_%.9g_%d.txt",N,t,pid());
  //snprintf(buf2, sizeof(buf2), "alpha_%g.txt",t);
  FILE * out_all = fopen(buf1,"wb");
  //FILE * interf = fopen(buf2,"w");  
    
  foreach(){
    fprintf(out_all,"%g %g %g %g %g %g %g %g\n",x,y,u.x[],u.y[],f[],c[],oxy[]);
  }
  //output_facets(f,interf);
  
  fclose(out_all);  //fclose(interf);
}
*/

event out_files_refine(t+=dt_file_refine; t<=t_end)
{
  if (t >= t_end_file){
    snprintf(buf1, sizeof(buf1), "Data_all/Data_all_%d_%.9g_%d.txt",N,t,pid());
    //snprintf(buf2, sizeof(buf2), "alpha_%g.txt",t);
    FILE * out_all = fopen(buf1,"wb");
    //FILE * interf = fopen(buf2,"w");  
    
    foreach(){
      fprintf(out_all,"%g %g %g %g %g %g %g %g\n",x,y,u.x[],u.y[],f[],c[],oxy[]);
    }
    //output_facets(f,interf);
  
    fclose(out_all);  //fclose(interf);
  }
}

event logstats (t += 0.1) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // Output i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}

event gfsview (t = 0.0; t += 10.0; t <= t_end) {
    char name_gfs[200];
    sprintf(name_gfs,"BioReactor-%0.1f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}
