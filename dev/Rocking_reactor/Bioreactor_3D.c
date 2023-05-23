#define FILTERED         // New for 3D run
#define mu(f) ( 1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)  )

#include "grid/octree.h" // New for 3D run, this needs to go first ahead of the embed (which searches for z-components)!
#include "embed.h"
#include "navier-stokes/centered.h"
//#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
// #include "/gpfs/scratch/mkim79/basilisk/henry_oxy2.h"
#include "henry_oxy2.h"
//#include "/gpfs/scratch/mkim79/basilisk/view2.h"
//#include "view2.h"
// #include "tag.h"
// #include "/gpfs/scratch/mkim79/basilisk/utils2.h"
#include "utils2.h"
#define _USE_MATH_DEFINES
#include <math.h>
// #include "view.h"

#define EMBED            1
#define CONTACT          0
#define OXYGEN           1
#define OXYGEN_CIRCLE    0
#define OXYGEN_AIR       1
#define TRACER           1
#define ACCELERATION     1
#define AMR              0
#define REMOVE_DROP      0
#define CFL_COND         0 
#define DUMP             0
#define NORMCAL          1
#define FIGURES          0
#define VIDEOS_original  0
#define VIDEOS_new       0   // this should be used with view.h
#define FIGURES_new      0
#define OUT_FILES        1
#define OUT_SPECIFIC_TIME 1
#define OUT_INTERFACE    1  // should be zero for multiple nodes (N>=1024)

//&&& Simulation setup &&&//
const double NN      = 64;     // resolution New for 3D run, lowered for testing
const double t_change= 6.8337; // change in rocking motion
const double th_cont = 90;      // contact angle
const double t_mix   = 0.3; //24.3;     // time for releasing tracers
const double t_dump  = 0.3; //24.3;     // save the dump file
const double t_end_file = 50; // t<t_end_file; saving freq=dt_file
const double t_end   = 50;    // final time
const double dt_file = 0.1519*7;   // saving frq before t_end_file
const double dt_video= 0.6074/10;  // 0.0281
const double dt_Fig  = 0.1519*7;  // same with the dt_file
const double t_spec_init = 25;  // a half of the t_end
const double t_spec_end  = 27;
const double dt_spec     = 0.01519;  // 1/10 of dt_file
const double dt_oxy  = 0.001;
const int    i_fig   = 1000;
const int    i_norm  = 5;
const double CFL_num = 0.01;    // CFL number
const double N_output= 128;    // output resolution
const double remove_minsize   = 20;
const double remove_threshold = 1.0e-4;

//&&& Reactor geometry &&&//
double LL = 1.0;        // width
double Ly = 0.286;      // 0.50; height
double Lz = 0.5;        // New for 3D run - not currently used, may be needed for confinement below
double y_init = 0.0; // -0.25*Ly; determine the filled liquid volume; y_init=0; a half volume
double L_piv  = 0.143;      // 0.143; distance to pivot point

//&&&&& Material properties &&&&&//
const double rho_w = 1.0e3;   // kg/m^3
const double rho_a = 1.225; // 0.01e3;   // kg/m^3
const double mu_w  = 1.0e-3;  // Pa-s @ 20C
const double mu_a  = 1.81e-5; // Pa-s - why multiplied 10?
const double grav  = 9.8;     // m/s^2
const double sigma = 0.0728;  // N/m
const double D_tracer_1 = 0.44e-9; // 0.44e-9; concentration = 1; diffusion in water
const double D_tracer_2 = 1.0e-30; // concentration = 0; diffusion in air
const double D_oxy_1 = 1.90e-9;    // 1.98e-5; m^2/s; phase = 1; diffusion in water
const double D_oxy_2 = 1.98e-5;    // 1.90e-9; m^2/s; phase = 0; diffusion in air
const double c_tracer_alpha = 1.0e30;    // concentration 1 to concentration 0
const double c_oxy_alpha    = 1./30;     // c_l=alpha*c_g; c1=alpha*c2; solubility; concentration 1 to concentration 0; air to water = 1/30;

//&&&&& Contact angle &&&&&//
#if CONTACT
vector h[];
h.t[left]  = contact_angle(th_cont*pi/180);
h.t[right] = contact_angle(th_cont*pi/180);
h.t[front] = contact_angle(th_cont*pi/180);  // New for 3D run
h.t[back]  = contact_angle(th_cont*pi/180);  // New for 3D run
#endif

double Th_max, T_per, R_tr, x_tr, y_tr, z_tr; // New for 3D run
scalar c[], oxy[], c1[], c2[], c3[];   // for tracer and oxygen transfer
scalar * stracers = {c,oxy,c1,c2,c3};
char buf1[100], buf2[100], buf3[100];
double (* gradient) (double, double, double) = minmod2;
double U0, Ub, Re_w, Re_a, We_w, Fr, rhor, mur, Pe_tracer_1, Pe_tracer_2, Pe_oxy_1, Pe_oxy_2, Th, Th_d, Th_2d, U_bio, w_bio, w_bio_st, T_per_st, T_bio, Th_max2, D_in_non, U_in_non;
int MINLEVEL, MAXLEVEL;
FILE * fp_stats, * fp_norm, * fp_stats2, * fp_stats3;

int main(int argc, char * argv[]){

  // Operating condition
  double L_bio = atof(argv[1]);   // m - reference length scale
  double ANGLE = atof(argv[2]);   // degree
  double RPM   = atof(argv[3]);   // RPM

  L0 = LL;                 // domian size
  origin (-L0/2., -L0/2., -L0/2.); // origin // New for 3D run
#if AMR == 0
  init_grid(NN);           // uniform grid
#endif

  // Mesh refinement
#if AMR
  init_grid(NN); // New for 3D run - we still want a decent grid to set the geometry and interface, even if it is lowered thereafter
  MINLEVEL = atof(argv[4]);
  MAXLEVEL = atof(argv[5]);
  double F_MAX = 1e-6;
  double U_MAX = 0;
#endif

  //&&& Tracer conditions &&&//
  R_tr  = 0.0084/L_bio;    // tracer radius
  x_tr  = 0.;              // tracer position
  y_tr  = -Ly*0.5*0.5;
  z_tr  = 0.;              // New for 3D run

  Th_max = ANGLE*pi/180;   // radian
  T_per  = 60./RPM;        // sec
  
  double H_bio,V_bio;
  H_bio  = L_bio*Ly;
  V_bio  = L_bio/4*(H_bio + 0.5*L_bio*tan(Th_max));
  U_bio  = V_bio/(H_bio*0.5)/T_per; // Characteristic velocity
  T_bio  = L_bio/U_bio;             // Characteristic time
  w_bio  = 2*pi/T_per;          // angular velocity rad/s
  w_bio_st = w_bio*T_bio;       // dimensionless w_bio = 2*pi
  T_per_st = T_per/T_bio;       // dimensionless rocking period
  U0     = w_bio_st*Th_max;     // initial velocity (clockwise)
  Ub     = 0.;                  // Boundary velocity

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
  
  //  Contact angle
#if CONTACT
  f.height = h;
#endif

#if TRACER  
  // D1 for c=1; D2 for c=0; c1=alpha*c2
  c.D1 = 1./Pe_tracer_1;
  c.D2 = 1./Pe_tracer_2;
  
  //c.alpha = 1.0e30;  // c1 = alpha*c2
  c.alpha = c_tracer_alpha;

  // slope limiter
  c.gradient = minmod2;

  // Other tracers: c1,c2,c3
  /*
  c1.D1 = 1./Pe_tracer_1;
  c1.D2 = 1./Pe_tracer_2;
  c1.alpha = c_tracer_alpha;
  c1.gradient = minmod2;
  c2.D1 = 1./Pe_tracer_1;
  c2.D2 = 1./Pe_tracer_2;
  c2.alpha = c_tracer_alpha;
  c2.gradient = minmod2;
  c3.D1 = 1./Pe_tracer_1;
  c3.D2 = 1./Pe_tracer_2;
  c3.alpha = c_tracer_alpha;
  c3.gradient = minmod2;
  */
#endif

#if OXYGEN
  oxy.D1 = 1./Pe_oxy_1;
  oxy.D2 = 1./Pe_oxy_2;
  oxy.alpha = c_oxy_alpha;

  // slope limiter
  oxy.gradient = minmod2;
#endif

  //&&&&& Boundary conditions &&&&&//
  u.n[left]  = dirichlet(0.);
  u.t[left]  = dirichlet(0.);
  u.n[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);
  u.n[front] = dirichlet(0.);  // New for 3D run
  u.t[front] = dirichlet(0.);  // New for 3D run
  u.n[back] = dirichlet(0.);   // New for 3D run
  u.t[back] = dirichlet(0.);   // New for 3D run


#if EMBED
  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);
#endif

#if CFL_COND
  CFL = CFL_num;
#endif
  
  char name[200],name2[200],name3[200],name4[200];

  sprintf(name, "logstats.dat");
  sprintf(name2,"normf.dat");
  sprintf(name3,"vol_frac_interf.dat");
  sprintf(name4,"tr_oxy.dat");
  fp_stats = fopen(name, "w");
  fp_norm  = fopen(name2,"w");
  fp_stats2= fopen(name3,"w");
  fp_stats3= fopen(name4,"w");

  fprintf(fp_norm, "i t Omega_liq_avg Omega_liq_rms Omega_liq_vol Omega_liq_max ux_liq_avg ux_liq_rms ux_liq_vol ux_liq_max uy_liq_avg uy_liq_rms uy_liq_vol uy_liq_max uz_liq_avg uz_liq_rms uz_liq_vol uz_liq_max \n"); // New for 3D run
  fprintf(fp_stats2, "i t f_liq_sum f_liq_interf posY_max posY_min \n");
  //fprintf(fp_stats3, "i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 c1_liq_sum c1_liq_sum2 c2_liq_sum c2_liq_sum2 c3_liq_sum c3_liq_sum2 \n");
  fprintf(fp_stats3, "i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 \n");

  // Iterations
  NITERMAX = 1000; 
  TOLERANCE = 5.0e-4;
  
  run();
  
  fclose(fp_stats); fclose(fp_norm); fclose(fp_stats2); fclose(fp_stats3);
}

//***** Initial conditions  *****//
event init (t = 0)
{
  // Initial phases: y>0 -> volume fraction=0 (air); y<0 -> volume fraction=1
  fraction (f, y_init-y);

  // solid
#if EMBED
  solid (cs,fs, intersection( -(y-0.5*Ly), -(-y-0.5*Ly) ) );
  // New for 3D: the above would need to be modified with any z-confinement should we wish to narrow the computational box
#endif
}

#if REMOVE_DROP
event remove_drop(i++){
  remove_droplets(f,remove_minsize,remove_threshold,false);  // remove droplets
  remove_droplets(f,remove_minsize,remove_threshold,true);   // remove bubbles
}
#endif

#if CFL_COND
event CFL_cond(i++){
  CFL = CFL_num;
}
#endif

#if TRACER
event tracer(t = t_mix){

  double h_tr;

  // tracer at the center of the initial liquid fill-out
  //fraction(c, -(sq(x-x_tr) + sq(y-y_tr) - sq(R_tr)) );

  // tracer at the center of the left-half side (same area)
  //fraction(c1, -(sq(x+0.25) + sq(y-y_tr) - sq(R_tr)) );

  // tracer at the center of the right-half side (same area)
  //fraction(c2, -(sq(x-0.25) + sq(y-y_tr) - sq(R_tr)) );

  h_tr = (M_PI*R_tr*R_tr);

  // tracer released as a line (same area)
  fraction(c, intersection( -(y-y_tr - 0.5*h_tr), -(-(y-y_tr + 0.5*h_tr)) ));
}
#endif

#if OXYGEN
event oxygen (t=t_mix; i++){

#if OXYGEN_CIRCLE
  fraction(oxy, -(sq(x-0) + sq(y-Ly*0.5*0.5) - sq(0.084*Ly) + sq(z-0)) );  // New for 3D run
#endif

#if OXYGEN_AIR
  foreach(){
    //if (abs(f[]) < 1e-4)
    if ((f[] == 0) && (cs[]==1))
      oxy[] = 1.;
  }
#endif
}
#endif

#if ACCELERATION
//&&&&& acceleration &&&&&&//
event acceleration(i++)
{
  // initial stage
  if (t >= t_change){
    Th    = Th_max*sin(w_bio_st*t);
    Th_d  = w_bio_st*Th_max*cos(w_bio_st*t);
    Th_2d = -w_bio_st*w_bio_st*Th_max*sin(w_bio_st*t);
  }
  // regular state
  else if (t < t_change){
    Th_max2 = (Th_max-0)/(t_change-0)*t;
    Th    = Th_max2*sin(w_bio_st*t);
    Th_d  = w_bio_st*Th_max2*cos(w_bio_st*t);
    Th_2d = -w_bio_st*w_bio_st*Th_max2*sin(w_bio_st*t);
  }
   
  face vector av = a;
  // 1st: gravitational force, 2nd: Coriolis force
  // 3rd: centrifugal force,   4th: azimuthal force, 5th: no traslational force
  foreach_face(x)
    av.x[] = -sin(Th)/(Fr*Fr) + 2*Th_d*face_value(u.y,0)	\
    + Th_d*Th_d*(x+L_piv*sin(Th)) + Th_2d*(y+L_piv*cos(Th));
  foreach_face(y)
    av.y[] = -cos(Th)/(Fr*Fr) - 2*Th_d*face_value(u.x,0)	\
    + Th_d*Th_d*(y+L_piv*cos(Th)) - Th_2d*(x+L_piv*sin(Th));
  foreach_face(z)  // New for 3D run                                  
    av.z[] = 0.0;
  a = av;
}
#endif

#if AMR
event adapt( t=0 ){  
  // Refine the domain inside bioreactor
  refine(level<MAXLEVEL && (  (y > -0.7*Ly) && (y < 0.7*Ly) ));
}
#endif

#if DUMP
event dump(t=t_dump){
  dump(file="dump");
  snprintf(buf1, sizeof(buf1), "Dump_%d_%g_%d.txt",N,t,pid());
  FILE * out_all = fopen(buf1,"w");
  foreach(){
    fprintf(out_all,"%g %g %g %g %g %g \n",x,y,u.x[],u.y[],f[],c[]);
  }
  fclose(out_all);
}
#endif

event logstats (t+=0.1; t <= t_end) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    if (pid() == 0){
      fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
      fflush(fp_stats);
    }
}

#if NORMCAL
event normcal (i+=i_norm){
  //timing s = timer_timing (perf.gt, i, perf.tnc, NULL);

    scalar ux_liq[],uy_liq[],uz_liq[],ux_liq_abs[],omega[],omega_liq[],oxy_liq[],f_liq[],posY[],c_liq[],c1_liq[],c2_liq[],c3_liq[]; // New for 3D run
    double omega_liq_avg,omega_liq_rms,omega_liq_vol,omega_liq_max;
    double ux_liq_avg,ux_liq_rms,ux_liq_vol,ux_liq_max,uy_liq_avg,uy_liq_rms,uy_liq_vol,uy_liq_max,uz_liq_avg,uz_liq_rms,uz_liq_vol,uz_liq_max; // New for 3D run
    double f_liq_sum,f_liq_interf,posY_max,posY_min,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2;
    
    vorticity (u, omega); // vorticity

    // only liquid velocity
    foreach(){
      ux_liq[]  = u.x[]*f[];
      uy_liq[]  = u.y[]*f[];
      uz_liq[]  = u.z[]*f[]; // New for 3D run
      oxy_liq[] = oxy[]*f[];
      omega_liq[] = omega[]*f[];
      f_liq[]   =  (1-cs[])*f[];
      c_liq[]   = c[]*f[];
      /*
      c1_liq[]  = c1[]*f[];
      c2_liq[]  = c2[]*f[];
      c3_liq[]  = c3[]*f[];
      */
    }

    position (f, posY, {0,1,0});  // (0,1,0) indicates the unit vector in the y-direction

    omega_liq_avg = normf(omega_liq).avg;
    omega_liq_rms = normf(omega_liq).rms;
    omega_liq_vol = normf(omega_liq).volume;
    omega_liq_max = normf(omega_liq).max;
    ux_liq_avg    = normf(ux_liq).avg;
    ux_liq_rms    = normf(ux_liq).rms;
    ux_liq_vol    = normf(ux_liq).volume;
    ux_liq_max    = normf(ux_liq).max;
    uy_liq_avg    = normf(uy_liq).avg;
    uy_liq_rms    = normf(uy_liq).rms;
    uy_liq_vol    = normf(uy_liq).volume;
    uy_liq_max    = normf(uy_liq).max;
    uz_liq_avg    = normf(uz_liq).avg;     // New for 3D run
    uz_liq_rms    = normf(uz_liq).rms;     // New for 3D run
    uz_liq_vol    = normf(uz_liq).volume;  // New for 3D run
    uz_liq_max    = normf(uz_liq).max;     // New for 3D run
    
    f_liq_sum     = statsf2(f_liq).sum;
    f_liq_interf  = interface_area(f);
    posY_max      = statsf(posY).max;
    posY_min      = statsf(posY).min;

    oxy_liq_sum   = statsf2(oxy_liq).sum;
    oxy_liq_sum2  = statsf2(oxy_liq).sum2;
    c_liq_sum     = statsf2(c_liq).sum;
    c_liq_sum2    = statsf2(c_liq).sum2;
    /*
    c1_liq_sum    = statsf(c1_liq).sum;
    c1_liq_sum2   = statsf2(c1_liq).sum2;
    c2_liq_sum    = statsf(c2_liq).sum;
    c2_liq_sum2   = statsf2(c2_liq).sum2;
    c3_liq_sum    = statsf(c3_liq).sum;
    c3_liq_sum2   = statsf2(c3_liq).sum2;
    */

   // i, timestep, no of cells, real time elapsed, cpu time
   if (pid() == 0){
      
      fprintf(fp_norm, "%i %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n",i,t,omega_liq_avg,omega_liq_rms,omega_liq_vol,omega_liq_max,ux_liq_avg,ux_liq_rms,ux_liq_vol,ux_liq_max,uy_liq_avg,uy_liq_rms,uy_liq_vol,uy_liq_max,uz_liq_avg,uz_liq_rms,uz_liq_vol,uz_liq_max); // New for 3D run
      fflush(fp_norm);

      fprintf(fp_stats2, "%i %g %g %g %g %g \n",i,t,f_liq_sum,f_liq_interf,posY_max,posY_min);
      fflush(fp_stats2);

      //fprintf(fp_stats3, "%i %g %g %g %g %g %g %g %g %g %g %g \n",i,t,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2);
      fprintf(fp_stats3, "%i %g %g %g %g %g \n",i,t,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2);
      fflush(fp_stats3);
   }
}
#endif

//event gfsview (t += 50.0) {
event gfsview (t += 1.0) {
    char name_gfs[200];
    sprintf(name_gfs,"Slices/Reactor-%1.0f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}

/*
event saveInterfaces (t += 10.0) {

    char nameInterfaces1[200];

    sprintf(nameInterfaces1,"Interfaces/interfacesLiquid-%0.1f.dat",t);

    FILE * fp1 = fopen(nameInterfaces1, "w");
    output_facets (f, fp1);	
    fclose(fp1);
}
*/

#if VIDEOS_original 
//***** Animations
event movies(t += dt_video; t<=t_end)
{
  scalar omega[], cc[], oxyy[];
  vorticity (u, omega);
  
  //output_ppm(omega, file = "vorticity.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
	//    min = -50, max = 50, linear = true, map = cool_warm);
  //output_ppm(u.x, file = "x-vel.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     min = -0.1, max = 0.1, linear = true, map = cool_warm);
  //output_ppm(u.y, file = "y-vel.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //  	     min = -0.1, max = 0.1, linear = true, map = cool_warm);  
  output_ppm(f, file = "vol_fraction.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
	     min = 0.0, max = 1.0, linear = true, map = cool_warm);
  //output_ppm(ft, file = "tracer.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     linear = true, map = cool_warm);
  //output_ppm(c, file = "tracer.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     min = 0, max = 1, linear = true, map = cool_warm);  
  //output_ppm(c1, file = "tracer1.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     min = 0, max = 1, linear = true, map = cool_warm);  
  //output_ppm(c2, file = "tracer2.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     min = 0, max = 1, linear = true, map = cool_warm);  
  //output_ppm(c3, file = "tracer3.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     min = 0, max = 1, linear = true, map = cool_warm);  
  //output_ppm(oxyy,file = "oxygen.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
//	     min = 0, max = 0.04, linear = true, map = cool_warm);
}
#endif

#if VIDEOS_new
event movies_upgrade(t += dt_video; t<=t_end)
{
  scalar omega[], ff[], cc[], oxyy[];
  char timestring[100];

  vorticity (u,omega);

  // vorticity
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("omega",map=cool_warm,min=-50.0,max=50.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("vorticity3.mp4");
  save("vorticity.png");

  // volume fraction
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("f",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  //cells();
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("volume_fraction3.mp4");
  save("volume_fraction.png");

  // tracer
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c",map=cool_warm,min=0.0,max=0.1);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("tracer.mp4");
  save("tracer.png");
/*
  // tracer1
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c1",map=cool_warm,min=0.0,max=0.1);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("tracer1.mp4");
  save("tracer1.png");

  // tracer2
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c2",map=cool_warm,min=0.0,max=0.1);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("tracer2.mp4");
  save("tracer2.png");

  // tracer3
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c3",map=cool_warm,min=0.0,max=0.1);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("tracer3.mp4");
  save("tracer3.png");
*/
  // oxygen
  #if OXYGEN
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("oxy",map=cool_warm,min=0.0,max=0.04);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("oxygen3.mp4");
  save("oxygen.png");
  #endif
}
#endif

#if FIGURES_new
event Figures_new(t=t_mix; t<=t_end; t += dt_Fig)
{
  scalar omega[];
  char timestring[100],figN1[100],figN2[100],figN3[100],figN4[100];
  
  vorticity (u,omega);
  
  snprintf(figN1, sizeof(figN1), "Fig_vor/vor_%d_%.9g.png",N,t);
  snprintf(figN2, sizeof(figN2), "Fig_vol/vor_%d_%.9g.png",N,t);
  snprintf(figN3, sizeof(figN3), "Fig_tr/vor_%d_%.9g.png",N,t);
  snprintf(figN4, sizeof(figN4), "Fig_oxy/vor_%d_%.9g.png",N,t);

  // vorticity
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("omega",map=cool_warm,min=-50.0,max=50.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN1);

  // volume fraction
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("f",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  //cells();
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN2);

  // tracer
  #if TRACER
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c",map=cool_warm,min=0.0,max=0.1);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN3);
  #endif
  
  // oxygen
  #if OXYGEN
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("oxy",map=cool_warm,min=0.0,max=0.04);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN4);
  #endif
}
#endif

#if FIGURES
event figures(i += i_fig; t <= t_end)
{
  scalar omega[];
  vorticity (u,omega);
  
  output_ppm(f, file = "vol_frac.png");
  #if TRACER
  output_ppm(c, file = "tracer.png");
  /*
  output_ppm(c1, file = "tracer1.png");
  output_ppm(c2, file = "tracer2.png");
  output_ppm(c3, file = "tracer3.png");
  */
  #endif
  //output_ppm(omega, file="vorticity.png");
  //output_ppm(u.x, file="x-vel.png");
  //output_ppm(u.y, file="y-vel.png");
  //output_ppm( cs, file="solid.png");
  #if OXYGEN
  output_ppm(oxy, file="oxygen.png");
  #endif
}
#endif

#if OUT_FILES
event out_files(t+=dt_file; t<=t_end_file)
{
  scalar omega[];
  vorticity(u,omega);

  ///*
  snprintf(buf1, sizeof(buf1), "Data_all/Data_all_%d_%.9g_%d.txt",N,t,pid());
  FILE * out_all = fopen(buf1,"wb");  

  fprintf(out_all,"x y z ux uy uz vol_frac tracer solid oxygen vorticity \n");
  foreach()
    //fprintf(out_all,"%g %g %g %g %g %g %g %g %g %g %g\n",x,y,u.x[],u.y[],f[],c[],cs[],oxy[],c1[],c2[],c3[]);
    fprintf(out_all,"%g %g %g %g %g %g %g %g %g %g \n",x,y,z,u.x[],u.y[],u.z[],f[],c[],cs[],oxy[],omega[]); // New for 3D run
  fclose(out_all);
  //*/
/*
  char name[80];
  sprintf(name,"Data_all/Output_%d_%.4g.dat",N,t);
  FILE * out_all2 = fopen(name,"w");
  output_field({u.x,u.y,f,c,cs,oxy},out_all2,n=N_output,linear=true);
  fclose(out_all2);
*/
}
#endif

#if OUT_SPECIFIC_TIME
event out_spec_time(t=t_spec_init; t<=t_spec_end; t+=dt_spec)
{
  scalar omega[];
  vorticity(u,omega);

  snprintf(buf2, sizeof(buf2), "Data_specific/Data_all_%d_%.9g_%d.txt",N,t,pid());
  FILE * out_all_spec = fopen(buf2,"wb");  

  fprintf(out_all_spec,"x y z ux uy uz vol_frac tracer solid oxygen vorticity \n");
  foreach()
    fprintf(out_all_spec,"%g %g %g %g %g %g %g %g %g %g %g \n",x,y,z,u.x[],u.y[],u.z[],f[],c[],cs[],oxy[],omega[]); // New for 3D run
  fclose(out_all_spec);

  // Only works for single node; for multiiple nodes, this should be turned off
  #if OUT_INTERFACE
    snprintf(buf3, sizeof(buf3), "Data_specific/Interf_%d_%.9g_%d.txt",N,t,pid());
    FILE * out_interf_spec = fopen(buf3,"wb");
    output_facets(f,out_interf_spec);   // Interface extraction
    fclose(out_interf_spec);
  #endif
}
#endif
