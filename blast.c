//    Code written in C (Basilisk.fr notation) to solve a problem using Basilisk code...

#include "grid/multigrid1D.h"
#include "green-naghdi-blast.h"
//#include "saint-venant-blast.h"


double h0 = 3000.;
//double slope = 1./19.85;


int main() {
  X0 = 0.;
  L0 = 400.*1000.;
  N = 1024*8;
  G = 9.81;
  //alpha_d = 1.;
  
  run();
}


double burst_p (double dist, double time)
{
  double thick = 17.0;
  double speed = 0.3915;
  double width = 30.0;

  double maxAmp1 = 96.0;
  double maxAmp2 = 96.0;
  double maxAmp3 = 5.0;
  double maxAmp4 = 34.0;
  double c = width/2.35482;
  double g;
  double currentRadius = time*speed;
  if (currentRadius > dist){
    g = 0. + maxAmp1* exp(-4.0*pow((currentRadius/c       ),2.0));
    g = g  + maxAmp2* exp(-2.0*pow((currentRadius/(1.80*c)),2.0));
    g = g  + maxAmp3* exp(-2.0*pow((currentRadius/(5.00*c)),8.0));
    g = g  + maxAmp4* exp(-2.0*   ((currentRadius/(12.0*c))    ));

    return g*exp(     -3.80*(currentRadius - dist)/(2.0*thick))
               *(1.0 - 6.00*(currentRadius - dist)/(2.0*thick));
  }
  else{
    return 0.0;
  }
}


p_burst[left] = 0.5*burst_p(0., t);

event calc_p (i++)
{
  foreach() {
    p_burst[] = 0.5*burst_p(x/1000., t);
  }
  boundary ({p_burst})
}






event init (i = 0)
{
  foreach() {
    zb[] = -h0; //x<0.0 ? 0 : 1.;
    h[] = -zb[]; //(hr+aa/2)-aa/2.*tanh((x-xd)/0.001); //x<0.0 ? 4.0 : 1.0299; //x<750 ? 20-zb[]: 15-zb[];
    u.x[] = 0.0;
  }
}


/*
event gnuplot (i += 5) {
  static FILE * fp = popen ("gnuplot", "w");

  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "p [0:300.*1000.][-12:5]'-' u 1:3:2 w filledcu lc 3 t ''\n", t);

  foreach()
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n\n");
  fprintf (stderr, "%.3f %.3f\n", t, statsf(u.x).max);
}
*/

event logfile (i++) {
  fprintf (stderr, "%g\n", t);
  static FILE * fp1 = fopen ("g1", "w");
  fprintf (fp1, "%g %g\n", t, interpolate (eta, 100.*1000, 0));
}

event output (t <= 1600.; t += 10.) {

  int NPP = 1024*8;
  for (int i = 0; i < NPP; i++) {
      double xz = X0 + (L0/NPP)*(i+1);
      fprintf (stdout, "%g %g\n",  xz, interpolate(eta,xz));
  }
//  fprintf (stdout, "\n");

}


