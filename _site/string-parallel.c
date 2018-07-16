#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// #define ONE_PERIOD_ONLY // uncomment this to do only one period, then print stats and exit. 
void get_forces(double *x, double *y, double *Fx, double *Fy, int N, double k, double r0) // pointers here are array parameters
{
  int i;
  double r;
  Fx[0]=Fy[0]=Fx[N]=Fy[N]=0; // have to set the end forces properly to avoid possible uninitialized memory shenanigans
  double U=0,E=0,T=0;
 
  for (i=1;i<N;i++)
  {
    // left force
    r = hypot(x[i]-x[i-1],y[i]-y[i-1]);
    Fx[i] = -(x[i]-x[i-1]) * k * (r-r0)/r;
    Fy[i] = -(y[i]-y[i-1]) * k * (r-r0)/r;

    // right force
    r = hypot(x[i]-x[i+1],y[i]-y[i+1]);
    Fx[i] += -(x[i]-x[i+1]) * k * (r-r0)/r;
    Fy[i] += -(y[i]-y[i+1]) * k * (r-r0)/r;
  }
}

 void evolve_leapfrog(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double dt) 
{
  int i;


  double Fx[N+1],Fy[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
                          // to avoid having to deal with malloc(). In any case memory allocation is faster than a bunch
                          // of square root calls in hypot().
  for (i=1;i<N;i++)
  {
    x[i] += vx[i]*dt/2;
    y[i] += vy[i]*dt/2;
  }
  get_forces(x,y,Fx,Fy,N,k,r0);
  for (i=1;i<N;i++)
  {
    vx[i] += Fx[i]/m*dt;
    vy[i] += Fy[i]/m*dt;
  }
  for (i=1;i<N;i++)
  {
    x[i] += vx[i]*dt/2;
    y[i] += vy[i]*dt/2;
  }
}


void evolve_euler(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double dt) 
{
  int i;


  double Fx[N+1],Fy[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
                          // to avoid having to deal with malloc(). In any case memory allocation is faster than a bunch
                          // of square root calls in hypot().
  get_forces(x,y,Fx,Fy,N,k,r0);
  for (i=1;i<N;i++)
  {
    x[i] += vx[i]*dt;
    y[i] += vy[i]*dt;
    vx[i] += Fx[i]/m*dt;
    vy[i] += Fy[i]/m*dt;
  }
}

// this function is around to go from Euler-Cromer to leapfrog, if we want second-order precision
void evolve_velocity_half(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double dt)
{
  int i;

  double Fx[N+1],Fy[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
  
  get_forces(x,y,Fx,Fy,N,k,r0);
  
  for (i=1;i<N;i++)
  {
    vx[i] += Fx[i]/m*dt/2;
    vy[i] += Fy[i]/m*dt/2;
  }
}

// Students might not be familiar with pass-by-reference as a trick for returning multiple values yet. 
// Ideally they should be coding this anyway, and there are a number of workarounds, in particular 
// just not using a function for this.
void get_energy(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double *E, double *T, double *U)
{
  *E=*T=*U=0;
  int i;
  double r;
  for (i=0;i<N;i++)
  {
    *T+=0.5*m*(vx[i]*vx[i] + vy[i]*vy[i]);
    r = hypot(x[i]-x[i+1],y[i]-y[i+1]);
    *U+=0.5*k*(r-r0)*(r-r0);
  }
  *E=*T+*U;
}

// does what it says on the tin
void evolve_euler_cromer(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double dt)
{
  int i;
  double Fx[N+1],Fy[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
  for (i=1;i<N;i++)
  {
    x[i] += vx[i]*dt;
    y[i] += vy[i]*dt;
  }
  
  get_forces(x,y,Fx,Fy,N,k,r0);
  
  for (i=1;i<N;i++)
  {
    vx[i] += Fx[i]/m*dt;
    vy[i] += Fy[i]/m*dt;
  }
}

// does what it says on the tin
void evolve_rk2(double *x, double *y, double *vx, double *vy, int N, double k, double m, double r0, double dt)
{
  int i;

  double Fx[N+1],Fy[N+1]; // this could be made faster by mallocing this once, but we leave it this way for students
  double xh[N+1],yh[N+1],vxh[N+1],vyh[N+1];
  vxh[0]=vyh[0]=vxh[N]=vyh[N]=0;
  get_forces(x,y,Fx,Fy,N,k,r0);

  for (i=0;i<=N;i++)
  {
    xh[i] = x[i] + vx[i]*dt/2;
    yh[i] = y[i] + vy[i]*dt/2;
    vxh[i] = vx[i] + Fx[i]/m*dt/2;
    vyh[i] = vy[i] + Fy[i]/m*dt/2;
  }
  
  get_forces(xh,yh,Fx,Fy,N,k,r0);
  
  for (i=0;i<=N;i++) // need two for loops -- can't interleave halfstep/fullstep updates (students have trouble with this sometimes!)
  {
    x[i] = x[i] + vx[i]*dt;
    y[i] = y[i] + vy[i]*dt;
    vx[i] = vx[i] + Fx[i]/m*dt;
    vy[i] = vy[i] + Fy[i]/m*dt;
  }
  
}

// function to encapsulate determining whether we need to shovel another frame to the animator. delay is the delay in msec.
// I usually aim to hide library calls, like clock() and CLOCKS_PER_SEC, from students in the beginning, since it's not really
// relevant to their development of computational skills, which are what I really care about.
int istime(int delay)
{
  static int nextdraw=0;
  if (clock() > nextdraw)
  {
    nextdraw = clock() + delay * CLOCKS_PER_SEC/1000.0;
    return 1;
  }
  return 0;
}

int main(int argc, char **argv)
{
  int i,N=80; //  number of links, not number of nodes!! Careful for the lurking fencepost errors
  int modenumber=2; // put in some defaults just in case 
  double t, dt=6e-5;
  double stiffness=10, density=1, length=1; // unstretched properties of original string
  double k, m, r0; // properties of single string
  double tension=1,Ls;
  int frame=0, frameskip;
  int Np=20;
  double amp_first=1e-2;
  double amp_last=0.2;
  if (argc < 12) // if they've not given me the parameters I need, don't just segfault -- tell the user what to do, then let them try again
  {
    fprintf(stderr,"Usage: <this> <N> <modenumber> <dt> <stiffness> <density> <length> <tension> <N_par> <amp_first> <amp_last>\nRunning with some defaults for now...\n"); 
  }
  else
  {
    N=atoi(argv[1]);
    modenumber=atoi(argv[2]);
    dt=atof(argv[3]);
    stiffness=atof(argv[4]);
    density=atof(argv[5]);
    length=atof(argv[6]);
    tension=atof(argv[7]);
    Np=atoi(argv[8]);
    amp_first=atof(argv[9]);
    amp_last=atof(argv[10]);
  }
  double amp_step=pow(amp_last/amp_first,1.0/(Np-1));
  double amplitude[Np];
  for (int p=0; p<Np; p++) amplitude[p]=amp_first * pow(amp_step,p);
  double pstep=2.0/Np;

  
  double x[(N+1)*Np], y[(N+1)*Np], vx[(N+1)*Np], vy[(N+1)*Np], E[Np], T[Np], U[Np];

  // compute microscopic properties from macroscopic ones 

  r0=length/N;
  m=density*length/N;
  k=stiffness*N/length;

  // figure out stretched length

  Ls=length + tension * length / stiffness;

  // make predictions based on what our freshman mechanics class taught us

  double density_stretched = density * length / Ls;
  double wavespeed = sqrt(tension/density_stretched);
  double period_predict = 2 * Ls / wavespeed / modenumber;
  double vym_last[Np]; for (int p=0;p<Np;p++) vym_last[p]=0;

  int monitor_node = N/modenumber/2; // this is the node that we'll be watching to see when a period has elapsed.
  int nperiods[Np]; for (int p=0;p<Np;p++) nperiods[p]=0;

  for (int p=0;p<Np;p++)
  {
//    printf("!Initializing string %d with amplitude %.3e\n",p,amplitude[p]);
    for (i=0;i<=N;i++) // remember, we have N+1 of these
    {
      x[i+(N+1)*p] = Ls*i/N - Ls/2;
      y[i+(N+1)*p] = amplitude[p]*sin(modenumber * M_PI * (x[i+(N+1)*p]+Ls/2) / Ls);
      vx[i+(N+1)*p]=0;
      vy[i+(N+1)*p]=0;
    }
  }

  double velampnow[Np], velamplast[Np];
  for (int p=0; p<Np; p++) velampnow[p]=velamplast[p]=0;
  // now, loop over time forever...
  for (t=0;1;t+=dt)
  {
    for (int p=0; p<Np; p++)
    {
      evolve_leapfrog(x+(N+1)*p,y+(N+1)*p,vx+(N+1)*p,vy+(N+1)*p,N,k,m,r0,dt);  // pointer hacking to point at the right parallel spring
      
      velamplast[p]=velampnow[p];
      velampnow[p]=0;
      for (i=0; i<N; i++)
      {
	velampnow[p] += vy[i+(N+1)*p]*sin(modenumber * M_PI * (x[i+(N+1)*p]+Ls/2) / Ls);
      }


	

    // "if we were going up, but now we're going down, then a period is complete"
    if (velamplast[p] > 0 && velampnow[p] < 0)
    {
      double tf=t + dt + velampnow[p] / (velamplast[p] - velampnow[p]) * dt;
      // now do extrapolation since we overshot
      nperiods[p]++;
      if (nperiods[p]==1) printf("!String with amplitude %e: measured period %e, predicted period %e, fractional difference %e\n",amplitude[p],tf,period_predict,fabs(1-tf/period_predict));
    }
    }
    frame++;
    if (istime(30))
    {
      for (int p=0; p<Np; p++)
      {
	for (i=0;i<=N;i++)
	{
          double mult;
	  printf("C %e 0.5 %e\n",0.5+y[i+(N+1)*p]/amplitude[p],0.5-y[i+(N+1)*p]/amplitude[p]); // use red/blue shading; this will make vibrations visible even if amp<<1
	  printf("c3 %f %f %f %f\n",x[i+(N+1)*p],(float)p*pstep-Np/2*pstep,y[i+(N+1)*p],length/N/2); // draw circles, with radius scaled to separation
	  if (i<N) printf("l3 %f %f %f %f %f %f\n",x[i+(N+1)*p],(float)(p)*pstep-Np/2*pstep,y[i+(N+1)*p],x[i+1+(N+1)*p],(float)p*pstep-Np/2*pstep,y[i+1+(N+1)*p]); // the if call ensures we don't drive off the array
	}
      }
      printf("T 0 -0.9\nTimesteps = %d\n",frame);
      printf("F\n"); // flush frame
    }
  }
}
