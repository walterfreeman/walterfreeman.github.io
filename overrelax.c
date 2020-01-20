#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void col(double color)
{
  printf("C %f %f %f\n",8+color,-color/2,-color/6);
}

main()
{
  int NFPS;
  int ws;
  double color;
  double eps,epsref;
  int N,i,j;
  double T[200][200],epshere[200][200];
  double A;
  double amp=0.3;
  double a;
  int step,delay,FPS;
  int waitsteps;
  N=80;
  a=(double)1/(double)N;
  A=1.90;
  waitsteps=0;
  FPS=60;
  NFPS=60;
  double x,y,epsmax;
  double navg;
  for (i=0;i<N;i++)
  {
    x=(double)(i)/(double)N;
    for (j=0;j<N;j++) 
    {
      y=(double)j/(double)N;
      T[i][j]=sin(x*M_PI*12)*sin(y*M_PI*9)*amp-drand48()*amp*1.9;
//      T[i][j]=0;
      if ((i-N/2)*(i-N/2)+(j-N/2)*(j-N/2) > N*N/16) T[i][j]=0;
      if (i==0.2*N && j>N/4 && j<3*N/4) T[i][j]=-0.3;
      if (i>0.4*N && i<0.9*N && j == 0.5*N) T[i][j]=0.3;

    }
  }

  printf("FPS %d\n",NFPS); 
  for (i=0;i<=N;i++)
    for (j=0;j<=N;j++)

  for (step=0;step<1000000;step++)
  {
    eps=0;epsmax=0;
    for (i=1;i<N;i++)
    {
      for (j=1;j<N;j++)
      {
      if (i==0.2*N && j>N/4 && j<3*N/4) continue;
      if (i>0.4*N && i<0.9*N && j == 0.5*N) continue;
        navg=0.25 * (T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1]);
        epshere[i][j] = fabs(T[i][j]-navg)/N/N;
        if (epshere[i][j] > epsmax) epsmax=epshere[i][j];
        eps += epshere[i][j];
        T[i][j] = T[i][j] - (T[i][j]-navg)*A;
      }
    }
    if (step == 0) epsref=epsmax;
    if (epsmax / epsref < 0.1) epsref=epsmax;
    if (step == 0) ws = NFPS*waitsteps; else ws=NFPS/FPS;
    for (delay=0;delay<ws;delay++)
    {

    for (x=epsref;x>epsref/1e6;x*=pow(10,-.3))
    {
      color=log(x/epsref);
      col(color);
//      printf("l %f %f %f %f\n",0.9,0.9+color/8,0.9,0.9+(color-0.7)/8);
//      printf("t %f %f\n%.0e\n",0.92,0.9+color/8,x);
    }

    for (i=0;i<N;i++)
    {
      col(-8);
      x=(double)i/(double)N-0.5;
      printf("l3 %f %f %f %f %f %f\n", x, 0.5, T[i][N],x+a, 0.5,T[i+1][N]);
      printf("l3 %f %f %f %f %f %f\n", 0.5, x, T[N][i],0.5, x+a,T[N][i+1]);
      for (j=0;j<N;j++)
      {
        y=(double)j/(double)N-0.5;
        color=log(0.5*(epshere[i][j]+epshere[i+1][j]))-log(epsref);
        if (j==0) color=-8;
        if (i==0) color=-8;
        col(color);  
        printf("q3 %e %e %e %e %e %e %e %e %e %e %e %e\n",x,y,T[i][j],
                                                          x+a,y,T[i+1][j],
                                                          x+a,y+a,T[i+1][j+1],
                                                          x,y+a,T[i][j+1]);
      }
    }
    printf("T -0.2 -0.9\n");
    printf("step %d\tepsilon %le\tepsmax %le\n",step,eps,epsmax);
    printf("F\n");
   }
}
  
  for (i=0;i<N;i++)
  for (j=0;j<N;j++) fprintf(stderr,"%d %d %f\n",i,j,T[i][j]);
  double sum=0;
  for (i=0;i<N;i++)
  {
    sum += (T[0][i] - T[1][i])*4; // units deg
  }
} 
