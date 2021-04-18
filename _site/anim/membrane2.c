#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double hypot3(double x,double y,double z)
{
	double value;
	value=pow(x*x+y*y+z*z,.5);
	return value;
}

double L2der(double s[], int i,int j, double m[], double r, double ro, double k,int Num)
{
	double value;
	value=(1/m[i])*((k*(r-ro)*(s[i-1+(j*Num)]-s[i+(j*Num)]))/r);
	return value;
}

double R2der(double s[], int i,int j, double m[], double r, double ro, double k,int Num)
{
	double value;
	value=(1/m[i])*((k*(r-ro)*(s[i+1+(j*Num)]-s[i+(j*Num)]))/r);
	return value;
}

double U2der(double s[], int i,int j, double m[], double r, double ro, double k,int Num)
{
	double value;
	value=(1/m[i])*((k*(r-ro)*(s[i+((j+1)*Num)]-s[i+(j*Num)]))/r);
	return value;
}

double D2der(double s[], int i,int j, double m[], double r, double ro, double k,int Num)
{
	double value;
	value=(1/m[i])*((k*(r-ro)*(s[i+((j-1)*Num)]-s[i+(j*Num)]))/r);
	return value;
}

double U(double k,double r,double ro, int i)
{
	double value;
	value=(.5)*k*(pow((r-ro),2.0));
	return value;
}

double Kin(double m[],double vy[],double vx[],double vz[],int i)
{
	double value;
	value=(.5)*m[i]*(vy[i]*vy[i]+vx[i]*vx[i]+vz[i]*vz[i]);
	return value;
}


int main(void)
{
	int i=0;
	int j=0;
	int q=0;
	int Num=20;
	int N=Num*Num;			//number of points
	double x[N];		
	double xstart[N];
	double x_half[N];
	double y[N];
	double ystart[N];
	double y_half[N];
	double z[N];
	double zstart[N];
	double z_half[N];
	double vz[N];
	double vzstart[N];
	double vz_half[N];
	double vx[N];
	double vxstart[N];
	double vx_half[N];
	double vy[N];
	double vystart[N];
int alreadydone=0;	
double vy_half[N];
	double m[N];			//mass of individual points

	double t=0;
	double dt=.0001;
	int frame_skip=1000;
	int frame=0;

	double energy=0;
	double Amp=.01;


	int xmode=1;
	int ymode=2;


	double md=1;			//mass density
	double L=1;			//Length of one side of square
	double a=1;			//stiffness
	double s=.1;			// % stretched
	double K=a/L;			//spring constant of whole
	double T=s*L*K;			//Tension
	double Ls=(T+L*K)/K;
	double dL=Ls/(Num-1);
	double k=K;
	double p=1;			//density
	double A=1;
	double M=p*A*L;

	double r=0;
	double ro=L/Num;

	double rstart=0;
	double vybefore=0;
	double vxbefore=0;
	double vzbefore=0;
	int w=(Num/(2*xmode))+(Num/(2*ymode))*Num;

	double actual_period=0;
	double past_t=0;
	

	for(i=0;i<Num;i++)
	{
	for(j=0;j<Num;j++)
	{
	x[i+(j*Num)]=i*dL;
	y[i+(j*Num)]=j*dL;

		z[i+(j*Num)]=Amp*(sin((xmode*M_PI*x[i+(j*Num)])/Ls))*(sin((ymode*M_PI*y[i+(j*Num)])/Ls));

		vx[i+(j*Num)]=0;
		vy[i+(j*Num)]=0;
		vz[i+(j*Num)]=0;

		vxstart[i+(j*Num)]=vx[i+(j*Num)];
		vx_half[i+(j*Num)]=0;

		vystart[i+(j*Num)]=vy[i+(j*Num)]; 
		vy_half[i+(j*Num)]=0;

		vzstart[i+(j*Num)]=vz[i+(j*Num)];
		vz_half[i+(j*Num)]=0;

		xstart[i+(j*Num)]=x[i+(j*Num)];
		x_half[i+(j*Num)]=x[i+(j*Num)];

		ystart[i+(j*Num)]=y[i+(j*Num)];
		y_half[i+(j*Num)]=y[i+(j*Num)];

		zstart[i+(j*Num)]=z[i+(j*Num)];
		z_half[i+(j*Num)]=z[i+(j*Num)];

		rstart=hypot3(x[i-1+(j*Num)]-x[i+(j*Num)],y[i-1+(j*Num)]-y[i+(j*Num)],z[i-1+(j*Num)]-z[i+(j*Num)]);

		m[i+(j*Num)]=(md*L*L)/(Num*Num);

	}
	}



for(t=0; 1==1; t=t+dt)
	{

		vybefore=vy[w];	
		vxbefore=vx[w];
		vzbefore=vz[w];

		for(i=1;i<Num-1;i++) //does not run for N=1
		{
		for(j=1;j<Num-1;j++)
		{
			x_half[i+(j*Num)]=x[i+(j*Num)]+vx[i+(j*Num)]*(dt/2.0);
			y_half[i+(j*Num)]=y[i+(j*Num)]+vy[i+(j*Num)]*(dt/2.0);
			z_half[i+(j*Num)]=z[i+(j*Num)]+vz[i+(j*Num)]*(dt/2.0);

			r=hypot3(x[(i-1)+(j*Num)]-x[i+(j*Num)],y[(i-1)+(j*Num)]-y[i+(j*Num)],z[(i-1)+(j*Num)]-z[i+(j*Num)]);
			vx_half[i+(j*Num)]=vx[i+(j*Num)]+L2der(x,i,j,m,r,ro,k,Num)*(dt/2.0);
			vy_half[i+(j*Num)]=vy[i+(j*Num)]+L2der(y,i,j,m,r,ro,k,Num)*(dt/2.0);
			vz_half[i+(j*Num)]=vz[i+(j*Num)]+L2der(z,i,j,m,r,ro,k,Num)*(dt/2.0);

			r=hypot3(x[(i+1)+(j*Num)]-x[i+(j*Num)],y[(i+1)+(j*Num)]-y[i+(j*Num)],z[(i+1)+(j*Num)]-z[i+(j*Num)]);
			vx_half[i+(j*Num)]=vx_half[i+(j*Num)]+R2der(x,i,j,m,r,ro,k,Num)*(dt/2.0);
			vy_half[i+(j*Num)]=vy_half[i+(j*Num)]+R2der(y,i,j,m,r,ro,k,Num)*(dt/2.0);
			vz_half[i+(j*Num)]=vz_half[i+(j*Num)]+R2der(z,i,j,m,r,ro,k,Num)*(dt/2.0);

			r=hypot3(x[i+((j+1)*Num)]-x[i+(j*Num)],y[i+((j+1)*Num)]-y[i+(j*Num)],z[i+((j+1)*Num)]-z[i+(j*Num)]);
			vx_half[i+(j*Num)]=vx_half[i+(j*Num)]+U2der(x,i,j,m,r,ro,k,Num)*(dt/2.0);
			vy_half[i+(j*Num)]=vy_half[i+(j*Num)]+U2der(y,i,j,m,r,ro,k,Num)*(dt/2.0);
			vz_half[i+(j*Num)]=vz_half[i+(j*Num)]+U2der(z,i,j,m,r,ro,k,Num)*(dt/2.0);

			r=hypot3(x[i+((j-1)*Num)]-x[i+(j*Num)],y[i+((j-1)*Num)]-y[i+(j*Num)],z[i+((j-1)*Num)]-z[i+(j*Num)]);
			vx_half[i+(j*Num)]=vx_half[i+(j*Num)]+D2der(x,i,j,m,r,ro,k,Num)*(dt/2.0);
			vy_half[i+(j*Num)]=vy_half[i+(j*Num)]+D2der(y,i,j,m,r,ro,k,Num)*(dt/2.0);
			vz_half[i+(j*Num)]=vz_half[i+(j*Num)]+D2der(z,i,j,m,r,ro,k,Num)*(dt/2.0);
		}
		}

		for(i=1;i<Num-1;i++)
		{
		for(j=1;j<Num-1;j++)
		{
			x[i+(j*Num)]=x[i+(j*Num)]+vx_half[i+(j*Num)]*dt;
			y[i+(j*Num)]=y[i+(j*Num)]+vy_half[i+(j*Num)]*dt;
			z[i+(j*Num)]=z[i+(j*Num)]+vz_half[i+(j*Num)]*dt;

			r=hypot3(x_half[(i-1)+(j*Num)]-x_half[i+(j*Num)],y_half[(i-1)+(j*Num)]-y_half[i+(j*Num)],z_half[(i-1)+(j*Num)]-z_half[i+(j*Num)]);
			vx[i+(j*Num)]=vx[i+(j*Num)]+L2der(x_half,i,j,m,r,ro,k,Num)*(dt/2.0);
			vy[i+(j*Num)]=vy[i+(j*Num)]+L2der(y_half,i,j,m,r,ro,k,Num)*(dt/2.0);
			vz[i+(j*Num)]=vz[i+(j*Num)]+L2der(z_half,i,j,m,r,ro,k,Num)*(dt/2.0);

			r=hypot3(x[(i-1)+(j*Num)]-x[i+(j*Num)],y[(i-1)+(j*Num)]-y[i+(j*Num)],z[(i-1)+(j*Num)]-z[i+(j*Num)]);
			energy=energy+U(k,r,ro,(i+j*Num));

			r=hypot3(x_half[(i+1)+(j*Num)]-x_half[i+(j*Num)],y_half[(i+1)+(j*Num)]-y_half[i+(j*Num)],z_half[(i+1)+(j*Num)]-z_half[i+(j*Num)]);
			vx[i+(j*Num)]=vx[i+(j*Num)]+R2der(x_half,i,j,m,r,ro,k,Num)*(dt/2.0);
			vy[i+(j*Num)]=vy[i+(j*Num)]+R2der(y_half,i,j,m,r,ro,k,Num)*(dt/2.0);
			vz[i+(j*Num)]=vz[i+(j*Num)]+R2der(z_half,i,j,m,r,ro,k,Num)*(dt/2.0);

			r=hypot3(x[(i+1)+(j*Num)]-x[i+(j*Num)],y[(i+1)+(j*Num)]-y[i+(j*Num)],z[(i+1)+(j*Num)]-z[i+(j*Num)]);
			energy=energy+U(k,r,ro,(i+j*Num));

			r=hypot3(x_half[i+((j+1)*Num)]-x_half[i+(j*Num)],y_half[i+((j+1)*Num)]-y_half[i+(j*Num)],z_half[i+((j+1)*Num)]-z_half[i+(j*Num)]);
			vx[i+(j*Num)]=vx[i+(j*Num)]+U2der(x_half,i,j,m,r,ro,k,Num)*(dt/2.0);
			vy[i+(j*Num)]=vy[i+(j*Num)]+U2der(y_half,i,j,m,r,ro,k,Num)*(dt/2.0);
			vz[i+(j*Num)]=vz[i+(j*Num)]+U2der(z_half,i,j,m,r,ro,k,Num)*(dt/2.0);

			r=hypot3(x[i+((j+1)*Num)]-x[i+(j*Num)],y[i+((j+1)*Num)]-y[i+(j*Num)],z[i+((j+1)*Num)]-z[i+(j*Num)]);
			energy=energy+U(k,r,ro,(i+j*Num));

			r=hypot3(x_half[i+((j-1)*Num)]-x_half[i+(j*Num)],y_half[i+((j-1)*Num)]-y_half[i+(j*Num)],z_half[i+((j-1)*Num)]-z_half[i+(j*Num)]);
			vx[i+(j*Num)]=vx[i+(j*Num)]+D2der(x_half,i,j,m,r,ro,k,Num)*(dt/2.0);
			vy[i+(j*Num)]=vy[i+(j*Num)]+D2der(y_half,i,j,m,r,ro,k,Num)*(dt/2.0);
			vz[i+(j*Num)]=vz[i+(j*Num)]+D2der(z_half,i,j,m,r,ro,k,Num)*(dt/2.0);

			r=hypot3(x[i+((j-1)*Num)]-x[i+(j*Num)],y[i+((j-1)*Num)]-y[i+(j*Num)],z[i+((j-1)*Num)]-z[i+(j*Num)]);
			energy=energy+U(k,r,ro,(i+j*Num))+Kin(m,vy,vx,vz,(i+j*Num));


			if(frame % frame_skip == 0)
			{
				for(i=0;i<Num-1;i++)
				{
				for(j=0;j<Num-1;j++)
				{
					printf("c3 %e %e %e 0.001\n",x[i+(j*Num)],y[i+(j*Num)],z[i+(j*Num)]); //draw the dots
					//printf("c3 %e %e %e 0.05\n",x[w],y[w],z[w]);	//Draws Special Circle
					printf("l3 %e %e %e %e %e %e\n",x[i+(j*Num)],y[i+(j*Num)],z[i+(j*Num)],x[(i+1)+(j*Num)],y[(i+1)+(j*Num)],z[(i+1)+(j*Num)]);		//Draws a Line to the Right
					printf("l3 %e %e %e %e %e %e\n",x[i+(j*Num)],y[i+(j*Num)],z[i+(j*Num)],x[i+((j+1)*Num)],y[i+((j+1)*Num)],z[i+((j+1)*Num)]);	//Draws a Line Upwards
				}
				}
				for(i=0;i<Num-1;i++)
				{
				j=Num-1;
				printf("l3 %e %e %e %e %e %e\n",x[i+(j*Num)],y[i+(j*Num)],z[i+(j*Num)],x[(i+1)+(j*Num)],y[(i+1)+(j*Num)],z[(i+1)+(j*Num)]);		//Draws a Line to the Right
				printf("c3 %e %e %e 0.001\n",x[i+(j*Num)],y[i+(j*Num)],z[i+(j*Num)]); //draw the dots
				printf("c3 %e %e %e 0.05\n",x[w],y[w],z[w]);
				}
				for(j=0;j<Num-1;j++)
				{
				i=Num-1;
				printf("l3 %e %e %e %e %e %e\n",x[i+(j*Num)],y[i+(j*Num)],z[i+(j*Num)],x[i+((j+1)*Num)],y[i+((j+1)*Num)],z[i+((j+1)*Num)]);	//Draws a Line Upwards
				printf("c3 %e %e %e 0.001\n",x[i+(j*Num)],y[i+(j*Num)],z[i+(j*Num)]); //draw the dots
				printf("c3 %e %e %e 0.05\n",x[w],y[w],z[w]);
				}
				printf("F\n"); // flush frame
			}
			frame++;
		
}
}
			//printf("!Total Energy %e\n",energy);
		if(vzbefore>=0)
		{
			if(vz[w]<0)
			{
				actual_period=t-past_t;
				printf("!Sign Change at %e, Period %e\n",t,actual_period);
				past_t=t;
			}
		}

			energy=0;
}
}
