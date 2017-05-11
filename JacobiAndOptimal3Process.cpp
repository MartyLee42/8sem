#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Pi 2*asin(1)
#define EPS 1E-10
#define a1 0.05
#define b1 10

double F(double x, double y);
void Jakobi(double *y, double *f, double *x,int n);
void IterMeth(double *f, double *x,double *xx, double *xxx, int n,double alpha);
void fill(double *y, double *f, int n,double alpha);
void output(double *y, int n,char *);
double max(double x, double y);
void print(double *y, int n);
double norm(double *y, int n);
double F2(double x, double y);

int main()
{
	double *y,*f,*x,*xx,*xxx;
	int n;
	double alpha;
	int i;

	scanf("%d",&n);
	scanf("%lf",&alpha);

	y=(double *)malloc(sizeof(double)*(n+1)*(n+1));
	f=(double *)malloc(sizeof(double)*(n+1)*(n+1));
	x=(double *)malloc(sizeof(double)*(n+1)*(n+1));
	xx=(double *)malloc(sizeof(double)*(n+1)*(n+1));
	xxx=(double *)malloc(sizeof(double)*(n+1)*(n+1));

	fill(y,f,n,alpha);
	Jakobi(y,f,x,n);
	output(y,n,"out1.txt");	

	printf("First norm = %lf\n", norm(y,n));

	fill(y,f,n,alpha);
	IterMeth(f,y,xx,xxx,n,alpha);	
	output(y,n,"out2.txt");

	printf("Second norm = %lf\n", norm(y,n));

}

void IterMeth(double *f, double *x, double *xx, double *xxx, int n,double alpha)
{
	double tau;
	double ksi;
	double g0;
	double g1;
	double tmp;

	double gamma0;
	

	int k;

	int i,j;

	gamma0=(b1-a1)/(b1+a1);

	g0=gamma0;
	g1=0;
	

	for(i=0; i<=n; i++)
	{
		for(j=0; j<=n; j++)
		{
			xxx[i*(n+1)+j]=alpha;
			xx[i*(n+1)+j]=alpha;
		}
	}

	k=0;
	while(k<n*n+10)
	{
		k++;
		for(i=1; i<n; i++)
		{
			x[i*(n+1)+1] = xx[i*(n+1)+1] + g0*g1*(xx[i*(n+1)+1] - xxx[i*(n+1)+1]) - (2.0/(a1+b1))*(1 + g0*g1)*( 0-xx[(i-1)*(n+1)+1] -xx[(i+1)*(n+1)+1]  - xx[i*(n+1)+2] +4*xx[i*(n+1)+1] - f[i*(n+1)+1]);
			for(j=2; j<n-1; j++)
			{
				x[i*(n+1)+j] = xx[i*(n+1)+j] + g0*g1*(xx[i*(n+1)+j] - xxx[i*(n+1)+j]) - (2.0/(a1+b1))*(1 + g0*g1)*( 0-xx[(i-1)*(n+1)+j] -xx[(i+1)*(n+1)+j] - xx[i*(n+1)+j-1] - xx[i*(n+1)+j+1] +4*xx[i*(n+1)+j] - f[i*(n+1)+j]);
			}
			x[i*(n+1)+n-1] = xx[i*(n+1)+n-1] + g0*g1*(xx[i*(n+1)+n-1] - xxx[i*(n+1)+n-1]) - (2.0/(a1+b1))*(1 + g0*g1)*( 0-xx[(i-1)*(n+1)+n-1] -xx[(i+1)*(n+1)+n-1] - xx[i*(n+1)+n-2] +4*xx[i*(n+1)+n-1] - f[i*(n+1)+n-1]);
		}


		
		tmp=g1;
		g1=g0;
		g0=1.0/(2.0/gamma0 - tmp);

		for(i=1; i<n; i++)
		{
			for(j=1; j<n; j++)
			{
				xxx[i*(n+1)+j]=xx[i*(n+1)+j];
				xx[i*(n+1)+j]=x[i*(n+1)+j];
			}
		}
		
		
	}

	return;
}


void Jakobi(double *y, double *f, double *x,int n)
{
	double h;
	int i,j;
	double norm;

	while(true)
	{
		norm=0;		

		for(i=1; i<n; i++)
		{

			x[i*(n+1)+1] = (f[i*(n+1)+1] + y[(i-1)*(n+1)+1] + y[(i+1)*(n+1)+1]  + y[i*(n+1)+2])/4.0;
			
			for(j=2; j<n-1; j++)
			{
				x[i*(n+1)+j] = (f[i*(n+1)+j] + y[(i-1)*(n+1)+j] + y[(i+1)*(n+1)+j]  + y[i*(n+1)+j-1]+y[i*(n+1)+j+1])/4.0;				
			}

			x[i*(n+1)+n-1] = (f[i*(n+1)+n-1] + y[(i-1)*(n+1)+n-1] + y[(i+1)*(n+1)+n-1]  + y[i*(n+1)+n-2])/4.0;
		}		

		for(i=1; i<n; i++)
		{
			for(j=1; j<n; j++)
			{
				norm=max(norm,fabs(y[i*(n+1)+j]-x[i*(n+1)+j]));
				y[i*(n+1)+j]=x[i*(n+1)+j];
			}
		}

		if(norm<EPS) break;

	}
	
	return;
}

double max(double x, double y)
{
	if(x>y) return x;
	return y;
}

void print(double *y, int n)
{
	for(int i=0; i<n; i++)
	{
		printf("%lf ",y[i]);
	}
	printf("\n");
}

void fill(double *y, double *f, int n,double alpha)
{
	int i,j;
	double h;

	h=(1.0)/n;

	for(i=0;i<=n; i++)
	{
		for(j=0; j<=n; j++)
		{
			y[i*(n+1)+j] = alpha;
			
		}
	}
	
	for(i=1;i<n; i++)
	{
		f[i*(n+1)+1]=h*h*F(i*h,h)+alpha;
		for(j=2; j<n-1; j++)
		{
			f[i*(n+1)+j]=h*h*F(i*h,j*h);
			
		}
		f[i*(n+1)+n-1]=h*h*F(i*h,(n-1)*h)+alpha;
	}

	return;
}

void output(double *y, int n,char *a)
{
	int i;
	int j;
	FILE *file;

	file=fopen(a,"w");

	for(i=0; i<=n; i++)
	{
		for(j=0; j<=n; j++)
		{
			fprintf(file,"%lf ",y[i*(n+1)+j]);
		}
		fprintf(file, "\n");
	}

	fclose(file);
}

double norm(double *y, int n)
{
	int i,j;
	double ans;
	double h;

	h=1.0/(n);
	ans=0.0;
	for(i=0; i<=n; i++)
	{
		for(j=0;j<=n; j++)
		{
			ans=max(ans,fabs(F2(i*h,j*h)-y[i*(n+1)+j]));
		}
	}
	return ans;
}

double F2(double x, double y)
{
	return sin(Pi*x)*sin(Pi*y);
}

double F(double x, double y)
{
	return 2*Pi*Pi*sin(Pi*x)*sin(Pi*y);
}
