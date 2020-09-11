#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#define N 10
#define Pi 3.1415926535897932384626433832795
using namespace std;

double scalar (double *a, double *b)
{
	int i;
	double scal=0;
	double h=1.0/N;

	for(i=1;i<N;i++)
	{
		scal+=h*a[i]*b[i];
	}

	return scal;
}


void Phi (double *phi, int n)
{
	int i;
	double h=1.0/N;

	for(i=1;i<N;i++)
	{
		phi[i]=sqrt(2.0)*sin(Pi*n*i*h);
	}
	
	return;
}


double CM (double *F, double *phi, int m)
{
	double Cm, lambda, h=1.0/N;	
	Phi(phi,m);
	lambda=(4/(h*h))*sin(Pi*m*h/2)*sin(Pi*m*h/2);

	Cm=scalar(F,phi)/lambda;

	return Cm;
}

int main()
{
	int i, j;
	double h=1.0/N, Cm;
	double *phi;
	double *Y;
	double *Y1;
	double *F;
	
	phi= new double[N];
	Y= new double[N+1];
	Y1= new double[N];
	F= new double[N];

	for(i=1;i<N;i++)
	{
		Y1[i]=0;
		Y[i]=rand()%10;
	}

	Y[0]=Y[N]=0;

	for(i=1;i<N;i++)
	{
		F[i]=(-Y[i+1]+2*Y[i]-Y[i-1])/(h*h);
	}
	
	F[0]=0;

	for(i=1;i<N;i++)
	{
		Cm = CM(F, phi, i);
		//Y1[i]=0;
		for(j=1;j<N;j++)
		{
			Y1[j]=Y1[j]+Cm*phi[j];
		}
	}

	for(i=1;i<N;i++)
	{
		cout<<"i = "<<i<<" Y = "<<Y[i]<<" Y1= "<<Y1[i]<<endl;
	}
	
	return 1;
} 

	



