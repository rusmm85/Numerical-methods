#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#define Pi  3.1415926535897932

using namespace std;

int main(void)
{
	int i, j, k, kmax=0, N=10, imax;
	double *res, *Y, *Matrix, lambda, norma=0, max=0, h=1.0/N, skal=0, skalmax=0;

	Y = new double[N+1];
	Matrix = new double[N+1];
	res = new double[N+1];
	
	for(k=1;k<N;k++)
	{
		norma=0;
		lambda=(4/(h*h))*sin(Pi*(2*k-1)/(4*N))*sin(Pi*(2*k-1)/(4*N));
		
		for(i=0;i<N+1;i++)
		{
			Y[i]=sin(Pi*i*(2*k-1)/(2*N));

		}
		
		for(i=1;i<N;i++)
		{
			Matrix[i]=(-1/(h*h))*(Y[i+1]-2*Y[i]+Y[i-1]);
		}
		
		Matrix[0]=0;
		Matrix[N]=(-2/(h*h))*(-Y[N]+Y[N-1]);
		
		for(i=0;i<N+1;i++)
		{
			res[i]=Matrix[i]-lambda*Y[i];
			norma+=res[i]*res[i];
		}
		
		norma=sqrt(norma)/fabs(lambda);

		if(norma>max)
		{
			max=norma;
			kmax=k;
		}
//cout<<"Norma = "<<norma<<endl;
	}
	
	cout<<"Norma = "<<max<<"  k = "<<kmax<<endl;
	
	for(i=0;i<N;i++)
	{
		
		
		for(k=i+1;k<N+1;k++)
		{	skal=0;
			for(j=0;j<N+1;j++)
			{
//printf ("%lf\n",sin(Pi*j*(2*i-1)/(2*N))*sin(Pi*j*(2*k-1)/(2*N)));
				
					
				if( (j==0) || (j==N))
				{
					skal+=(1.0/2.0)*sin(Pi*j*(2*i-1)/(2*N))*sin(Pi*j*(2*k-1)/(2*N));
				} else {
						skal+=sin(Pi*j*(2*i-1)/(2*N))*sin(Pi*j*(2*k-1)/(2*N));
					}
			}
		skal=skal/N;
		cout<<"(Y"<<i<<",Y"<<k<<") = "<< skal<<endl;
		
		if(skalmax<skal)
		{
			skalmax=skal;
			imax=i;
			kmax=k;
		}
		}
	}

	cout<<"skalmax = "<< skalmax<<"   imax = "<<imax<<"   kmax = "<<kmax<<endl;
	return 1;
}
