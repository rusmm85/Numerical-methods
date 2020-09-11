#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#define PI 3.1415926535897932384626433832795

using namespace std;

double func(double x)
{
	return sin(x);
}

double power (double x, double a){
	for (int i=0;i<a;i++)
	x*=x;
	return x;
}
double lagr(double*values, double*points, double x, int N)
{
	double y=0;
	double p=1;

	for(int i=0;i<N;i++)
	{
		p=1;
		for(int j=0;j<N;j++)
		{
			if(j!=i)
			{
				p*=(x-points[j])/(points[i]-points[j]);
			}

		}
		y+=values[i]*p;
	}
	return y;
}

double * gauss(double **a, double *y, int n)
{
  double *x, max;
  int k, index;
  const double eps = 0.00001;  
  x = new double[n];
  k = 0;
  while (k < n)
  {
    
    max = abs(a[k][k]);
    index = k;
    for (int i = k + 1; i < n; i++)
    {
      if (abs(a[i][k]) > max)
      {
        max = abs(a[i][k]);
        index = i;
      }
    }
   
    if (max < eps)
    {
     
      cout << "Ðåøåíèå ïîëó÷èòü íåâîçìîæíî èç-çà íóëåâîãî ñòîëáöà ";
      cout << index << " ìàòðèöû A" << endl;
      return 0;
    }
    for (int j = 0; j < n; j++)
    {
      double temp = a[k][j];
      a[k][j] = a[index][j];
      a[index][j] = temp;
    }
    double temp = y[k];
    y[k] = y[index];
    y[index] = temp;
  
    for (int i = k; i < n; i++)
    {
      double temp = a[i][k];
      if (abs(temp) < eps) continue; 
      for (int j = 0; j < n; j++)
        a[i][j] = a[i][j] / temp;
      y[i] = y[i] / temp;
      if (i == k)  continue;
      for (int j = 0; j < n; j++)
        a[i][j] = a[i][j] - a[k][j];
      y[i] = y[i] - y[k];
    }
    k++;
  }
  
  for (k = n - 1; k >= 0; k--)
  {
    x[k] = y[k];
    for (int i = 0; i < k; i++)
      y[i] = y[i] - a[i][k] * x[k];
  }
  return x;
}


double* gausss(double** a, double* y, double n)
{
  double *x;
  
  x = gauss(a, y, n);
  
  return x;
}


int main(void)
{
	ofstream fout;
	fout.open("file2.txt");
	int i, j, N=6;
	double *points1, *points2, *points3, *points4, *values1, *values2, *values3, *values4, *v1, *v3, *x, *y, *L1, *L2, **matrix1, **matrix2;
	double *g1, *g3;
	double a=-4, b=4, h1=(b-a)/(N-1), h2=(b-a)/(3*(N-1)), h;

	points1 = new double[N];
	points2 = new double[3*(N-1)+1];
	points3 = new double[N];
	points4 = new double[3*(N-1)+1];
	values1 = new double[N];
	values2 = new double[3*(N-1)+1];
	values3 = new double[N];
	values4 = new double[3*(N-1)+1];  
	g1 = new double[3*(N-1)+1];
	g3 = new double[3*(N-1)+1];
	v1 = new double[N];
	v3 = new double[N];
	matrix1 = new double*[N];
	matrix2 = new double*[N];
	L1 = new double[3*(N-1)+1];
	L2 = new double[3*(N-1)+1];

	for(i=0;i<3*(N-1)+1;i++)
	{
		points2[i] = a+h2*i;
		values2[i] = func(points2[i]);
		//cout<<"points2 = "<<points2[i]<<endl;
	}

	for(i=0;i<N;i++)
	{
		matrix1[i] = new double[N];
		matrix2[i] = new double[N];

		points1[i] = a+h1*i; cout<<"points1 = "<<points1[i]<<endl;
        	points3[N-1-i]=0.5*(b+a)+0.5*(b-a)*cos((2*i+1)*PI/(2*N));
		//cout<< points1[i]<<endl;
	}//cout<<endl;

	for(i=0;i<N-1;i++)
	{
		//cout<< points2[i]<<endl;
		points4[3*i]=points3[i];
		h=(points3[i+1]-points3[i])/3;
		points4[3*i+1]=points3[i]+h;
		points4[3*i+2]=points3[i]+2*h;
	}
	points4[3*(N-1)]=points3[N-1];	
	
	for(i=0;i<3*(N-1)+1;i++)
		values4[i]=func(points4[i]);
	/*for(i=0;i<3*(N-1)+1;i++)
	{	
		cout<<points3[i]<<endl;
	}*/

	for(i=0;i<N;i++)
	{
		values1[i]=func(points1[i]); cout<<values1[i]<<endl;  
		values3[i]=func(points3[i]);
		
		v1[i]=values1[i];
		v3[i]=values3[i];
	}

	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			matrix1[i][j]=pow(points1[i],j); //cout<<matrix1[i][j]<<" ";
			matrix2[i][j]=pow(points3[i],j);
		}
	//cout<<endl;
	}		

	x=gausss(matrix1, v1, N);	//for(i=0;i<N;i++){cout<<"x"<<i<<" = "<< x[i]<<endl;}
	y=gausss(matrix2, v3, N);
	//cout<< lagr(values1, points1, 2, N);
	for(i=0;i<3*(N-1)+1;i++)
	{	
		g1[i]=g3[i]=0;		
			
		for(j=0;j<N;j++)
		{
			g1[i]+=x[j]*pow(points2[i],j);
			g3[i]+=y[j]*pow(points4[i],j);
		}
		//cout<<"g1 "<<i<<" = "<<g1[i]<<endl;
	}

	fout<<"Points L  "<<"Points G  "<<"L1        "<<"G1        "<<"f1        "<<"L2        "<<"G2        "<<"f2        "<<endl<<endl;
	fout<<endl;
	for(i=0;i<3*(N-1)+1;i++)
	{
		fout<<points2[i]<<" "<<points4[i]<<" "<<lagr(values1, points1, points2[i], N)<<"  "<<g1[i]<<" "<<values2[i]<<" "<<lagr(values3, points3, points4[i], N)<<"  "<<g3[i]<<" "<<values4[i]<<" "<<endl;
	}

	fout.close();


return 1;
}






