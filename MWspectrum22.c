#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>

#define min(x, y)  (((x) < (y)) ? x : y)
#define max(x, y)  (((x) > (y)) ? x : y)
#define ITLIM	41
#define	ISLIM 20
#define	JLIM 40001
int main(int argc, char *argv[]) {
	FILE *opt;
	double mu[JLIM];
	double a[JLIM]; // phi, the allele frequency spectrum
	double Fmatrix[ITLIM][ISLIM];
	double kmatrix[ITLIM][ISLIM];
	int t,i,j,r,ne;
	int it,itmax,is,islow,ismax;
	double uterm,u,c;
	double Svec[51]={0.0,20.0,40.0,60.0,80.0,100.0,120.0,140.0,160.0,180.0,200.0,220.0,240.0,260.0,280.0,300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0};
	double theta,s,S;
	double sum1,sum2,prod;
	void fillmu(double mu[JLIM],int ne,double s);
	double product1(double mu[JLIM],int i,int j);
	double product2(double mu[JLIM],int i,int j);

	itmax=40;
	islow=1;
	ismax=50;
	opt=fopen("MWspectrum22.out","w");
	ne=2500;
	for (it=1;it<=itmax;it++) {
		theta=it;
		for(is=islow;is<=ismax;is++) {
			S=Svec[is];
			u=theta/(2.0*ne);
			uterm=theta/(1.0-u);
			s=S/ne;
			fillmu(mu,ne,s);
			for(i=1;i<2*ne;i++) a[i]=0.0;
			a[2*ne]=1.0;
			for(i=2*ne-1;i>=1;i--) {
				sum1=0.0;
				for (j=1;j<=min(i,2*ne-i);j++) {
					sum1+=uterm*product1(mu,i,j)*a[i+j]*(i+j)/(i*(2*ne-i)); }
				sum2=0.0;
				for (j=1+min(i,2*ne-i);j<=2*ne-i;j++) {
					sum2+=uterm*product2(mu,i,j)*a[i+j]*(i+j)/(i*(2*ne-i)); }
				a[i]=sum1+sum2;
				}
			c=0.0;
			for(i=1;i<=2*ne;i++)
				c+=i*a[i];
			for(i=1;i<=2*ne;i++)
				a[i]*=2*ne/c;  // this is the normalization step
			for(j=1;j<2*ne;j++)
				fprintf(opt,"%g\t",log(a[j]));
			fprintf(opt,"%g\n",log(a[j]));
			}
		} 
	exit(0);
	}

	double product1(double mu[20001],int i,int j) {
		int r;
		double prod;

		prod=1.0;
		for (r=1;r<=j;r++) prod*=(mu[i+r]/mu[j-r+1]);
		return prod;
		} 

	double product2(double mu[20001],int i, int j) {
		int r;
		double prod;

		prod=1.0;
		for (r=1;r<=i;r++) prod*=(mu[j+r]/mu[i-r+1]);
		return prod;
		} 

	void fillmu(double mu[10001],int ne,double s) {
		int j;
		for(j=1;j<=2*ne;j++) 
			mu[j]=1.0+s*j/(2.0*ne);
		return;
		}
		
