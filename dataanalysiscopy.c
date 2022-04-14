#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include <ctype.h>

#define ITLIM 41
#define ISLIM 52
#define JLIM 5001
#define ALLELELIMIT 201
#define SAMLIM 2001
#define STATLIM 10  
double loga[ITLIM][ISLIM][JLIM];
double blog[ITLIM][ISLIM][2*SAMLIM];
double temploga[5][JLIM],templogb[5][JLIM];
int itmax=40,ismax=50,ne=2500;
double logfactorial[SAMLIM];
int main(int argc, char *argv[]) {
	FILE *ipt,*spt;
	int alleleno,sampleno,maxsample,maxallele[SAMLIM],count[SAMLIM][ALLELELIMIT],i,j,length;
	int idum,samsize,kmax=100,sam[SAMLIM],statno;
	double **outmatrix,**b;
	char line[401],num[12]; 
	char name[80];
	void readloga(int itmax,int ismax,int ne);
	void filllogfactorial(double logfactorial[SAMLIM],int n);
	double **matrix(long nrl, long nrh, long ncl, long nch);
	void statistics(double **outmatrix,int sampleno,int sam[ALLELELIMIT],int maxallele,int samsize,double **b,int *idum);
	void fillb(double **b,int kmax,int n);
	
	ipt=fopen(argv[1],"r");
	sampleno=0;
	while(1) {
		fgets(line,400,ipt);
		if(feof(ipt)) {
			printf("End of input file reached\n");
			break;
			}
		else {
			sampleno++;
			length=strlen(line);
			alleleno=i=0;
			while(i<length) {
				for(j=0;j<10;j++) 
					num[j]='\0';
				while(isspace(line[i]))
					i++;
				j=0;
				while(isdigit(line[i])) {
					num[j]=line[i];
					i++;
					j++;
					}
				alleleno++;
				if(isdigit(num[0]))
					count[sampleno][alleleno]=atoi(num);
				maxallele[sampleno]=alleleno-1; 
				}  // end while i
			} // end else
		} // while (1)
	maxsample=sampleno;
	printf("%d samples read from %s\n",maxsample,argv[1]);
	for(sampleno=1;sampleno<=maxsample;sampleno++) {
		for(alleleno=1;alleleno<=maxallele[sampleno];alleleno++)
			printf("%d ",count[sampleno][alleleno]);
		printf("\n");
		}

	readloga(itmax,ismax,ne);
	filllogfactorial(logfactorial,201);
	outmatrix=matrix(1,maxsample,1,STATLIM); 
	idum=-time(NULL);
	b = matrix(1,kmax,1,201); 
	fillb(b,kmax,201);
	for(sampleno=1;sampleno<=maxsample;sampleno++) {
			for(samsize=0,alleleno=1;alleleno<=maxallele[sampleno];alleleno++) {
				sam[alleleno-1]=count[sampleno][alleleno];
				samsize+=sam[alleleno-1];
				}  // end for alleleno
				samsize/=2;
				statistics(outmatrix,sampleno,sam,maxallele[sampleno],samsize,b,&idum);
				printf("Finished analyzing sample %d\n",sampleno);
			}  // end for sampleno
		sprintf(name,"%s.output",argv[1]);
		spt=fopen(name,"w");
		for(sampleno=1;sampleno<=maxsample;sampleno++) 
			for(statno=1;statno<=9;statno++) {
				fprintf(spt,"%g",outmatrix[sampleno][statno]);
				if(statno<9) fprintf(spt,"\t");
					else fprintf(spt,"\n");
				}
		exit(1);
	}

void fillblog(int samsize) {
	int it,is,j,i;
	double x,term,tempb[2*SAMLIM],bsum;

  for(it=1;it<=itmax;it++)
    for(j=1;j<=2*ne;j++)
      loga[it][0][j]=log(it/(1.0+j));
  for(it=1;it<=itmax;it++)
    for(is=0;is<=ismax; is++) {
      for(j=0;j<2*samsize;j++) {
        bsum=0.0;
        for(i=1;i<2*ne;i++) {
          x=i/(2.0*ne);
          term=exp(loga[it][is][i]+logfactorial[2*samsize]-logfactorial[j]-logfactorial[2*samsize-j]+j*log(x)+(2*samsize-j)*log(1-x));
          bsum+=term;
          }
        tempb[j]=bsum;
        }
      bsum=0.0;
      for(j=1;j<2*samsize;j++)
        bsum+=tempb[j];
      for(j=1;j<2*samsize;j++)
        blog[it][is][j]=log(tempb[j]/bsum);
      } // end, for it, is
	}
	
void readloga(int itmax,int ismax,int ne) {
	FILE *apt;
  int it,is,j;

	apt=fopen("MWspectrum22.out","r"); // unnormalized spectrum
  for(it=1;it<=itmax;it++)
    for(is=1;is<=ismax; is++) {
      for(j=1;j<=2*ne;j++) {
        fscanf(apt,"%lf",&loga[it][is][j]);
        } 
      } 
  } 

void filllogfactorial(double logfactorial[SAMLIM],int n) {
  int i;

  logfactorial[0]=logfactorial[1]=0.0;
  for (i=2;i<=2*n;i++)
    logfactorial[i]=log(i)+logfactorial[i-1];
  return;
  }

double loglikelihood(int spec[JLIM],int specmax,double blog[2*SAMLIM]) {
  int j;
  double sum;

  sum=0.0;
  for(j=1;j<=specmax;j++)
    if(spec[j]) {
      sum+=spec[j]*blog[j];
      }
  return sum;
  } 

void MLE(double *thetapt,double *spt,int spec[ALLELELIMIT],int maxspec) {
	int it,is,tval,sval;
	int itmax=40,ismax=50;
	double logL,logLmax=-10000000.0;
	double loglikelihood(int spec[JLIM],int specmax,double blog[2*SAMLIM]);
	double Svec[51]={0.0,20.0,40.0,60.0,80.0,100.0,120.0,140.0,160.0,180.0,200.0,
		220.0,240.0,260.0,280.0,300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,
		520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,
		820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0};

  tval=sval=0;
  for(it=1;it<=itmax;it++)
    for(is=0;is<=ismax;is++) {
    logL=loglikelihood(spec,maxspec,blog[it][is]);
    if(logL>logLmax) {
      tval=it;
      sval=is;
      logLmax=logL;
      }
    } 
	*thetapt=tval;
	*spt=Svec[sval];
	}

void spectrum(int spec[ALLELELIMIT],int *spt,int sam[ALLELELIMIT],int maxallele) {
	int i;

	for(i=1;i<ALLELELIMIT;i++)
		spec[i]=0;
	*spt=0;
	for(i=0;i<maxallele;i++) {
		if(sam[i]>*spt) *spt=sam[i];
		if(sam[i]) spec[sam[i]]++;
		}
	}
	
void statistics(double **outmatrix,int sampleno,int sam[ALLELELIMIT],int maxallele,int samsize,double **b,int *idum) {
	int k,specmax,i,spec[ALLELELIMIT],maxrep=1000,repno,Pcount;
	double F,Pval,theta,s;
	void spectrum(int spec[ALLELELIMIT],int *spt,int sam[ALLELELIMIT],int maxallele);
	double generateF(int k,int n,double **b,int *idum);
	void MLE(double *thetapt,double *spt,int spec[ALLELELIMIT],int maxspec);
	void fillblog(int samsize);
	
	fillblog(samsize);
	spectrum(spec,&specmax,sam,maxallele);
	F=0.0;
	for(i=0;i<maxallele;i++)
		F+=pow(sam[i]/(2.0*samsize),2.0);
	k=0;
	for(i=1;i<=specmax;i++)
		k+=spec[i];
	Pcount=0;
	for(repno=1;repno<=maxrep;repno++) // replicates are for the EW test
		if(generateF(k,2*samsize,b,idum)<F) Pcount++;
	Pval=(double)Pcount/maxrep;
	MLE(&theta,&s,spec,specmax);
	outmatrix[sampleno][1]=k;
	outmatrix[sampleno][2]=spec[1];
	outmatrix[sampleno][3]=spec[2];
	outmatrix[sampleno][4]=spec[3];
	outmatrix[sampleno][5]=F;
	outmatrix[sampleno][6]=Pval;
	outmatrix[sampleno][7]=theta;
	outmatrix[sampleno][8]=s;
	outmatrix[sampleno][9]=k*F;
	outmatrix[sampleno][10]=samsize;
	printf("sampleno=%d: ",sampleno);
	printf("{%d,%d,%d,%d,%g,%g,%g,%g,%g,%d}\n",k,spec[1],spec[2],spec[3],F,Pval,theta,s,k*F,samsize);
	} // end statistics

void fillb(double **b,int kmax,int n) {
	int i,j;

  for (j=1; j<=n; j++)
    b[1][j] = 1.0 / j;
  for (i=2; i<=kmax; i++)  {
    b[i][i] = 1.0;
    for (j=i; j<n; j++)
      b[i][j+1] = (i * b[i-1][j] + j * b[i][j]) / (j + 1.0);
    }
	}

double generateF(int k,int n,double **b,int *idum)  {
  double ran1(int *idum),cum,F,*ranvec;
  int i,l,nleft,*r;
	int *ivector(long nl,long nh);
	double *vector(long nl,long nh);

	r=ivector(1,k);
	ranvec=vector(1,k);
  for (i=1; i<=k-1; i++)
    ranvec[i] = ran1(idum);
  nleft = n;
  for (l=1; l<k; l++)  {
    cum = 0.0;
    for (i=1; i<=nleft; i++) {
      cum += b[k-l][nleft-i] / (i * b[k-l+1][nleft]);
      if (cum >= ranvec[l]) break;
      }
    r[l] = i;
    nleft -= i;
    }
  r[k] = nleft;
	F=0.0;
	for(i=1;i<=k;i++) {
		F+=pow((double) r[i]/n,2.0);;
		}
	return F;
  }

#define NR_END 1
double *vector(long nl, long nh)
{
  double *v;
	void nrerror(char error_text[]);

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

int *ivector(long nl, long nh)
{
  int *v;
	void nrerror(char error_text[]);

  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl+NR_END;
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.141592654
#define min(x, y)  (((x) < (y)) ? x : y)
#define max(x, y)  (((x) > (y)) ? x : y)

double ran1(int *idum)  {
    int j;
    int k;
    static int iy=0;
    static int iv[NTAB];
    double temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;
	void nrerror(char error_text[]);

  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m;
}

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

