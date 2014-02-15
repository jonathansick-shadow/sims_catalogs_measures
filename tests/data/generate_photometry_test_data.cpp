#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define LIGHTSPEED 299792458.0
#define NM2M 1.0e-9
#define ERGSETC2JANSKY 1.0e23


main(){

double lambda_min=100.0,lambda_max=30000.0,lambda_step=10.0;

int nfilters=5,ii;

char *filters[5];

for(ii=0;ii<nfilters;ii++)filters[ii]=new char[10];

sprintf(filters[0],"u");
sprintf(filters[1],"g");
sprintf(filters[2],"r");
sprintf(filters[3],"i");
sprintf(filters[4],"z");

double mu[5];

mu[0]=lambda_min+200.0*lambda_step;
mu[1]=mu[0]+300.0*lambda_step;
mu[2]=mu[1]+700.0*lambda_step;
mu[3]=mu[2]+700.0*lambda_step;
mu[4]=mu[3]+500.0*lambda_step;

double sigma[5];

sigma[0]=100.0*lambda_step;
sigma[1]=200.0*lambda_step;
sigma[2]=200.0*lambda_step;
sigma[3]=150.0*lambda_step;
sigma[4]=200.0*lambda_step;

double _sb,norm,_phi,ll;

FILE *output;
char bname[100],aname[100];

int el=int((lambda_max-lambda_min)/lambda_step)+1;

double **phi,**sb,*lambda;

phi=new double*[nfilters];
sb=new double*[nfilters];

lambda=new double[el];
int i,j;
for(i=0;i<nfilters;i++){
    phi[i]=new double[el];
    sb[i]=new double[el];
}

for(ii=0;ii<nfilters;ii++){
    norm=0.0;
    sprintf(bname,"test_bandpass_%s.dat",filters[ii]);
    sprintf(aname,"test_phi_%s.dat",filters[ii]);
    
    output=fopen(bname,"w");
    
    i=0;
    for(ll=lambda_min;ll<lambda_max+1.0;ll+=lambda_step){
        if(i>=el){
            printf("WARNING i overstepped sb %d %d %e %e\n",i,el,ll,lambda_max);
            exit(1);
        }
        
        _sb=exp(-0.5*(ll-mu[ii])*(ll-mu[ii])/(sigma[ii]*sigma[ii]));
        norm+=lambda_step*_sb/ll;
        
        fprintf(output,"%.18e %.18e\n",ll,_sb);
        sb[ii][i]=_sb;
        lambda[i]=ll;
        i++;
        
    }
    
    fclose(output);

    output=fopen(aname,"w");
    i=0;
    for(ll=lambda_min;ll<lambda_max+1.0;ll+=lambda_step){
        if(i>=el){
            printf("WARNING i overstepped el %d %d %e %e\n",i,el,ll,lambda_max);
            exit(1);
        }
        
        _sb=exp(-0.5*(ll-mu[ii])*(ll-mu[ii])/(sigma[ii]*sigma[ii]));
        _phi=_sb/(ll*norm);
        
        fprintf(output,"%.18e %.18e\n",ll,_phi);
        phi[ii][i]=_phi;
        i++;
        
    }
    fclose(output);
}

int n_sed=10;
char sedname[100];

double **sed;
sed=new double*[n_sed];
for(ii=0;ii<n_sed;ii++)sed[ii]=new double[el];

for(ii=0;ii<n_sed;ii++){
    
    sprintf(sedname,"test_sed_%d.dat",ii);
    
    output=fopen(sedname,"w");
    i=0;
    for(ll=lambda_min;ll<lambda_max+1.0;ll+=lambda_step){
        if(i>=el){
            printf("WARNING i overstepped in sed %d %d\n",i,el);
        }
        
        if(ii<5){
            sed[ii][i]=1.6+atan((ll-ii*3000.0)/2000.0);
        }
        else{
            //sed[ii][i]=2.0*exp(-ll/(ii*500.0))+0.5
            
            sed[ii][i]=1.6-atan((ll-exp(1.9*log(ii-4.0))*1000.0)/2000.0);
            
        }  
        
        fprintf(output,"%.18e %.18e\n",ll,sed[ii][i]);
        
        i++;
        
    }
    
    fclose(output);
}

double fnu,flux;
int ifilter,ised;

output=fopen("test_magnitudes.dat","w");

for(ifilter=0;ifilter<nfilters;ifilter++){
    for(ised=0;ised<n_sed;ised++){
        

        
        flux=0.0;
        for(i=0;i<el;i++){
            fnu=sed[ised][i]*lambda[i]*lambda[i]*NM2M*ERGSETC2JANSKY/LIGHTSPEED;
            flux+=fnu*phi[ifilter][i]*lambda_step;
        }
        
        fprintf(output,"filter %s sed %d mag %.18e\n",
        filters[ifilter],ised,-2.5*log(flux)/log(10.0)+8.9);
        
    }
}

fclose(output);

}
