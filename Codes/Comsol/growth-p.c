
/* insert libraries defining standard functions */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#define pi 3.141592653589793
/* ran2 is a function that draws a number between 0 and 1 
at "random" from a uniform distribution.  The syntax for
usage is ran2(&idum). */
long int idum;

/* set maximum size of system */
/* SIDEMAX is maximum allowable number of spins per side */

#define tSIDEMAX 15001
#define nSIDEMAX 10001
#define bondcut 1.2
#define rcut 1.122462048309373
#define Nmax 10000
#define epsilon 1.
double x[nSIDEMAX][2];
double pf[nSIDEMAX][2];
double rbond[nSIDEMAX];
double rbond_d[nSIDEMAX];
double sigma[nSIDEMAX];
int bond[nSIDEMAX];
int parent[nSIDEMAX];
int flag[nSIDEMAX];

double boxl,boxly,sigmac;
double grate;
double spring;
double Diff;
double dt;
int nsteps,stat,Nint,Nout;

char* prefix;

/* list functions for later use */
void read_input();
void initialize();
void integrate();
void pairforce();

int periodic(int);
void sample();

float ran2(long *idum);

FILE* samplefile;
FILE* output;

/* BEGIN MAIN PROGRAM */
int main(int argc, char **argv) {
    int step;
    
    output=fopen("position.xyz","w");
    samplefile=fopen("sample.txt","w");

    /* read input parameters from file "input.txt" */
    read_input();
    sample();
    
    for(step=0; step<=nsteps; step++){
        integrate();
        if(step%stat==0) {
            sample();
            /* Output step, total number of particles (mothers+daughtors), number that have left domain, and their difference */
            printf("%d %d %d %d\n",step,Nint,Nout,Nint-Nout);
            fprintf(samplefile,"%d %d %d %d\n",step,Nint,Nout,Nint-Nout);
            fflush(samplefile);
        }

        
    }
    
}


/* read input parameters from file "input.txt" */
void read_input() {
    char line[101];
    int max_length = 100;
    FILE * input_file;
    input_file = fopen("input.txt","r");
    // Number of integration steps
    fgets(line,max_length,input_file);
    sscanf(line,"%d",&nsteps);
    // Output frequency
    fgets(line,max_length,input_file);
    sscanf(line,"%d",&stat);
    // Initial number of mother particles
    fgets(line,max_length,input_file);
    sscanf(line,"%d",&Nint);
    // Mature particle diameter
    fgets(line,max_length,input_file);
    sscanf(line,"%lf",&sigmac);
    // X-domain size
    fgets(line,max_length,input_file);
    sscanf(line,"%lf",&boxl);
    // Y-domain size
    fgets(line,max_length,input_file);
    sscanf(line,"%lf",&boxly);
    // Growth rate
    fgets(line,max_length,input_file);
    sscanf(line,"%lf",&grate);
    // Mother-daughtor spring stiffness
    fgets(line,max_length,input_file);
    sscanf(line,"%lf",&spring);
    // Diffusion constant
    fgets(line,max_length,input_file);
    sscanf(line,"%lf",&Diff);
    // Integartion timestep
    fgets(line,max_length,input_file);
    sscanf(line,"%lf",&dt);

    printf("Initial N %d\n",Nint);
    printf("Growth Rate %lf\n",grate);
    
    initialize();
    
}

/* Intialize particle positions and flags*/
void initialize() {
    int n;
    double r1,r2,etax,etay;
    /* use time on clock to initialize random number generator */
    idum=-time(NULL);
    Nout=0;
    
    /* initialize variables */
    for (n=0; n<Nmax; n++) {
        // vector of ids for mother-daughtor pairs
        bond[n]=0;
        // Bond length vector
        rbond[n]=0.0;
        // Sigma vector
        sigma[n]=0.;
        // 1 for mother 0 for daughtor
        parent[n]=0.;
        // Flag to integrate 1 stil in domain
        // 0 outside of the domain
        flag[n]=0;
    }

    /* Initialize particle positions */
    for (n=0; n<Nint; n++) {
        
        // Random initial mother position
        x[2*n][0]=boxl*ran2(&idum);
        x[2*n][1]=boxl*ran2(&idum);
        
        // Mother-daughtor flag
        parent[2*n]=1;
        
        r1=ran2(&idum);
        r2=ran2(&idum);
        
        // Mother-daughtor pair flags
        bond[2*n]=2*n+1;
        bond[2*n+1]=2*n;

        etax=sqrt(-2.*log(r1))*cos(2.*pi*r2)/50.;
        etay=sqrt(-2.*log(r1))*sin(2.*pi*r2)/50.;
        
        // Random initial daughtor position
        x[2*n+1][0]=x[2*n][0]+etax;
        x[2*n+1][1]=x[2*n][1]+etay;
        
        // Bond vector
        rbond_d[2*n]=sqrt(etax*etax+etay*etay);
        rbond_d[2*n+1]=sqrt(etax*etax+etay*etay);

        // Initial bond equilibrium position
        rbond[2*n]=0.0;
        rbond[2*n+1]=0.0;
        
        sigma[2*n]=1.;
        sigma[2*n+1]=0.;
        flag[2*n]=1;
        flag[2*n+1]=1;
        printf("%lf %lf\n",x[2*n][0],x[2*n][1]);
    }
    
    Nint=2*Nint;
    
}

void integrate(){
    int n,nn;
    double dbondx,dbondy,rbondt,r1,r2,etax,etay,force,bfx,bfy;
    
    // Compute vector of pair forces
    pairforce();
    
    for(n=0; n<Nint; n++){
        // Check if inside domain
        if(flag[n]==0) continue;
        
        dbondx=x[n][0]-x[bond[n]][0];
        dbondy=x[n][1]-x[bond[n]][1];
        
        rbondt=(dbondx*dbondx+dbondy*dbondy);
        rbondt=sqrt(rbondt);
        
        rbond_d[n]=rbondt;
        rbond_d[bond[n]]=rbondt;
        
        // Bond force
        force=-spring*(rbondt-rbond[n]);
        
        r1=ran2(&idum);
        r2=ran2(&idum);
        
        // Random force
        etax=sqrt(-2.*log(r1))*cos(2.*pi*r2);
        etay=sqrt(-2.*log(r1))*sin(2.*pi*r2);
        
        // Exponentially repulsive wall force
        bfx=-exp(-(boxl-x[n][0])/0.1)/0.1;
        bfx+=exp(-(x[n][0])/0.1)/0.1;
        bfy=exp(-(x[n][1])/0.1)/0.1;

        // Displacement
        double ddx=(pf[n][0]+bfx)*dt+dt*force*dbondx/rbondt+sqrt(2.*Diff*dt)*etax;
        double ddy=(pf[n][1]+bfy)*dt+dt*force*dbondy/rbondt+sqrt(2.*Diff*dt)*etay;
        
        // First order Euler integration
        x[n][0]+=ddx;
        x[n][1]+=ddy;
        
        // Particles escape if y>Boxly
        if(x[n][1]>boxly) {
            flag[n]=0;
            flag[bond[n]]=0;
            Nout+=2;
        }
    }
    
    // Grow daughtor by increasing \sigma and equil bond displacement
    for(n=0; n<Nint; n++) {
        if(flag[n]==0) continue;
        rbond[n]+=dt*grate;
        if(sigma[n]<sigmac) sigma[n]+=dt*grate;
    }
    
    int Nintt=Nint;
    
    // Once daughtor is matured
    // sever bond and initialize new daughtors
    for(n=0; n<Nintt; n++) {
        if(rbond_d[n]>bondcut) {
            Nint+=1;
            flag[Nint-1]=1;
            
            r1=ran2(&idum);
            r2=ran2(&idum);

            etax=sqrt(-2.*log(r1))*cos(2.*pi*r2)/50.;
            etay=sqrt(-2.*log(r1))*sin(2.*pi*r2)/50.;
            
            x[Nint-1][0]=x[n][0]+etax;
            x[Nint-1][1]=x[n][1]+etay;
            
            rbond[n]=0.0;
            rbond[Nint-1]=0.0;
            rbond_d[n]=0.0;
            rbond_d[Nint-1]=0.0;
            
            parent[n]=1;
            
            bond[n]=Nint-1;
            bond[Nint-1]=n;
            
            sigma[Nint-1]=0.0;
            
        }
    }
    
    
    
}

/* Compute vector of pair forces */
void pairforce(){
    double dx,dy,sig,sig12,sig6,r,rsq,rsqin,rsqin8,rsqin14,fmag;
    int n,nn;
    
    for(n=0; n<Nint; n++){
        if(flag[n]==0) continue;
        
        pf[n][0]=0.;
        pf[n][1]=0.;

        for(nn=0; nn<Nint; nn++){
            // Skip if bonded or outside domain
            if(n==nn || bond[n]==nn) continue;
            if(flag[nn]==0) continue;
            
            dx=x[n][0]-x[nn][0];
            dy=x[n][1]-x[nn][1];
                
            rsq=(dx*dx+dy*dy);
            r=sqrt(rsq);
            sig=(sigma[n]+sigma[nn])/2.;
            
            if (r<rcut*sig){
                sig6=sig*sig*sig*sig*sig*sig;
                sig12=sig6*sig6;
                rsqin=1./rsq;
                rsqin8=rsqin*rsqin*rsqin*rsqin;
                rsqin14=rsqin8*rsqin*rsqin*rsqin;
                
                fmag=4.*epsilon*(12.*sig12*rsqin14-6.*sig12*rsqin8);
                pf[n][0]+=fmag*dx;
                pf[n][1]+=fmag*dy;
            }
        
        }
    
    }
}


void sample(){
  int n;

  fprintf(output,"%d\nSpins\n",Nmax);
  for (n=0; n<Nint; n++) {
		fprintf(output,"O %lf %lf %lf\n",x[n][0],x[n][1],parent[n]/10.+(1-flag[n]));
  }
    for (n=Nint; n<Nmax; n++) {
        fprintf(output,"O 0 0 0\n");
    }

  fflush(samplefile);
    
 
}




/* ran2 is a function that draws a number between 0 and 1
at "random" from a uniform distribution.  The syntax for
usage is ran2(&idum). */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2 = (*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy+= IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
		   

