
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
#define nSIDEMAX 100001
#define bondcut 1.2
#define rcut 1.122462048309373
#define Nmax 10000
#define nneigh 20

#define epsilon 1.
double x[nSIDEMAX][2];
double pf[nSIDEMAX][2];
double rbond[nSIDEMAX];
double rbond_d[nSIDEMAX];
double sigma[nSIDEMAX];
double sigmaz[nSIDEMAX];
double msd[nSIDEMAX][2];
int bond[nSIDEMAX];
int parent[nSIDEMAX];
int lineage[nSIDEMAX];
int flag[nSIDEMAX];

double vpos[nSIDEMAX][2];
double skin=0.5;
double skinsq;
int nlist[nSIDEMAX];
int list[nSIDEMAX][nneigh];
int vNint;

double boxl,boxly,sigmac;
double grate;
double spring;
double Diff,press;
double dt;
int nsteps,stat,Nint,Nout;

char* prefix;

/* list functions for later use */
void read_input();
void initialize();
void integrate(double);
void pairforce();
void pairforcev();
void verletlist();

int periodic(int);
void sample(double);

float ran2(long *idum);

FILE* samplefile;
FILE* msdfile;
FILE* output;

/* BEGIN MAIN PROGRAM */
int main(int argc, char **argv) {
    int step;
    int n;
    double dx,dy,vr;
    
    output=fopen("position.xyz","w");
    samplefile=fopen("sample.txt","w");
    msdfile=fopen("msd.txt","w");

    /* read input parameters from file "input.txt" */
    read_input();
    sample(0);
    
    skinsq=skin*skin/4.;
    
    verletlist();
    
    for(step=0; step<=nsteps; step++){
        press=0;
        integrate(step);
        if(step%stat==0) {
            sample(step);
            /* Output step, total number of particles (mothers+daughtors), number that have left domain, and their difference and pressure on y boundary */

            printf("%lf %d %d %d %lf\n",step*dt,Nint,Nout,Nint-Nout,press);
            fprintf(samplefile,"%lf %d %d %d %lf\n",step*dt,Nint,Nout,Nint-Nout,press);
            fflush(samplefile);
        }
        
        int vflag=0;
        // Construct Verlet list if new particles
        if(vNint!=Nint){
            verletlist();
        }
        // Construct Verlet list if new displacement>skin/2
        else{
            for(n=0; n<Nint; n++){
                dx=vpos[n][0]-x[n][0];
                dy=vpos[n][1]-x[n][1];
                vr=dx*dx+dy*dy;
                if(vr>skinsq) {
                    vflag=1;
                    break;
                }

            }
            if(vflag==1) verletlist();
        }
        if(Nint>Nmax) break;
        
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
        // Lineage indicator
        lineage[n]=0;
        // mean squared displacement
        msd[n][0]=0;
        msd[n][1]=0;

    }

    /* Initialize particle positions */
    for (n=0; n<Nint; n++) {
        // Random initial mother position
        x[2*n][0]=boxl*ran2(&idum);
        x[2*n][1]=boxl*ran2(&idum);
        // Mother flag
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
        rbond[2*n]=0.01;
        rbond[2*n+1]=0.01;
        
        sigma[2*n]=1.;
        sigma[2*n+1]=0.01;
        flag[2*n]=1;
        flag[2*n+1]=1;
        
        msd[2*n][0]=x[2*n][0];
        msd[2*n][1]=0;
        
        lineage[2*n]=n;
        lineage[2*n+1]=n;
        
    }
    
    Nint=2*Nint;
    
}

void integrate(double step){
    int n,nn;
    double dbondx,dbondy,rbondt,r1,r2,etax,etay,force,bfx,bfy;
    
    // Compute vector of pair forces using verlet list
    pairforcev();
    
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
        press+=bfy;
        
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
        //rbond[n]+=dt*grate;
        //if(sigma[n]<sigmac) sigma[n]+=dt*grate;
        
        if(rbond[n]<1.) rbond[n]+=dt*grate*rbond[n]*4.5;
        else rbond[n]+=dt*grate;
        if(sigma[n]<sigmac) sigma[n]+=dt*grate*sigma[n]*4.5;
    }
    
    int Nintt=Nint;
    
    // Once daughtor is matured
    // sever bond and initialize new daughtors
    for(n=0; n<Nintt; n++) {
        if(rbond_d[n]>bondcut && sigma[n]>=sigmac && sigma[bond[n]]>=sigmac) {

            Nint+=1;
            flag[Nint-1]=1;
            
            r1=ran2(&idum);
            r2=ran2(&idum);

            etax=sqrt(-2.*log(r1))*cos(2.*pi*r2)/50.;
            etay=sqrt(-2.*log(r1))*sin(2.*pi*r2)/50.;
            
            x[Nint-1][0]=x[n][0]+etax;
            x[Nint-1][1]=x[n][1]+etay;
            
            lineage[Nint-1]=lineage[n];
            
            rbond[n]=0.01;
            rbond[Nint-1]=0.01;
            rbond_d[n]=0.01;
            rbond_d[Nint-1]=0.01;
            
            parent[n]=1;
            
            msd[n][0]=x[n][0];
            msd[n][1]=step;
            
            bond[n]=Nint-1;
            bond[Nint-1]=n;
            
            sigma[Nint-1]=0.01;
            
        }
    }
        
}

/* Construct Verlet list */
void verletlist(){
    
    int n,nn;
    double dx,dy,r,rsq;
    
    vNint=Nint;
    
    for(n=0; n<Nint; n++){
        nlist[n]=0;
        vpos[n][0]=x[n][0];
        vpos[n][1]=x[n][1];
        for(nn=0; nn<nneigh; nn++) list[n][nn]=0;
    }
    
    for(n=0; n<Nint; n++){
        if(flag[n]==0) continue;
        
        for(nn=0; nn<Nint; nn++){
            if(n==nn || bond[n]==nn) continue;
            if(flag[nn]==0) continue;
        
            dx=x[n][0]-x[nn][0];
            dy=x[n][1]-x[nn][1];
                
            rsq=(dx*dx+dy*dy);
            r=sqrt(rsq);
            
            if(r<rcut*sigmac+skin){
                list[n][nlist[n]]=nn;
                nlist[n]++;
            }
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

/* Compute vector of pair forces with Verlet List*/
void pairforcev(){
    double dx,dy,sig,sig12,sig6,r,rsq,rsqin,rsqin8,rsqin14,fmag;
    int n,nn,m;
    
    for(n=0; n<Nint; n++){
        if(flag[n]==0) continue;
        
        pf[n][0]=0.;
        pf[n][1]=0.;

        for(m=0; m<nlist[n]; m++){
            nn=list[n][m];
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

void sample(double step){
  int n;

  fprintf(output,"%d\nSpins\n",Nmax);
  for (n=0; n<Nint; n++) {
      if(flag[n]==1) fprintf(output,"O %lf %lf %lf %d\n",x[n][0],x[n][1],(sigma[n])/10.,lineage[n]);
      else fprintf(output,"O 0 0 0\n");
      
      if(parent[n]==1 && flag[n]==1) fprintf(msdfile,"%d %lf %lf %lf %lf\n",n,msd[n][1],step-msd[n][1],x[n][1],(x[n][0]-msd[n][0])*(x[n][0]-msd[n][0]));
      
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
		   

