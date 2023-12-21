#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

#define rho 0.8442
#define V N*1./rho
#define N          108
#define TEMP       0.728
#define CALC       1
#define BLOCKS     100    // num part para calc msd
#define TEMPO_MAX  20000
#define DT         0.01
#define numMSD 200

/*    Funcoes   */

void inicializacao(void);
void forca_LJ(void);
void integra_vel_verlet(void);
void calcula_gr(void);
void imprime_gr(void);
void inicializa_msd(void);
void measure_msd(long int steps);


/*Variaveis Globais   */

double rx[N],ry[N],rz[N],vx[N],vz[N],vy[N],ax[N],ay[N],az[N];
double gr[N],ax1[N],ay1[N],az1[N],ene_pot,ene_kin;
double vcmx,vcmy,vcmz,fs;
double tempo;
int ngr;
int M,n,itt;
double L = pow(V,1./3.);


/***************Parametros para calcular MSD******************/	

double msd[numMSD], msdx[numMSD],msdy[numMSD],msdz[numMSD];
double vacf[numMSD],vacfx[numMSD],vacfy[numMSD],vacfz[numMSD];
double rx0[N],ry0[N],rz0[N],vx0[N],vy0[N],vz0[N];
/*********************************/	
int main (void)
/*********************************/	
{

FILE *fpdados;

int i;

inicializacao();
forca_LJ();

fpdados = fopen("energias.dat","w");
tempo =0.0; itt = 0;

for(itt=0;itt<TEMPO_MAX;itt++)
{

integra_vel_verlet();

fprintf(fpdados,"%.4f %f %f %f\n",tempo,ene_kin/N,ene_pot/N,(ene_kin+ene_pot)/N);
    
if((itt%CALC)==0)
{
calcula_gr();
}

measure_msd(itt);

tempo = tempo + DT;
 
  }
  
measure_msd(itt);

fclose(fpdados);
 
return 0;
}

/*********************************/	
void inicializacao (void)
/*********************************/	
{

int i,j;
double dx,dy,dz,v2,vcmx1,vcmy1,vcmz1;
  

ngr = 0;
vcmx = 0.0;
vcmy = 0.0;
vcmz = 0.0;
fs = 0.0;
v2 = 0.0;

for(i=0; i<N; i++)
  {
gr[i]=0.0;
  }

n=ceil(pow(N,1./3.));
M=pow(n,3);

dx=pow(V/M,(1./3.));
dy=pow(V/M,(1./3.));
dz=pow(V/M,(1./3.));
  
for(i=0; i<N; i++)
 {

rx[i] = (i%n)*dx;
ry[i] = ((i/n)%n)*dy;
rz[i] = ((i/(n*n))%n)*dz; 
      	
vx[i] = (1.*rand()/RAND_MAX)-0.5;
vy[i] = (1.*rand()/RAND_MAX)-0.5;
vz[i] = (1.*rand()/RAND_MAX)-0.5;
	
vcmx += vx[i];
vcmy += vy[i];
vcmz += vz[i];

v2 = v2 + vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];


  }
vcmx = vcmx/(N);
vcmy = vcmy/(N);
vcmz = vcmz/(N);

v2 = v2/(N);
    
fs = sqrt(3.0*TEMP/v2);

  
for(i=1; i<N; i++)
   {
vx[i] = (vx[i]-vcmx)*fs;
vy[i] = (vy[i]-vcmy)*fs;
vz[i] = (vz[i]-vcmz)*fs;
       
   }

// inicializa o msd
inicializa_msd();

return;

}

 /*********************************/	
void forca_LJ(void)
/*********************************/	
{

int i,j;
double dx,dy,dz,r2,r6i,r2i,f;

double rc,rc2,ecut;

rc = (L)/2;
rc2 = rc*rc;
ecut = 4*(pow(rc2,-6)-pow(rc2,-3));
ene_pot = 0.0;

for(i=0;i<N;i++)
 {
ax[i]= 0.0;
ay[i]= 0.0;
az[i]= 0.0;
 }
  
  
for(i=0;i<(N-1);i++)
{
for(j=(i+1);j<N;j++) 
  {

dx = rx[i]-rx[j];
dy = ry[i]-ry[j];
dz = rz[i]-rz[j];
      
dx=dx-rint(dx/L)*L;
dy=dy-rint(dy/L)*L;
dz=dz-rint(dz/L)*L;
      
r2 = dx*dx + dy*dy + dz*dz;
         
if(r2 < rc2)
{
      
r2i = 1.0/r2;
r6i = r2i*r2i*r2i;
f = 48.0*r6i*r2i*(r6i - 0.5);
      
// componente x das aceleracoes
ax[i] = ax[i] + f*dx;
ax[j] = ax[j] - f*dx;

// componente y das aceleracoes
ay[i] = ay[i] + f*dy;
ay[j] = ay[j] - f*dy;
      
// componente z das aceleracoes
az[i] = az[i] + f*dz;
az[j] = az[j] - f*dz;
      
ene_pot += 4.0*r6i*(r6i-1.0)-ecut;

      }
    }

  }

return;
  
}
	
void integra_vel_verlet(void)
/*********************************/	
{

int i;
double dth;

dth = 0.5*DT;

for(i=0;i<N;i++)
{

//atualizao as velocidade ate DT/2
    
vx[i] = vx[i] + ax[i]*dth;
vy[i] = vy[i] + ay[i]*dth;
vz[i] = vz[i] + az[i]*dth;

// atualizao as posicoes ate DT
rx[i] = rx[i] + vx[i]*DT;
ry[i] = ry[i] + vy[i]*DT;
rz[i] = rz[i] + vz[i]*DT;
  }


forca_LJ();

ene_kin = 0.0;

vcmx = 0.0; 
vcmy = 0.0,
vcmz = 0.0;

for(i=0;i<N;i++)
 {
//atualizao as velocidade ate DT
vx[i] = vx[i] + ax[i]*dth;
vy[i] = vy[i] + ay[i]*dth;
vz[i] = vz[i] + az[i]*dth;
    
vcmx += vx[i];
vcmy += vy[i];
vcmz += vz[i];
    
ene_kin = ene_kin + (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])/2;
  }
  
 
vcmx = vcmx/(1.0*N);
vcmy = vcmy/(1.0*N);
vcmz = vcmz/(1.0*N);
  
  return;

}
/*********************************/	

/*********************************/	
void calcula_gr(void)
/*********************************/	
{

int i,j,ig;
double delg,dx,dy,dz,r2,r1;

//lagura do bin
delg = L/(2.0*N);
//contador no numero de chamadas da funcao
ngr = ngr + 1;

for(i=0;i<(N-1);i++)
{
for(j=(i+1);j<N;j++)
  {
//minimum image convention
dx = rx[i]-rx[j];
dx=dx-rint(dx/L)*L;
dy = ry[i]-ry[j];
dy=dy-rint(dy/L)*L;
dz = rz[i]-rz[j];
dz=dz-rint(dz/L)*L;

r2 = dx*dx + dy*dy + dz*dz;
r1 = sqrt(r2);
 
if(r1<(0.5*L))
    {
ig = rint(r1/delg);
gr[ig] = gr[ig] + 2; //atualizo i E j
    }
  }
}

return;

}

/*********************************/	
void imprime_gr(void)
/*********************************/	
{

FILE *fpgr;

int i;
double delg,r,vb,nid;

fpgr=fopen("gr.dat","w");
delg = L/(2.0*N);

for(i=0;i<N;i++)
{

r = delg*(i+0.5);
vb=(pow(((i+1)*delg),3)-pow((i*delg),3));
nid = (4./3.)*M_PI*vb*rho;
gr[i] = gr[i]/(1.0*ngr*N*nid);
fprintf(fpgr,"%d %f %f\n",i,r,gr[i]);
}

fclose(fpgr);
 
return;

}

/*********************************/	
void inicializa_msd(void)
/*********************************/	
{

int i;



//guardo os vetores das "posições" e "velocidades" iniciais
for(i=0;i<N;i++)
 {
rx0[i] = rx[i];
ry0[i] = ry[i];
rz0[i] = rz[i];
vx0[i] = vx[i];
vy0[i] = vy[i];
vz0[i] = vz[i]; 
  }

 
for(i=0;i<numMSD;i++) 
  {
msd[i] = 0.0;
msdx[i] = 0.0;
msdy[i] = 0.0;
msdz[i] = 0.0;
vacf[i] = 0.0;
vacfx[i] = 0.0;
vacfy[i] = 0.0;
vacfz[i] = 0.0;

  }

  return;
}
/*********************************/	
void measure_msd(long int steps)
/*********************************/	
{

int i,ref;
double mmsdx,mmsdy,mmsdz;
double mvacfx,mvacfy,mvacfz;


//calculo o tempo dentro de um bloco
ref = steps%numMSD;
//novos r0 e v0
if(ref==0)
   {
//redefino os pontosiniciais para msd
for(i=0;i<N;i++) 
{
rx0[i] = rx[i];
ry0[i] = ry[i] ;
rz0[i] = rz[i];

vx0[i] = vx[i];
vy0[i] = vy[i];
vz0[i] = vz[i];
    }

  }
   
  
mmsdx = 0.0; mvacfx = 0.0;
mmsdy = 0.0; mvacfy = 0.0;
mmsdz = 0.0; mvacfz = 0.0;
  
for(i=0;i<N;i++)
   {
    
    
mmsdx = mmsdx + (rx[i]-rx0[i])*(rx[i]-rx0[i]);
mmsdy = mmsdy + (ry[i]-ry0[i])*(ry[i]-ry0[i]);
mmsdz = mmsdz + (rz[i]-rz0[i])*(rz[i]-rz0[i]);
    
mvacfx = mvacfx + vx[i]*vx0[i];
mvacfy = mvacfy + vy[i]*vy0[i];
mvacfz = mvacfz + vz[i]*vz0[i];

  }
  
mmsdx = mmsdx/(1.0*N);
mmsdy = mmsdy/(1.0*N);
mmsdz = mmsdz/(1.0*N);
  
msdx[ref] = msdx[ref] + mmsdx;
msdy[ref] = msdy[ref] + mmsdy;
msdz[ref] = msdz[ref] + mmsdz;  
    
msd[ref] = msd[ref] + mmsdx + mmsdy + mmsdz;

mvacfx = mvacfx/(1.0*N);
mvacfy = mvacfy/(1.0*N);
mvacfz = mvacfz/(1.0*N);
  
vacfx[ref] = vacfx[ref] + mvacfx;
vacfy[ref] = vacfy[ref] + mvacfy;  
vacfz[ref] = vacfz[ref] + mvacfz;  

vacf[ref] = vacf[ref] + mvacfx + mvacfy +mvacfz;
  
/******************CALCULA MSD*******************/  

FILE *fpmsd;

double diff,dif;

diff = 0.0;  
dif =0.0;
fpmsd = fopen("fmsd.dat","w");
fprintf(fpmsd,"#part time msd vcaf diff\n");

for(i=0;i<numMSD;i++)
  
   {


fprintf(fpmsd,"%f\t%f\t%f\t%f\n",i*DT,msd[i]/BLOCKS,vacf[i]/BLOCKS,diff);
    
    
/******integro o vacf para obter Diff********/
    
diff = diff + vacf[i]*DT/(3*BLOCKS);
  }
  
  
 fprintf(fpmsd,"#diff = %f\n",diff);

fclose(fpmsd);  
  
  

return;
}
