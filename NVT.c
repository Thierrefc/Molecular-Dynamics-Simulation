#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>

#define Nc 24

void Initialize_Positions(int Np, double V, double *x,double *y,double *xo,double *yo,double *xd, double *yd, double *vx,double *vy, double Lx, double Ly, double dt);

void Initialize_Velocities(int Np, double dt, double To, double *vx,double *vy, double *T, double *Ec);

void Compute_Forces(int Np, double V, double rho, double *T, double *fx,double *fy,double *x,double *y,double *Ep, double Rc2, double Lx, double Ly, double *P, double ParametrosPotencial[3][3]);

void Integrate_Velocity_Verlet_pt1(int Np, double *x, double *y, double *vx, double *vy, double *fx, double *fy, double dt, double *Ec, double *Ep, double Rc2, double Lx, double Ly, int *counter, double *P, double *xd, double *yd, double *MSD);

void Integrate_Velocity_Verlet_pt2(int Np, double *x, double *y, double *vx, double*vy, double *fx, double *fy, double dt, double *Ec, double *Ep, double *T, double Rc2, double Lx, double Ly, int *counter, double *P, double *xd, double *yd, double *Rd, double Tand, double ColisionFrequency);

void Initialize_System(int Np, double To, double rho, double V,double *x, double *y, double *xo, double *yo, double *xd, double *yd,  double *vx, double *vy, double Lx, double Ly, double *fx, double *fy, double *Ep, double *Ec, double *P, double *T, double Rc2, double dt, double ParametrosPotencial[3][3]);

void Calculate_Radial_Distribution(int Np, double *g, double *x, double *y, double Lx, double Ly, double delg, double Rc, double Rc2);

//void Normalize_Radial_Distribution

double gauss (double mu,double sigma);

void Simulation(int Np, double rho, double To, double Tand, double ColisionFrequency, double t, double dt, double Transient, double CycleTime, double TotalTime, double RecordingTimeInterval, double Nbins, double Samples_gr, double gr_dt, double toAverageParameters, double StartDiffusion, double ParametrosPotencial[3][3], double ParameterRegister[6][Nc], int i);

char* GetFileName(char FileType[100],double rho, double T, int order);

void Calculate_MSD(int Np, double *xd, double *yd, double *MSD, double *vx, double *vy, double dt);

void main(){

  FILE *pd;

  pd=fopen("DiagramaDeFase.dat","w");
  fprintf(pd,"Np rho Tandersen Ep D T P\n");
  fclose(pd);

  
  int Np = 800; // Número de Particulas
  double rho = .242;
  double To = 0.15;//Temperatura Inicial
  double Tand = 0.15;
  double Tkeep = Tand;//Temperatura do Termostato
  double rhokeep = rho;//Densidade para reinicializar o loop
  double ColisionFrequency = 0.01; //Freauência de Colisões com o Reservatório
  double t = 0.0;
  double dt = 0.01;
  double Transient = 0.0;
  double CycleTime = 0.0; //Sera utilizado para o diagrama de fases
  double TotalTime = 15000.0;
  double RecordingTimeInterval = .5; //Parâmetros de tempo
  int Nbins = 100; //bins da g(r) "Nhist"
  int Samples_gr=200; //Medidas de g(r)
  double gr_dt = dt; //Intervalo de tempo entre cada medida da g(r)
  double toAverageParameters = 13000.0;
  double StartDiffusion = 5000.0;

  double ParameterRegister[6][Nc]; //[D,Ep,E,T,P,rho]
  

  double h1 = 2.1895, c1 = 0.8199, w1 = 0.0420;
  double h2 = 9.6240, c2 = 0.7947, w2 = 0.7197;
  double h3 = -3.8685, c3 = 1.1664, w3 = 0.240; //Parâmetros do Potencial
  double ParametrosPotencial[3][3];
  ParametrosPotencial[0][0]=h1;
  ParametrosPotencial[0][1]=c1;
  ParametrosPotencial[0][2]=w1;
  ParametrosPotencial[1][0]=h2;
  ParametrosPotencial[1][1]=c2;
  ParametrosPotencial[1][2]=w2;
  ParametrosPotencial[2][0]=h3;
  ParametrosPotencial[2][1]=c3;
  ParametrosPotencial[2][2]=w3;

  printf("h1 = %f H2 = %f H3 = %f C1= %f  C2 = %f  c3 = %f  w1 = %f  w2 = %f  w3 = %f \n",h1,h2,h3,c1,c2,c3,w1,w2,w3);

  double rhoVector[24];
  int i;

  rhoVector[0]=0.25;
  rhoVector[1]=0.2625;
  rhoVector[2]=0.275;
  rhoVector[3]=0.3;
  rhoVector[4]=0.3125;
  rhoVector[5]=0.325;
  rhoVector[6]=0.35;
  rhoVector[7]=0.3625;
  rhoVector[8]=0.375;
  rhoVector[9]=0.38125;
  rhoVector[10]=0.3875;
  rhoVector[11]=0.39375;
  rhoVector[12]=0.4;
  rhoVector[13]=0.405;
  rhoVector[14]=0.41;
  rhoVector[15]=0.415;
  rhoVector[16]=0.42;
  rhoVector[17]=0.425;
  rhoVector[18]=0.45;
  rhoVector[19]=0.475;
  rhoVector[20]=0.5;
  rhoVector[21]=0.175;
  rhoVector[22]=0.2;
  rhoVector[23]=0.225;
  for(i=0;i<Nc;i++){
    
    pd=fopen("DiagramaDeFase.dat","a");
    printf("%d\n",i);

    //Simulação
    Simulation(Np,rhoVector[i],Tand,Tand,ColisionFrequency,0.0,dt,Transient,CycleTime,TotalTime,RecordingTimeInterval,Nbins,Samples_gr,gr_dt,toAverageParameters,StartDiffusion,ParametrosPotencial,ParameterRegister,i);

    //Imprimindo parâmetros no arquivo de saída
    fprintf(pd,"%d %lf %lf %lf %lf %lf %lf\n",Np,ParameterRegister[5][i],ParameterRegister[0][i],ParameterRegister[1][i],ParameterRegister[2][i],ParameterRegister[3][i],ParameterRegister[4][i]);

    
    //printf("%d rho = %lf T = %lf  P = %lf\n",i,ParameterRegister[5][i],ParameterRegister[3][i],ParameterRegister[4][i]);
    
    //incremento em parâmetro 1
    //rho=rho+0.025;
    //rho=rho+0.05;
    //loop
    /*if(i%25==24){
      //incremento no parâmetro 2 e reset do parâmetro 1
      rho=rhokeep;
      Tand=Tand+0.05;
      To=Tand;
      }*/

 fclose(pd);
  }//Fim do loop das Simulações


fclose(pd);

}

void Simulation(int Np, double rho, double To, double Tand, double ColisionFrequency, double t, double dt, double Transient, double CycleTime, double TotalTime, double RecordingTimeInterval, double Nbins, double Samples_gr, double gr_dt, double toAverageParameters, double StartDiffusion, double ParametrosPotencial[3][3], double ParameterRegister[6][Nc], int i){

  FILE *oup,*fin,*dif,*rd;
  char *oupname,*finname,*rdname,*difname;

  //oupname=GetFileName("Output",rho,Tand,i);
  //oup=fopen(oupname,"w");
  oup=fopen("Output.dat","w");
  finname=GetFileName("FinalParameters",rho,Tand,i);
  fin=fopen(finname,"w");

  difname=GetFileName("DiffusionCoefficient",rho,Tand,i);
  dif=fopen(difname,"w");

  /*fprintf(fin,"#Np = %d\n#Lx = %lf\n#Ly = %lf\n#BinWidth = %d\n#Rc = %lf\n",Np,Lx,Ly,BinWidth,Rc);
  fprintf(fin,"x y vx vy\n");
  fprintf(dif,"t MSD D\n");*/


  
  double V = Np*1.0/rho;
  double L = pow(V,(1.0/2.0));
  double Lx=L, Ly=L;
  double Rc=5.0; //Raio de Corte
  unsigned int ii,j,step,q,dd;
  double *xo,*yo,*x,*y,*vx,*vy,*fx,*fy; //Variaveis de movimento
  double *xd,*yd,*xdo,*ydo,*MSD,msd; //Variaveis da difusao
  double D=0.0;
  double *Ep,*Ec,*P,*T,*E; //Parâmetros de medida
  double Epav=0., Epstd=0.0, Ecav=0.0, Ecstd=0.0, Pav=0.0, Pstd=0.0, Tav=0.0, Tstd=0.0, Eav=0.0, Estd=0.0, Dav=0.0, Dstd=0.0, avCounter=0.0;
  double *g; //Parâmetros da g(r)
  int MeasurementStepInterval = RecordingTimeInterval/dt;
  int MeasurementCounter = 0;
  int Nsteps = TotalTime/dt;
  int Samples = Transient/dt + (TotalTime - Transient)/RecordingTimeInterval; //Paraetros de tempo
  int *ColisionCounter; //Contagem de colisões com o reservatório térmico 
  double Rc2=Rc*Rc; //Raio de corte ao quadrado
  double BinWidth=L/(2*Nbins);
  double gr_to=TotalTime - (Samples_gr - 1)*gr_dt;
  int grSamplesCounter;
  double gr_t=0.0;
  double rb,vol,pi=3.1415;
  double toup=25000.0;


  size_t size = Np*sizeof(double);

  x= (double*)malloc(size);
  y= (double*)malloc(size);
  // z= (double*)malloc(size);
  xo=(double*)malloc(size);
  yo=(double*)malloc(size);
  //zo=(double*)malloc(size);
  xd=(double*)malloc(size);
  yd=(double*)malloc(size);
  //zd=(double*)malloc(size);
  MSD=(double*)malloc(size);
  vx=(double*)malloc(size);
  vy=(double*)malloc(size);
  //vz=(double*)malloc(size);
  fx=(double*)malloc(size);
  fy=(double*)malloc(size);
  // fz=(double*)malloc(size);

  Ep=(double*)malloc(sizeof(double));
  Ec=(double*)malloc(sizeof(double));
  ColisionCounter=(int*)malloc(sizeof(int));
  //grSamplesCounter=(int*)malloc(sizeof(int));
  T=(double*)malloc(sizeof(double));
  g= (double*)malloc(sizeof(double)*Nbins);
  P= (double*)malloc(sizeof(double));


  
  fprintf(fin,"#Np = %d\n#Lx = %lf\n#Ly = %lf\n#BinWidth = %lf\n#Rc = %lf\n#dt = %lf\n#rho = %lf\n#Tand = %lf\n",Np,Lx,Ly,BinWidth,Rc,dt,rho,Tand);
  fprintf(fin,"x y vx vy\n");
  fprintf(dif,"t MSD D\n");
  
  for(ii=0;ii<Nbins;ii++){
    g[ii]=0;
  }

  //printf("Lx = %lf Ly =  %lf\n",Lx,Ly);
  Initialize_System(Np,To,rho,V,x,y,xo,yo,xd,yd,vx,vy,Lx,Ly,fx,fy,Ep,Ec,P,T,Rc2,dt,ParametrosPotencial);
   //printf("Sistema inicializado  ");
  for(step=0;step<Nsteps;step++){
    //printf("time = %f\n",t);
    
    Integrate_Velocity_Verlet_pt1(Np,x,y,vx,vy,fx,fy,dt,Ec,Ep,Rc2,Lx,Ly,ColisionCounter,P,xd,yd,MSD);
    if(t>=StartDiffusion){
      Calculate_MSD(Np,xd,yd,MSD,vx,vy,dt);
    }
    Compute_Forces(Np,V,rho,T,fx,fy,x,y,Ep,Rc2,Lx,Ly,P,ParametrosPotencial);
    Integrate_Velocity_Verlet_pt2(Np,x,y,vx,vy,fx,fy,dt,Ec,Ep,T,Rc2,Lx,Ly,ColisionCounter,P,xd,yd,MSD,Tand,ColisionFrequency);
    //fprintf(oup,"%lf %lf %lf\n",t,*T,*P);

    if(t>=StartDiffusion){
        D=0.0;
      //Calculate_MSD(Np,xd,yd,MSD,vx,vy,dt);
      for(dd=0;dd<Np;dd++){
	D=D+MSD[dd];
      }
      msd=D/(Np);
      D=D/(4*(t-StartDiffusion)*Np);
      if(D>10){
	D=0;
      }
      //fprintf(dif,"%lf %lf %lf\n",t,msd,D);
      if(t>=toAverageParameters){
	
	fprintf(dif,"%lf %lf %lf\n",t,msd,D);
	//if(t>=toup){
	//fprintf(oup,"%lf %lf %lf\n",t,*T,*P);
	//}
	Epav=Epav+*Ep;
	Ecav=Ecav+*Ec;
	Eav=Eav+(*Ec)+(*Ep);
	Pav=Pav+*P;
	Tav=Tav+*T;
	avCounter++;
	Dav=Dav+D;
      }

      
    }
    



    t=t+dt;
  }



  Epav=Epav/(avCounter*Np);
  Ecav=Ecav/(avCounter*Np);
  Tav=Tav/(avCounter);
  Pav=Pav/(avCounter);
  Eav=Eav/(avCounter*Np);
  Dav=Dav/avCounter;

  ParameterRegister[0][i]=Tand;
  ParameterRegister[1][i]=Epav;
  ParameterRegister[2][i]=Dav;
  ParameterRegister[3][i]=Tav;
  ParameterRegister[4][i]=Pav;
  ParameterRegister[5][i]=rho;

  //printf("Imprimindo no arquivo PD  ");
  //fprintf(pd,"%d %lf %lf %lf %lf %lf %lf\n",Np,rho,Ecav,Epav,Eav,Tav,Pav);
  /* for(ii=0;ii<Nbins;ii++){
    rb=(pow(((ii+1)*BinWidth),2)-pow((ii*BinWidth),2));
    vol=pi*rb*rho;
    g[ii]=g[ii]/(grSamplesCounter*vol*Np);
    fprintf(rd,"%f %f\n",(ii+0.5)*(BinWidth),g[ii]);
  }
  */
  for(i=0;i<Np;i++){
    fprintf(fin,"%lf %lf %lf %lf\n",x[i],y[i],vx[i],vy[i]);
    //printf(vf,"%lf %lf\n",vx[i],vy[i]);
    }

  //fclose(rd);
  fclose(dif);
  //fclose(rf);
  fclose(fin);
  //fclose(rdist);
  fclose(oup);
  
  /* free(x); free(y); free(z); printf("free(x,y,z)  ");
 free(vx); free(vy); free(vz);  printf("free(vx,vy,vz)  ");
 free(fx); free(fy); free(fz);  printf("free(fx,fy,fz)  ");
 free(xo); free(yo); free(zo);  printf("free(xo,yo,zo)  ");
 free(xd); free(yd); free(zd); free(MSD);  printf("free(xd,yd,zd,MSD)  ");
 free(g);  printf("free(g)  ");
 free(Ep); free(Ec); free(P); free(T); printf("free(Ep,Ec,P,T)  ");
 free(ColisionCounter); free(grSamplesCounter);  printf("free(ColisionCounter,grSamplesCounter)  ");*/




}



void Initialize_Positions(int Np, double V, double *x,double *y,double *xo,double *yo, double *xd, double *yd, double *vx,double *vy,  double Lx, double Ly, double dt){

        unsigned int i;
        int n,N;
	double l;
	//FILE *ro;

	//ro=fopen(filename,"w");
	//fprintf(ro,"xo yo\n");

	n=ceil(pow(Np,1./2.));
	N=pow(n,2);
	l=pow(V/N,(1./2.));
	//printf("%d  l=%f\n",n,l);

        for(i=0;i<Np;i++)
        {
                x[i]=(i%n)*l;
                y[i]=(((i/n)%n))*l;
		//z[i]=(((i/n)/n)%n)*l;
		
		
                xo[i]=x[i]-vx[i]*dt;
                yo[i]=y[i]-vy[i]*dt;
		//zo[i]=z[i]-vz[i]*dt;

		xo[i]=fmod(xo[i],Lx);
		yo[i]=fmod(yo[i],Ly);
		//zo[i]=fmod(zo[i],Lz);
		xo[i]=xo[i]-rint(xo[i]/Lx)*Lx;
		yo[i]=yo[i]-rint(yo[i]/Ly)*Ly;
		//zo[i]=zo[i]-rint(zo[i]/Lz)*Lz;

		x[i]=fmod(x[i],Lx);
		y[i]=fmod(y[i],Ly);
		//z[i]=fmod(z[i],Lz);
		x[i]=x[i]-rint(x[i]/Lx)*Lx;
		y[i]=y[i]-rint(y[i]/Ly)*Ly;
		//z[i]=z[i]-rint(z[i]/Lz)*Lz;

		xd[i]=0.0;
		yd[i]=0.0;

		//	fprintf(ro,"%lf %lf \n",xo[i],yo[i]);

        }
	//fclose(ro);

}

void Initialize_Velocities(int Np, double dt, double To, double *vx,double *vy, double *T, double *Ec){

  unsigned i;
  double Vcx=0.,Vcy=0.,V2=0.,fs,V20=0.,fs0=0.;
  //FILE *vo;

  // vo=fopen(filename,"w");
  //fprintf(vo,"vxo vyo\n");

  for(i=0;i<Np;i++)
  {
    vx[i]=(1.*rand()/RAND_MAX);
    vy[i]=(1.*rand()/RAND_MAX);
    //vz[i]=(1.*rand()/RAND_MAX);

      Vcx+=vx[i];
      Vcy+=vy[i];
      //Vcz+=vz[i];
      V20+=vx[i]*vx[i] + vy[i]*vy[i];
  }

  Vcx/=Np;
  Vcy/=Np;
  // Vcz/=Np;
  V20/=Np;
  fs0=sqrt(2.*To/V20);

  for(i=0;i<Np;i++)
  {
        vx[i]=(vx[i]-Vcx);
        vy[i]=(vy[i]-Vcy);
	//vz[i]=(vz[i]-Vcz);
  }

   for(i=0;i<Np;i++){

     V2 +=vx[i]*vx[i] + vy[i]*vy[i];
   }

   V2/=Np;
   fs=sqrt(2.*To/V2);
   *Ec=0.0;

   for(i=0;i<Np;i++)
  {
        vx[i]=(vx[i])*fs;
        vy[i]=(vy[i])*fs;
	//vz[i]=(vz[i])*fs;
	*Ec=*Ec+vx[i]*vx[i]/2.0+vy[i]*vy[i]/2.0;
	*T=*Ec*2.0/(2.0*Np);
	//fprintf(vo,"%lf %lf %lf\n",vx[i],vy[i],sqrt(vx[i]*vx[i]+vy[i]*vy[i]));
	
	}

   //fclose(vo);
   // printf("V20= %f e fs0= %f     V2=%f  e  fs=%f\n",V20,fs0,V2,fs);
}


void Compute_Forces(int Np, double V, double rho, double *T, double *fx,double *fy,double *x,double *y,double *Ep, double Rc2, double Lx, double Ly, double *P, double ParametrosPotencial[3][3])
{
  unsigned i,j,a=0;
  double  x1,y1,x2,y2,dx,dy,r2,r2i,r6,ff,Energy=0.,Ecut,pr,rr,rri,Rc,h1,h2,h3,c1,c2,c3,w1,w2,w3,InstantaneousEnergy;
  //FILE *checkforces;

  h1=ParametrosPotencial[0][0]; c1=ParametrosPotencial[0][1]; w1=ParametrosPotencial[0][2];
  h2=ParametrosPotencial[1][0]; c2=ParametrosPotencial[1][1]; w2=ParametrosPotencial[1][2];
  h3=ParametrosPotencial[2][0]; c3=ParametrosPotencial[2][1]; w3=ParametrosPotencial[2][2];

  //printf("h1 = %f H2 = %f H3 = %f C1= %f  C2 = %f  c3 = %f  w1 = %f  w2 = %f  w3 = %f \n",h1,h2,h3,c1,c2,c3,w1,w2,w3);
  Rc=sqrt(Rc2);
  //checkforces=fopen("CheckForces.dat","a");

  pr=0.;
  Ecut=4*(pow(Rc2,-6)-pow(Rc2,-3)) + h1*exp(-((Rc-c1)/w1)*((Rc-c1)/w1)) + h2*exp(-((Rc-c2)/w2)*((Rc-c2)/w2)) + h3*exp(-((Rc-c3)/w3)*((Rc-c3)/w3)) ;
  /*Zerar forças antes de cada iteração*/
  for(i=0;i<Np;i++)
  {
        fx[i] = 0.;
        fy[i] = 0.;
	//fz[i] = 0.;
  }

   for(i=0;i<Np;i++)
  {
    x1=x[i]; y1=y[i];

        for(j=i+1;j<Np;j++)
        {
                x2=x[j];

                dx=x1 - x2;
                dx=fmod(dx,Lx);
	        dx=dx-rint(dx/Lx)*Lx;

                y2=y[j];

                dy=y1 - y2;
                dy=fmod(dy,Ly);
                dy=dy-rint(dy/Ly)*Ly;

		
                r2=dx*dx + dy*dy;
		rr=pow(r2,1./2.);


                if(r2<Rc2)
                {
		  rri=1/rr;
		  r2i=1/r2 ;
		  r6=r2i*r2i*r2i;
		  ff=48*r2i*r6*(r6-0.5)+rri*h1*exp(-((rr-c1)/w1)*((rr-c1)/w1))*2*((rr-c1)/(w1*w1))+rri*h2*exp(-((rr-c2)/w2)*((rr-c2)/w2))*2*((rr-c2)/(w2*w2))+rri*h3*exp(-((rr-c3)/w3)*((rr-c3)/w3))*2*((rr-c3)/(w3*w3));
		  pr=pr+ff*dx*dx+ff*dy*dy;
		  fx[i] += ff*dx;   fy[i] += ff*dy;
		  fx[j] -= ff*dx;   fy[j] -= ff*dy;
		  //InstantaneousEnergy=4*r6*(r6-1) + h1*exp(-((rr-c1)/w1)*((rr-c1)/w1)) + h2*exp(-((rr-c2)/w2)*((rr-c2)/w2)) + h3*exp(-((rr-c3)/w3)*((rr-c3)/w3)) - Ecut;
		  Energy += 4*r6*(r6-1) + h1*exp(-((rr-c1)/w1)*((rr-c1)/w1)) + h2*exp(-((rr-c2)/w2)*((rr-c2)/w2)) + h3*exp(-((rr-c3)/w3)*((rr-c3)/w3)) - Ecut;
		  //fprintf(checkforces,"%lf %lf %lf\n",rr,ff,InstantaneousEnergy);
            }
        }
   }
  *Ep=Energy;
  *P=pr*1./(2.*V)+rho*(*T);
  //fclose(checkforces);
}



void Integrate_Velocity_Verlet_pt1(int Np, double *x, double *y, double *vx, double*vy, double *fx, double *fy, double dt, double *Ec, double *Ep, double Rc2, double Lx, double Ly, int *counter, double *P, double *xd, double *yd, double *MSD){

  int i;
  double ec,temp,sigma,rnd,ccc=0.,q,xkeep,ykeep;
  //FILE *checkbc;

  //checkbc=fopen("CheckBoundaryConditions.dat","a");

  ec=0.0;

  for(i=0;i<Np;i++){

    x[i]=x[i]+vx[i]*dt+0.5*fx[i]*dt*dt;
    y[i]=y[i]+vy[i]*dt+0.5*fy[i]*dt*dt;
    //z[i]=z[i]+vz[i]*dt+0.5*fz[i]*dt*dt;

    //Diffusion
    //xd[i]=xd[i]+vx[i]*dt;
    //yd[i]=yd[i]+vy[i]*dt;
    //zd[i]=zd[i]+vz[i]*dt;
    //MSD[i]=xd[i]*xd[i]+yd[i]*yd[i];

    //xkeep=x[i]; ykeep=y[i];
    x[i]=fmod(x[i],Lx);
    y[i]=fmod(y[i],Ly);
    //z[i]=fmod(z[i],Lz);
    x[i]=x[i]-rint(x[i]/Lx)*Lx;
    y[i]=y[i]-rint(y[i]/Ly)*Ly;
    //z[i]=z[i]-rint(z[i]/Lz)*Lz;
    
    /* if(x[i]!=xkeep){
      fprintf(checkbc,"%d %lf %lf %lf X\n",i,xkeep,fmod(xkeep,Lx),x[i]);
    }

    if(y[i]!=ykeep){
      fprintf(checkbc,"%d %lf %lf %lf Y\n",i,ykeep,fmod(ykeep,Ly),y[i]);
      }*/
    
    vx[i]=vx[i]+0.5*fx[i]*dt;
    vy[i]=vy[i]+0.5*fy[i]*dt;
    //vz[i]=vz[i]+0.5*fz[i]*dt;

  }
  //fclose(checkbc);
}

// compute_forces(fx,fy,x,y,Ep,Rc2,Lx,Ly,P);

void Integrate_Velocity_Verlet_pt2(int Np, double *x, double *y, double *vx, double *vy, double *fx, double *fy, double dt, double *Ec, double *Ep, double *T, double Rc2, double Lx, double Ly, int *counter, double *P, double *xd, double *yd, double *Rd, double Tand, double ColisionFrequency){

  int i;
  double ec=0.0,q,rnd,ccc=0.0,temp,sigma;
  for(i=0;i<Np;i++){
    vx[i]=vx[i]+0.5*fx[i]*dt;
    vy[i]=vy[i]+0.5*fy[i]*dt;
    //vz[i]=vz[i]+0.5*fz[i]*dt;

    ec=ec+(vx[i]*vx[i]+vy[i]*vy[i])*0.5;

  }

   //Andersen
  //temp=ec*2./(2.*Np);
     sigma=sqrt(Tand);

   for(i=0;i<Np;i++){
     rnd=(1.*rand()/RAND_MAX);
     if(rnd<ColisionFrequency){
       q=gauss(0.,sigma);
       vx[i]=gauss(0.,sigma);
       vy[i]=gauss(0.,sigma);
       //ccc=ccc+1.;
     }


   }


   //*counter=ccc;
   *Ec=ec;
   *T=ec/Np;
}

void Calculate_MSD(int Np, double *xd, double *yd, double *MSD, double *vx, double *vy, double dt){

  unsigned int i;
  FILE *arq;

  arq=fopen("DiffRegister.dat","a");

  for(i=0;i<Np;i++){

    xd[i]=xd[i]+vx[i]*dt;
    yd[i]=yd[i]+vy[i]*dt;
    MSD[i]=xd[i]*xd[i]+yd[i]*yd[i];
    //fprintf(arq,"%lf\n",MSD[i]);
  }

  fclose(arq);
}

void Initialize_System(int Np, double To, double rho, double V,double *x, double *y, double *xo, double *yo, double *xd, double *yd, double *vx, double *vy, double Lx, double Ly, double *fx, double *fy, double *Ep, double *Ec, double *P, double *T, double Rc2, double dt, double ParametrosPotencial[3][3]){

  //Initialize_Positions(Np,V,x,y,z,xo,yo,zo,vx,vy,vz,Lx,Ly,Lz,dt);
  FILE *oup,*rf,*vf,*rdist,*dif,*rd,*pd;
  oup=fopen("Output.dat","w");
  rf=fopen("FinalPositions.dat","w");
  vf=fopen("FinalVelocities.dat","w");
  rd=fopen("RadialDistribution.dat","w");
  dif=fopen("DiffusionCoefficient.dat","w");

  fprintf(pd,"Np rho To Ep D T P \n");
  fclose(pd);
  fprintf(rf,"x y\n");
  fclose(rf);
  fprintf(vf,"vx vy v2\n");
  fclose(vf);
  fprintf(rd,"r g(r)\n");
  fclose(rd);
  fprintf(dif,"t MSD D\n");
  fclose(dif);

  Initialize_Velocities(Np,dt,To,vx,vy,T,Ec);

  Initialize_Positions(Np,V,x,y,xo,yo,xd,yd,vx,vy,Lx,Ly,dt);
  
  Compute_Forces(Np,V,rho,T,fx,fy,x,y,Ep,Rc2,Lx,Ly,P,ParametrosPotencial);
}

void Calculate_Radial_Distribution(int Np, double *g, double *x, double *y, double Lx, double Ly, double delg, double Rc, double Rc2){

  unsigned i,j,ig,iig;
  double dx,dy,r,r2;

  //grSamplesCounter=grSamplesCounter+1;
  
  for(i=0;i<Np;i++){
    for(j=i+1;j<Np;j++){
      dx=x[i]-x[j];
      dy=y[i]-y[j];
      //dz=z[i]-z[j];
      dx=dx-rint(dx/Lx)*Lx;
      dy=dy-rint(dy/Ly)*Ly;
      //dz=dz-rint(dz/Lz)*Lz;

      r2=(dx*dx+dy*dy);
      // if(r2<=Rc2){
	
	r=sqrt(r2);
	ig=rint(r/delg);
	iig=rint(ig);
	g[iig]=g[iig]+2;



	//}
    }
  }

  //return(grSamplesCounter);

}


double gauss (double mu,double sigma){

  double v1,v2,r,l;
  do{

    v1 = 2.0 * rand() / (double)RAND_MAX - 1.0;
    v2 = 2.0 * rand() / (double)RAND_MAX - 1.0;
    r=v1*v1+v2*v2;
  }while(r>=1.);

  l=v1*sqrt(-2.*log(r)/r);
  l=mu+sigma*l;
  return l;
}

char* GetFileName(char FileType[100],double rho, double T, int order){

  static char filename[100],strrho[50],strT[50],ord[10];

  strcpy(filename,FileType);
  sprintf(ord,"%d",order);
  strcat(filename,ord);

  if(rho<1){
    if(rho<0.1){
      if(rho<0.01){
	sprintf(strrho,"__rho_0_00%d",(int)rint(rho*10000));
      }
      else{
	sprintf(strrho,"__rho_0_0%d",(int)rint(rho*1000));
      }
    }
    else{
      sprintf(strrho,"__rho_0_%d",(int)rint(rho*100));
    }

  }
  else{

    sprintf(strrho,"__rho_%d_%d",(int)floor(rho),(int)rint((rho-floor(rho))*100.0));

  }

  if(T<1){

    if(T<0.1){
      if(T<0.01){
	sprintf(strT,"__T_0_00%d.dat",(int)rint(T*10000));
      }
      else{
	sprintf(strT,"__T_0_0%d.dat",(int)rint(T*1000));
      }
    }
    else{//T>0.1
      sprintf(strT,"__T_0_%d.dat",(int)rint(T*100));
    }
  }  
  else{
    sprintf(strT,"__T_%d_%d.dat",(int)floor(T),(int)rint((T-floor(T))*100.0));
  }
  

  strcat(filename,strrho);
  strcat(filename,strT);
  return(filename);

}
