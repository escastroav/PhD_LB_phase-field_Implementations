#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1;
const int Ly=64;

const int Q=9;

const double tau=1.4;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
const double U_Cs2 = 3.0;
const double Cs2 = 1.0/U_Cs2;
const double TresUmUsobre2tau=U_Cs2*(1-1.0/(2*tau));
const double W0 = 10.0;

class LB{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
  double g[Lx][Ly][Q], gnew[Lx][Ly][Q]; // f[ix][iy][i]
public:
  LB(void);
  double phi(int ix,int iy,bool UseNew);
  double rho(int ix,int iy,bool UseNew);
  double dx_phi(int ix,int iy,bool UseNew);
  double dy_phi(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew,double Fx);
  double Jy(int ix,int iy,bool UseNew,double Fy);
  double Fi(double Ux0,double Uy0,double Fx,double Fy,int i);
  double Gi(double phi0,double dxphi,double dyphi,int i);
  double feq(double rho0,double Ux0,double Uy0,int i);
  double geq(double phi0,double Ux0,double Uy0,int i);
  void Collision(double gx,double gy);
  void Advection(void);
  void Init(double phi0,double rho0,double Ux0,double Uy0);
  //void ImposeFields(void);
  void Print(const char * NombreArchivo,double gx,double gy);
};
LB::LB(void){
  //Cargar los pesos
  w[0]=4/9.0; 
  w[1]=w[2]=w[3]=w[4]=1/9.0;
  w[5]=w[6]=w[7]=w[8]=1/36.0;
  //Cargar los vectores
  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;

  V[0][5]=1;  V[0][6]=-1; V[0][7]=-1; V[0][8]=1;
  V[1][5]=1;  V[1][6]=1;  V[1][7]=-1; V[1][8]=-1;
}
double LB::rho(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]; else suma+=f[ix][iy][i];
  return suma;
}
double LB::phi(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=gnew[ix][iy][i]; else suma+=g[ix][iy][i];
  return suma;
}
double LB::dx_phi(int ix,int iy,bool UseNew){
  int i; double sum;
  int jx, jy;
  if(iy==0||iy==Ly-1)
    return 0;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    jy=(Ly+iy+V[1][i])%Ly;
    sum+=w[i]*V[0][i]*phi(jx,jy,UseNew);
  }
  return sum*U_Cs2;
}
double LB::dy_phi(int ix,int iy,bool UseNew){
  int i; double sum;
  int jx, jy;
  if(iy==0||iy==Ly-1)
    return 0;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    jy=(Ly+iy+V[1][i])%Ly;
    sum+=w[i]*V[1][i]*phi(jx,jy,UseNew);
  }
  return sum*U_Cs2;
}
double LB::Jx(int ix,int iy,bool UseNew,double Fx){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[0][i]; else suma+=f[ix][iy][i]*V[0][i];
  return suma+0.5*Fx;
}
double LB::Jy(int ix,int iy,bool UseNew,double Fy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[1][i]; else suma+=f[ix][iy][i]*V[1][i];
  return suma+0.5*Fy;
}
double LB::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2);
}
double LB::geq(double phi0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  //return phi0*w[i]*(1+U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2);
  return phi0*w[i]*(1+U_Cs2*UdotVi);
}
double LB::Fi(double Ux0,double Uy0,double Fx,double Fy,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], FdotVi=Fx*V[0][i]+Fy*V[1][i], UdotF=Ux0*Fx+Uy0*Fy;
  return TresUmUsobre2tau*w[i]*(FdotVi-UdotF+3*UdotVi*FdotVi);
}
double LB::Gi(double phi0,double dxphi,double dyphi,int i){
  if(dxphi==0&&dyphi==0)
    return 0;
  double grad = sqrt(dxphi*dxphi + dyphi*dyphi);
  double nx = dxphi / grad, ny = dyphi / grad;
  double NdotVi = nx*V[0][i]+ny*V[1][i];
  double lambda = 4.0*phi0*(1.0 - phi0)/W0;
  return TresUmUsobre2tau*w[i]*lambda*NdotVi;
}
void LB::Collision(double gx,double gy){
  int ix,iy,i; double phi0=0,rho0=0,Ux0=0,Uy0=0; double Fx,Fy;
  double dxphi0,dyphi0;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macrosc�picas
      rho0=rho(ix,iy,false);  
      phi0=phi(ix,iy,false);  
      dxphi0=dx_phi(ix,iy,false);
      dyphi0=dy_phi(ix,iy,false);
      Fx=gx*rho0;   Fy=gy*rho0; 
      Ux0=Jx(ix,iy,false,Fx)/rho0;  Uy0=Jy(ix,iy,false,Fy)/rho0;
      for(i=0;i<Q;i++){
	      fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Ux0,Uy0,i)+Fi(Ux0,Uy0,Fx,Fy,i);
	      gnew[ix][iy][i]=UmUtau*g[ix][iy][i]+Utau*geq(phi0,Ux0,Uy0,i)+Gi(phi0,dxphi0,dyphi0,i);}
    }
}
void LB::Advection(void){
  double D = 1.0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++){
	      f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
	      g[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=gnew[ix][iy][i];
        g[ix][0][0]=D*gnew[ix][0][0];
        g[ix][0][1]=D*gnew[ix][0][1];
        g[ix][0][2]=D*gnew[ix][0][2];
        g[ix][0][3]=D*gnew[ix][0][3];
        g[ix][0][4]=D*gnew[ix][0][4];
        g[ix][0][5]=D*gnew[ix][0][5];
        g[ix][0][6]=D*gnew[ix][0][6];
        g[ix][0][7]=D*gnew[ix][0][7];
        g[ix][0][8]=D*gnew[ix][0][8];
        g[ix][Ly-1][0]=D*gnew[ix][Ly-1][0];
        g[ix][Ly-1][1]=D*gnew[ix][Ly-1][1];
        g[ix][Ly-1][2]=D*gnew[ix][Ly-1][2];
        g[ix][Ly-1][3]=D*gnew[ix][Ly-1][3];
        g[ix][Ly-1][4]=D*gnew[ix][Ly-1][4];
        g[ix][Ly-1][5]=D*gnew[ix][Ly-1][5];
        g[ix][Ly-1][6]=D*gnew[ix][Ly-1][6];
        g[ix][Ly-1][7]=D*gnew[ix][Ly-1][7];
        g[ix][Ly-1][8]=D*gnew[ix][Ly-1][8];
      
      }
}
void LB::Init(double phi0,double rho0,double Ux0,double Uy0){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++){
        //phi0 = exp(-(iy-Ly/2)*(iy-Ly/2)/50.0);
        phi0 = 0.5-0.5*tanh((iy-Ly/2)*2.0);
	      f[ix][iy][i]=feq(rho0,Ux0,Uy0,i);
	      g[ix][iy][i]=geq(phi0,Ux0,Uy0,i);}
}
/**void LB::ImposeFields(void){
  int i,ix,iy; double rho0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      //Paredes
      if(iy==0 || iy==Ly-1) 
	for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,0,0,i);
    }
}**/
void LB::Print(const char * NombreArchivo,double gx,double gy){
  ofstream MiArchivo(NombreArchivo); double phi0,rho0,Ux0,Uy0; double Fx,Fy;
  int ix=0;
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true); 
      phi0=phi(ix,iy,true);  
      Fx=gx*rho0;   Fy=gy*rho0; 
      Ux0=Jx(ix,iy,true,Fx)/rho0;  Uy0=Jy(ix,iy,true,Fy)/rho0;
      MiArchivo<<iy<<" "<<phi0<<" "<<rho0<<endl;
    }
  MiArchivo.close();
}


int main(void){
  LB Aire;
  int t,tmax=1000;
  double phi0=1.0, rho0=1.0, g=0.0;
  
  Aire.Init(phi0,rho0,0,0);
  
  for(t=0;t<tmax;t++){
    Aire.Collision(g,0);
    //Aire.ImposeFields();
    Aire.Advection();
  }
  
  Aire.Print("AD_FL1000.dat",g,0);

  return 0;
}


