#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1;
const int Ly=64;

const int Q=9;

const double tau=1.2;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
const double U_Cs2 = 3.0;
const double Cs2 = 1.0/U_Cs2;
const double TresUmUsobre2tau=U_Cs2*(1-1.0/(2*tau));

class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew,double Fx);
  double Jy(int ix,int iy,bool UseNew,double Fy);
  double Fi(double Ux0,double Uy0,double Fx,double Fy,int i);
  double feq(double rho0,double Ux0,double Uy0,int i);
  void Collision(double gx,double gy);
  void Advection(void);
  void Init(double rho0,double Ux0,double Uy0);
  void ImposeFields(void);
  void Print(const char * NombreArchivo,double gx,double gy);
};
LatticeBoltzmann::LatticeBoltzmann(void){
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
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]; else suma+=f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew,double Fx){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[0][i]; else suma+=f[ix][iy][i]*V[0][i];
  return suma+0.5*Fx;
}
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew,double Fy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[1][i]; else suma+=f[ix][iy][i]*V[1][i];
  return suma+0.5*Fy;
}
double LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2);
}
double LatticeBoltzmann::Fi(double Ux0,double Uy0,double Fx,double Fy,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], FdotVi=Fx*V[0][i]+Fy*V[1][i], UdotF=Ux0*Fx+Uy0*Fy;
  return TresUmUsobre2tau*w[i]*(FdotVi-UdotF+3*UdotVi*FdotVi);
}
void LatticeBoltzmann::Collision(double gx,double gy){
  int ix,iy,i; double rho0,Ux0,Uy0; double Fx,Fy;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macrosc�picas
      rho0=rho(ix,iy,false);  
      Fx=gx*rho0;   Fy=gy*rho0; 
      Ux0=Jx(ix,iy,false,Fx)/rho0;  Uy0=Jy(ix,iy,false,Fy)/rho0;
      for(i=0;i<Q;i++)
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Ux0,Uy0,i)
	                +Fi(Ux0,Uy0,Fx,Fy,i);
    }
}
void LatticeBoltzmann::Advection(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}
void LatticeBoltzmann::Init(double rho0,double Ux0,double Uy0){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[ix][iy][i]=feq(rho0,Ux0,Uy0,i);
}
void LatticeBoltzmann::ImposeFields(void){
  int i,ix,iy; double rho0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      //Paredes
      if(iy==0 || iy==Ly-1) 
	for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,0,0,i);
    }
}
void LatticeBoltzmann::Print(const char * NombreArchivo,double gx,double gy){
  ofstream MiArchivo(NombreArchivo); double rho0,Ux0,Uy0; double Fx,Fy;
  int ix=0;
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);  
      Fx=gx*rho0;   Fy=gy*rho0; 
      Ux0=Jx(ix,iy,true,Fx)/rho0;  Uy0=Jy(ix,iy,true,Fy)/rho0;
      MiArchivo<<iy<<" "<<Ux0<<endl;
    }
  MiArchivo.close();
}


int main(void){
  LatticeBoltzmann Aire;
  int t,tmax=100000;
  double RHOinicial=1.0, g=0.01;
  
  Aire.Init(RHOinicial,0,0);
  
  for(t=0;t<tmax;t++){
    Aire.Collision(g,0);
    Aire.ImposeFields();
    Aire.Advection();
  }
  
  Aire.Print("Aire.dat",g,0);

  return 0;
}


