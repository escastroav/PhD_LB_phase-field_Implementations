#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1;
const int Ly=100;

const int Q=9;

const double tau=1.2;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
const double U_Cs2 = 3.0;
const double Cs2 = 1.0/U_Cs2;
const double TresUmUsobre2tau=U_Cs2*(1-1.0/(2*tau));

class LB{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
  double rho_g = 1.0, rho_l = 0.5;
  double rho_0 = rho_g;
public:
  LB(void);
  double rho(int ix,int iy,bool UseNew);
  double Psi(int ix, int iy, bool UseNew);
  double Psi_rho(double rho);
  double G(double rho_l,double rho_g);
  double FSCx(int ix,int iy, bool UseNew);
  double FSCy(int ix,int iy, bool UseNew);
  double Jx(int ix,int iy,bool UseNew,double FSCx,double Fgx);
  double Jy(int ix,int iy,bool UseNew,double FSCy,double Fgy);
  double Fi(double Ux0,double Uy0,double Fx,double Fy,int i);
  double feq(double rho0,double Ux0,double Uy0,int i);
  void Collision(double gx,double gy);
  void Advection(void);
  void Init(double Ux0,double Uy0);
  void ImposeFields(void);
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
double LB::Psi(int ix, int iy, bool UseNew){
  return rho_0*(1.0 - exp(-rho(ix,iy,UseNew)/rho_0));
}
double LB::Psi_rho(double rho){
  return rho_0*(1.0 - exp(-rho/rho_0));
}
double LB::G(double rho_l,double rho_g){
  return 2*(rho_g - rho_l)/(Psi_rho(rho_l) - Psi_rho(rho_g)); 
}
double LB::FSCx(int ix,int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[0][i]*Psi(jx,jy,UseNew);}
  return -G(rho_l,rho_g)*Psi(ix,iy,UseNew)*sum;
}
double LB::FSCy(int ix,int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[1][i]*Psi(jx,jy,UseNew);}
  return -G(rho_l,rho_g)*Psi(ix,iy,UseNew)*sum;
}
double LB::Jx(int ix,int iy,bool UseNew,double FSCx,double Fgx){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[0][i]; else suma+=f[ix][iy][i]*V[0][i];
  return suma+tau*FSCx+0.5*Fgx;
}
double LB::Jy(int ix,int iy,bool UseNew,double FSCy,double Fgy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[1][i]; else suma+=f[ix][iy][i]*V[1][i];
  return suma+tau*FSCy+0.5*Fgy;
}
double LB::Fi(double Ux0,double Uy0,double Fx,double Fy,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], FdotVi=Fx*V[0][i]+Fy*V[1][i], UdotF=Ux0*Fx+Uy0*Fy;
  return TresUmUsobre2tau*w[i]*(FdotVi-UdotF+3*UdotVi*FdotVi);
}

double LB::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2);
}
void LB::Collision(double gx,double gy){
  int ix,iy,i; double rho0,Ux0,Uy0; double FSCx0,FSCy0,Fgx0,Fgy0;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macrosc�picas
      rho0=rho(ix,iy,false);  
      FSCx0=FSCx(ix,iy,false);   FSCy0=FSCy(ix,iy,false); 
      Fgx0=gx*rho0; Fgy0=gy*rho0;
      Ux0=Jx(ix,iy,false,FSCx0,Fgx0)/rho0;  Uy0=Jy(ix,iy,false,FSCy0,Fgy0)/rho0;
      for(i=0;i<Q;i++)
	      fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Ux0,Uy0,i)+Fi(Ux0,Uy0,Fgx0,Fgy0,i);
    }
}
void LB::Advection(void){
  double D = 1.0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++){
	      f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
        //bounceback iy = 0
       f[ix][0][1]=D*fnew[ix][0][3];
       f[ix][0][2]=D*fnew[ix][0][4];
       f[ix][0][3]=D*fnew[ix][0][1];
       f[ix][0][4]=D*fnew[ix][0][2];
       f[ix][0][5]=D*fnew[ix][0][7];
       f[ix][0][6]=D*fnew[ix][0][8];
       f[ix][0][7]=D*fnew[ix][0][5];
       f[ix][0][8]=D*fnew[ix][0][6];
       //bounceback iy = Ly-1
       f[ix][Ly-1][1]=D*fnew[ix][Ly-1][3];
       f[ix][Ly-1][2]=D*fnew[ix][Ly-1][4];
       f[ix][Ly-1][3]=D*fnew[ix][Ly-1][1];
       f[ix][Ly-1][4]=D*fnew[ix][Ly-1][2];
       f[ix][Ly-1][5]=D*fnew[ix][Ly-1][7];
       f[ix][Ly-1][6]=D*fnew[ix][Ly-1][8];
       f[ix][Ly-1][7]=D*fnew[ix][Ly-1][5];
       f[ix][Ly-1][8]=D*fnew[ix][Ly-1][6];
      }
}
void LB::Init(double Ux0,double Uy0){
  double W = 2.0; double rho0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      rho0 = rho_l+(rho_g-rho_l)*0.5*tanh(2*(double)(iy-Ly/4)/W)-(rho_g-rho_l)*0.5*tanh(2*(double)(iy-3*Ly/4)/W);
      for(int i=0;i<Q;i++)
	      f[ix][iy][i]=feq(rho0,Ux0,Uy0,i);
  }
}
void LB::ImposeFields(void){
  int i,ix,iy; double rho0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      //Paredes
      if(iy==0 || iy==Ly-1) 
	for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,0,0,i);
    }
}
void LB::Print(const char * NombreArchivo,double gx,double gy){
  ofstream MiArchivo(NombreArchivo); double rho0,Ux0,Uy0; 
  double FSCx0,FSCy0,Fgx0,Fgy0;
  int ix=0;
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);  
      FSCx0=FSCx(ix,iy,false);   FSCy0=FSCy(ix,iy,false); 
      Fgx0=gx*rho0;   Fgy0=gy*rho0; 
      Ux0=Jx(ix,iy,true,FSCx0,Fgx0)/rho0;  Uy0=Jy(ix,iy,true,FSCy0,Fgy0)/rho0;
      MiArchivo<<iy<<" "<<rho0<<" "<<Ux0<<endl;
    }
  MiArchivo.close();
}


int main(void){
  LB Aire;
  int t,tmax=50;
  double g=0;

  Aire.Init(0,0);
  
  for(t=0;t<tmax;t++){
    Aire.Collision(g,0);
    //Aire.ImposeFields();
    Aire.Advection();
  }
  
  Aire.Print("out2tau.dat",g,0);

  return 0;
}


