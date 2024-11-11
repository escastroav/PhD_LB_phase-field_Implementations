#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1;
const int Ly=64;

const int Q=9;

const double tau=0.8;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
const double U_Cs2 = 3.0;
const double Cs2 = 1.0/U_Cs2;
const double TresUmUsobre2tau=U_Cs2*(1-1.0/(2*tau));

class LB{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i] for Navier-Stokes
  double h[Lx][Ly][Q], hnew[Lx][Ly][Q]; // f[ix][iy][i] for Conservative Alen-Cahn
  double W=3;
  double M=0;
  double tauh=M*U_Cs2+0.5;
  double Utauh=1.0/tauh;
  double UmUtauh=1-Utauh;
  double rho_l=1,rho_g=1;
public:
  LB(void);
  double rho(int ix,int iy,bool UseNew);
  double phi(int ix,int iy,bool UseNew);
  double dphi_dx(int ix, int iy, bool UseNew);
  double dphi_dy(int ix, int iy, bool UseNew);
  double norm_gr(double gr_x, double gr_y);
  double nx(double gr_x, double norm);
  double ny(double gr_y, double norm);
  double Jx(int ix,int iy,bool UseNew,double Fx);
  double Jy(int ix,int iy,bool UseNew,double Fy);
  double Fi(double Ux0,double Uy0,double Fx,double Fy,int i);
  double feq(double rho0,double Ux0,double Uy0,int i);
  double heq(double phi0,double Ux0,double Uy0,double nx0,double ny0,int i);
  void Collision(double gx,double gy);
  void Advection(void);
  void Init(double rho0,double Ux0,double Uy0);
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
double LB::phi(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=hnew[ix][iy][i]; else suma+=h[ix][iy][i];
  return suma;
}
double LB::dphi_dx(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[0][i]*phi(jx,jy,UseNew);
  }
  return sum*U_Cs2;
}
double LB::dphi_dy(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[1][i]*phi(jx,jy,UseNew);
  }
  return sum*U_Cs2;
}
double LB::nx(double gr_x, double norm){
  return gr_x / norm;
}
double LB::ny(double gr_y, double norm){
  return gr_y / norm;
}
double LB::norm_gr(double gr_x, double gr_y){
  return (gr_x == 0 && gr_y == 0) ? 1 : sqrt(gr_x*gr_x + gr_y*gr_y);
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
double LB::heq(double phi0,double Ux0,double Uy0,double nx0,double ny0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0,  NdotVi=nx0*V[0][i]+ny0*V[1][i];
  return phi0*w[i]*(1+U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2)+w[i]*M*4*phi0*(1.0-phi0)*NdotVi*U_Cs2/W;
}
double LB::Fi(double Ux0,double Uy0,double Fx,double Fy,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], FdotVi=Fx*V[0][i]+Fy*V[1][i], UdotF=Ux0*Fx+Uy0*Fy;
  return TresUmUsobre2tau*w[i]*(FdotVi-UdotF+3*UdotVi*FdotVi);
}
void LB::Collision(double gx,double gy){
  int ix,iy,i; double rho0,Ux0,Uy0; double Fx,Fy;
  double phi0,gr_x,gr_y,norm0,nx0,ny0;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macrosc�picas
      rho0=rho(ix,iy,false);  phi0=phi(ix,iy,false);
      gr_x=dphi_dx(ix,iy,false);  gr_y=dphi_dy(ix,iy,false);       
      norm0=norm_gr(gr_x,gr_y);  nx0=nx(gr_x,norm0);  ny0=ny(gr_y,norm0);
      Fx=gx*rho0;   Fy=gy*rho0; 
      Ux0=Jx(ix,iy,false,Fx)/rho0;  Uy0=Jy(ix,iy,false,Fy)/rho0;
      for(i=0;i<Q;i++)
	      fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Ux0,Uy0,i)+Fi(Ux0,Uy0,Fx,Fy,i);
	      hnew[ix][iy][i]=UmUtau*h[ix][iy][i]+Utau*heq(phi0,Ux0,Uy0,nx0,ny0,i);
    }
}
void LB::Advection(void){
  double D = 1.0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++){
	      f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
	      h[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=hnew[ix][iy][i];
        //bounceback iy = 0
       f[ix][0][1]=D*fnew[ix][0][3];h[ix][0][1]=D*hnew[ix][0][3];
       f[ix][0][2]=D*fnew[ix][0][4];h[ix][0][2]=D*hnew[ix][0][4];
       f[ix][0][3]=D*fnew[ix][0][1];h[ix][0][3]=D*hnew[ix][0][1];
       f[ix][0][4]=D*fnew[ix][0][2];h[ix][0][4]=D*hnew[ix][0][2];
       f[ix][0][5]=D*fnew[ix][0][7];h[ix][0][5]=D*hnew[ix][0][7];
       f[ix][0][6]=D*fnew[ix][0][8];h[ix][0][6]=D*hnew[ix][0][8];
       f[ix][0][7]=D*fnew[ix][0][5];h[ix][0][7]=D*hnew[ix][0][5];
       f[ix][0][8]=D*fnew[ix][0][6];h[ix][0][8]=D*hnew[ix][0][6];
       //bounceback iy = Ly-1
       f[ix][Ly-1][1]=D*fnew[ix][Ly-1][3];h[ix][Ly-1][1]=D*hnew[ix][Ly-1][3];
       f[ix][Ly-1][2]=D*fnew[ix][Ly-1][4];h[ix][Ly-1][2]=D*hnew[ix][Ly-1][4];
       f[ix][Ly-1][3]=D*fnew[ix][Ly-1][1];h[ix][Ly-1][3]=D*hnew[ix][Ly-1][1];
       f[ix][Ly-1][4]=D*fnew[ix][Ly-1][2];h[ix][Ly-1][4]=D*hnew[ix][Ly-1][2];
       f[ix][Ly-1][5]=D*fnew[ix][Ly-1][7];h[ix][Ly-1][5]=D*hnew[ix][Ly-1][7];
       f[ix][Ly-1][6]=D*fnew[ix][Ly-1][8];h[ix][Ly-1][6]=D*hnew[ix][Ly-1][8];
       f[ix][Ly-1][7]=D*fnew[ix][Ly-1][5];h[ix][Ly-1][7]=D*hnew[ix][Ly-1][5];
       f[ix][Ly-1][8]=D*fnew[ix][Ly-1][6];h[ix][Ly-1][8]=D*hnew[ix][Ly-1][6];
      }
}
void LB::Init(double rho0,double Ux0,double Uy0){
  double phi0,gr_x,gr_y,norm0,nx0,ny0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      phi0 =0.5-0.5*tanh(2*(double)(iy-Ly/2)/W); 
      rho0 =(rho_l-rho_g)*phi0+rho_g;
      gr_x=dphi_dx(ix,iy,false);  gr_y=dphi_dy(ix,iy,false);       
      norm0=norm_gr(gr_x,gr_y);  nx0=nx(gr_x,norm0);  ny0=ny(gr_y,norm0);
      for(int i=0;i<Q;i++){
	      f[ix][iy][i]=feq(rho0,Ux0,Uy0,i);
	      h[ix][iy][i]=heq(phi0,Ux0,Uy0,nx0,ny0,i);
      }
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
  ofstream MiArchivo(NombreArchivo); double rho0,Ux0,Uy0; double Fx,Fy;
  double phi0,gr_x,gr_y,norm0,nx0,ny0;
  int ix=0;
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);  phi0=phi(ix,iy,true);
      gr_x=dphi_dx(ix,iy,true);  gr_y=dphi_dy(ix,iy,true);       
      norm0=norm_gr(gr_x,gr_y);  nx0=nx(gr_x,norm0);  ny0=ny(gr_y,norm0);
      Fx=gx*rho0;  Fy=gy*rho0; 
      Ux0=Jx(ix,iy,true,Fx)/rho0;  Uy0=Jy(ix,iy,true,Fy)/rho0;
      MiArchivo<<iy<<" "<<phi0<<" "<<rho0<<" "<<nx0<<" "<<ny0<<endl;
    }
  MiArchivo.close();
}


int main(void){
  LB Geier;
  int t,tmax=1;
  double RHOinicial=1.0, g=0;
  
  Geier.Init(RHOinicial,0,0);
  
  for(t=0;t<tmax;t++){
    Geier.Collision(g,0);
    //Geier.ImposeFields();
    Geier.Advection();
  }
  
  Geier.Print("Geier.dat",g,0);

  return 0;
}


