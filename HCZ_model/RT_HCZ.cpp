#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=64;
const int Ly=128;

const int Q=9;

const double U_Cs2 = 3.0;
const double Cs2 = 1.0/U_Cs2; // Cs2 = RT

class LB{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
  double g[Lx][Ly][Q], gnew[Lx][Ly][Q]; // f[ix][iy][i]
  //Multiphase parameters
  double kappa = 1e-2;
  double phi_l = 0.024, phi_h = 0.251;
  double rho_l = 0.04, rho_h = 0.12;
  double a = 12*Cs2, b = 4.0;
  double tau_l = 0.56, tau_h = 0.56;
  double nu_l = Cs2*(tau_l - 0.5);
  double nu_h = Cs2*(tau_h - 0.5);
public:
  LB(void);
  double phi(int ix,int iy,bool UseNew);
  //Differencial operators
  double Grad_x_phi(int ix, int iy, bool UseNew);
  double Grad_y_phi(int ix, int iy, bool UseNew);
  double Laplacian_phi(int ix, int iy, bool UseNew);
  //Vector fields from distributions f and g
  double Fsx(int ix, int iy, bool UseNew);
  double Fsy(int ix, int iy, bool UseNew);
  double Jx(double Fsx, double Fx, int ix,int iy,bool UseNew); //Jx = rho*Cs^2*U_x
  double Jy(double Fsy, double Fy, int ix,int iy,bool UseNew); //Jy = rho*Cs^2*U_y
  //Scalar physical fields interpolated from phi
  double rho(int ix, int iy, bool UseNew);
  double nu(double phi);
  double tau(double phi);
  //EOS and Psi computations
  double p_th(double phi);  //EOS Carnahan-Starling
  double d_p_th_d_phi(double phi);
  double Psi_phi(int ix, int iy, bool UseNew);
  double Psi_rho(int ix, int iy, bool UseNew);
  double d_Psi_d_phi(double phi);
  double Grad_Psi_phi_x(int ix, int iy, bool UseNew);
  double Grad_Psi_phi_y(int ix, int iy, bool UseNew);
  double Grad_Psi_rho_x(int ix, int iy, bool UseNew);
  double Grad_Psi_rho_y(int ix, int iy, bool UseNew);
  double p(double phi,double rho,double gr_x,double gr_y,double Ux0,double Uy0,int ix,int iy,bool UseNew);
  //Forcing terms and equilibria
  double Gamma(double Ux0,double Uy0,int i);
  double Fi(double tau,double Ux0,double Uy0,double gr_x,double gr_y,int i);
  double Gi(double tau,double Ux0,double Uy0,double gr_x, double gr_y,double Fsx,double Fsy,double Fx,double Fy,int i);
  double feq(double Ux0, double Uy0, double phi, int i);
  double geq(double p, double rho, double Ux0, double Uy0, int i);
  //Automata main steps
  void Collision(double gx,double gy);
  void Advection(void);
  void Init(double Ux0,double Uy0,double gx,double gy);
  void ImposeFields(double gx, double gy);
  void Print(const char * NombreArchivo,double gx,double gy);
};
LB::LB(void){
  //Load the weights
  w[0]=4/9.0; 
  w[1]=w[2]=w[3]=w[4]=1/9.0;
  w[5]=w[6]=w[7]=w[8]=1/36.0;
  //Load D2Q9 vectors
  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;

  V[0][5]=1;  V[0][6]=-1; V[0][7]=-1; V[0][8]=1;
  V[1][5]=1;  V[1][6]=1;  V[1][7]=-1; V[1][8]=-1;
}
double LB::phi(int ix,int iy,bool UseNew){
  int i; double sum;
  for(sum=0,i=0;i<Q;i++)
    if(UseNew) sum+=fnew[ix][iy][i]; else sum+=f[ix][iy][i];
  return sum;
}
double LB::Grad_x_phi(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[0][i]*phi(jx,jy,UseNew);
  }
  return sum*U_Cs2;
}
double LB::Grad_y_phi(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[1][i]*phi(jx,jy,UseNew);
  }
  return sum*U_Cs2;
}
double LB::Laplacian_phi(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0) sum+=w[i]*(phi_l-phi(ix,iy,UseNew));
    else if(iy==Ly-1) sum+=w[i]*(phi_h-phi(ix,iy,UseNew));
    else sum+=w[i]*(phi(jx,jy,UseNew)-phi(ix,iy,UseNew));
  }
  return sum*2*U_Cs2; 
}
double LB::Fsx(int ix, int iy, bool UseNew){
 int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[0][i]*Laplacian_phi(jx,jy,UseNew);
  } 
  return kappa*phi(ix,iy,UseNew)*sum*U_Cs2;
}
double LB::Fsy(int ix, int iy, bool UseNew){
 int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[1][i]*Laplacian_phi(jx,jy,UseNew);
  }
  return kappa*phi(ix,iy,UseNew)*sum*U_Cs2;
}
double LB::rho(int ix, int iy, bool UseNew){
  return rho_l+(phi(ix, iy, UseNew)-phi_l)*(rho_h-rho_l)/(phi_h-phi_l);  
}
double LB::nu(double phi){
  return nu_l+(phi-phi_l)*(nu_h-nu_l)/(phi_h-phi_l);  
}
double LB::tau(double phi){
  return tau_l+(phi-phi_l)*(tau_h-tau_l)/(phi_h-phi_l);  
}
double LB::Jx(double Fsx, double Fx, int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=gnew[ix][iy][i]*V[0][i]; else suma+=g[ix][iy][i]*V[0][i];
  return suma+0.5*Cs2*(Fsx+Fx);
}
double LB::Jy(double Fsy, double Fy, int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=gnew[ix][iy][i]*V[1][i]; else suma+=g[ix][iy][i]*V[1][i];
  return suma+0.5*Cs2*(Fsy+Fy);
}
double LB::p_th(double phi){
  double phi_p = b*phi*0.25;
  double num=1+phi_p+phi_p*phi_p-phi_p*phi_p*phi_p;
  double den=(1-phi_p)*(1-phi_p)*(1-phi_p); 
  return phi*Cs2*num/den - a*phi*phi;
}
double LB::d_p_th_d_phi(double phi){
  double phi_p = b*phi*0.25;
  double num=1+phi_p+phi_p*phi_p-phi_p*phi_p*phi_p;
  double d_num=1+2*phi_p-3*phi_p*phi_p;
  double den=(1-phi_p)*(1-phi_p)*(1-phi_p);
  double d_den=(1-phi_p)*(1-phi_p)*(1-phi_p)*(1-phi_p);
  return Cs2*d_num/den + num/d_den - 2*a*phi;
}
double LB::Psi_phi(int ix, int iy, bool UseNew){
  double phi0;
  if(iy==0) phi0=phi_l; else if(iy==Ly-1) phi0=phi_h; else phi0=phi(ix,iy,UseNew);
  return p_th(phi0) - Cs2*phi(ix,iy,UseNew);
}
double LB::Psi_rho(int ix, int iy, bool UseNew){
  double rho0;
  if(iy==0) rho0=rho_l; else if(iy==Ly-1) rho0=rho_h; else rho0=rho(ix,iy,UseNew);
  return p_th(rho0) - Cs2*phi(ix,iy,UseNew);
}
double LB::d_Psi_d_phi(double phi){
  return d_p_th_d_phi(phi) - Cs2;
}
double LB::Grad_Psi_phi_x(int ix, int iy, bool UseNew)
{
 int i; double sum;
 int jx, jy;
 for(sum=0,i=0;i<Q;i++){
   jx=(Lx+ix+V[0][i])%Lx;
   //if(ix==0||ix==Lx-1) jx=ix;
   jy=(Ly+iy+V[1][i])%Ly;
   if(iy==0||iy==Ly-1) jy=iy;
   sum+=w[i]*V[0][i]*Psi_phi(jx,jy,UseNew);
 }
 return sum*U_Cs2;
}
double LB::Grad_Psi_phi_y(int ix, int iy, bool UseNew)
{
 int i; double sum;
 int jx, jy;
 for(sum=0,i=0;i<Q;i++){
   jx=(Lx+ix+V[0][i])%Lx;
   //if(ix==0||ix==Lx-1) jx=ix;
   jy=(Ly+iy+V[1][i])%Ly;
   if(iy==0||iy==Ly-1) jy=iy;
   sum+=w[i]*V[1][i]*Psi_phi(jx,jy,UseNew);
 }
 return sum*U_Cs2;
}
double LB::Grad_Psi_rho_x(int ix, int iy, bool UseNew)
{
 int i; double sum;
 int jx, jy;
 for(sum=0,i=0;i<Q;i++){
   jx=(Lx+ix+V[0][i])%Lx;
   //if(ix==0||ix==Lx-1) jx=ix;
   jy=(Ly+iy+V[1][i])%Ly;
   if(iy==0||iy==Ly-1) jy=iy;
   sum+=w[i]*V[0][i]*Psi_rho(jx,jy,UseNew);
 }
 return sum*U_Cs2;
}
double LB::Grad_Psi_rho_y(int ix, int iy, bool UseNew)
{
 int i; double sum;
 int jx, jy;
 for(sum=0,i=0;i<Q;i++){
   jx=(Lx+ix+V[0][i])%Lx;
   //if(ix==0||ix==Lx-1) jx=ix;
   jy=(Ly+iy+V[1][i])%Ly;
   if(iy==0||iy==Ly-1) jy=iy;
   sum+=w[i]*V[1][i]*Psi_rho(jx,jy,UseNew);
 }
 return sum*U_Cs2;
}
double LB::p(double phi, double rho, double gr_x, double gr_y, double Ux0, double Uy0, int ix,int iy,bool UseNew){
  int i; double sum;
  double UdotGr = Ux0*gr_x + Uy0*gr_y;
  for(sum=0,i=0;i<Q;i++)
    if(UseNew) sum+=gnew[ix][iy][i]; else sum+=g[ix][iy][i];
  return sum-0.5*UdotGr;
}
double LB::Gamma(double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return w[i]*(1+U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2);
}
double LB::Fi(double tau,double Ux0,double Uy0, double gr_x, double gr_y, int i){
  double VUdotGr=(V[0][i]-Ux0)*gr_x+(V[1][i]-Uy0)*gr_y;
  double Gu=Gamma(Ux0,Uy0,i);
  double UmU2tau = 1.0 - 1.0/(2.0*tau);
  return U_Cs2*UmU2tau*VUdotGr*Gu;
}
double LB::Gi(double tau,double Ux0,double Uy0,double gr_x, double gr_y,double Fsx,double Fsy,double Fx,double Fy,int i){
  double VU_x=(V[0][i]-Ux0),VU_y=(V[1][i]-Uy0);
  double Gu=Gamma(Ux0,Uy0,i);
  double Gu_G0=Gamma(Ux0,Uy0,i)-Gamma(0,0,i);
  double VUdotFG=VU_x*(Fsx+Fx)+VU_y*(Fsy+Fy);
  double VUdotGr=VU_x*gr_x + VU_y*gr_y; 
  double UmU2tau = 1.0 - 1.0/(2.0*tau);
  return UmU2tau*(VUdotFG*Gu-VUdotGr*Gu_G0);
}
double LB::feq(double Ux0, double Uy0, double phi, int i){
  return phi*Gamma(Ux0,Uy0,i);
}
double LB::geq(double p, double rho, double Ux0, double Uy0, int i){
  double Gu=Gamma(Ux0,Uy0,i);
  double Gu_G0=Gamma(Ux0,Uy0,i)-Gamma(0,0,i);
  return w[i]*p + rho * Gu_G0 * Cs2; // include RT=Cs2 besides rho!!!! 
}
void LB::Collision(double gx,double gy){
  int ix,iy,i; double phi0,gr_phi_x0,gr_phi_y0,gr_rho_x0,gr_rho_y0;
  double Fsx0, Fsy0, Fx, Fy;
  double rho0, rhoRT, nu0, tau0, Ux0,Uy0,p0;
  double UmUtau, Utau;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Compute macroscopic fields
      phi0=phi(ix,iy,false); rho0=rho(ix,iy,false); 
      gr_phi_x0=Grad_Psi_phi_x(ix,iy,false); gr_phi_y0=Grad_Psi_phi_y(ix,iy,false);
      gr_rho_x0=Grad_Psi_rho_x(ix,iy,false); gr_phi_y0=Grad_Psi_rho_y(ix,iy,false);
      tau0=tau(phi0);  nu0=nu(phi0); rhoRT = rho0*Cs2; UmUtau=1.0-1.0/tau0; Utau=1.0/tau0;
      Fsx0=Fsx(ix,iy,false); Fsy0=Fsy(ix,iy,false);
      Fx=gx*rho0; Fy=gy*rho0;
      Ux0=Jx(Fsx0,Fx,ix,iy,false)/rhoRT;  Uy0=Jy(Fsy0,Fy,ix,iy,false)/rhoRT;
      p0=p(phi0,rho0,gr_rho_x0,gr_rho_y0,Ux0,Uy0,ix,iy,false);
      for(i=0;i<Q;i++){
    	  fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(Ux0,Uy0,phi0,i)-Fi(tau0,Ux0,Uy0,gr_phi_x0,gr_phi_y0,i);
    	  gnew[ix][iy][i]=UmUtau*g[ix][iy][i]+Utau*geq(p0,rho0,Ux0,Uy0,i)+Gi(tau0,Ux0,Uy0,gr_rho_x0,gr_rho_y0,Fsx0,Fsy0,Fx,Fy,i);
      }
    }
}
void LB::Advection(void){
  double D = 1.0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++){
	      f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
	      g[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=gnew[ix][iy][i];
        //bounceback iy = 0
      f[ix][0][1]=D*fnew[ix][0][3];g[ix][0][1]=D*gnew[ix][0][3];
      f[ix][0][2]=D*fnew[ix][0][4];g[ix][0][2]=D*gnew[ix][0][4];
      f[ix][0][3]=D*fnew[ix][0][1];g[ix][0][3]=D*gnew[ix][0][1];
      f[ix][0][4]=D*fnew[ix][0][2];g[ix][0][4]=D*gnew[ix][0][2];
      f[ix][0][5]=D*fnew[ix][0][7];g[ix][0][5]=D*gnew[ix][0][7];
      f[ix][0][6]=D*fnew[ix][0][8];g[ix][0][6]=D*gnew[ix][0][8];
      f[ix][0][7]=D*fnew[ix][0][5];g[ix][0][7]=D*gnew[ix][0][5];
      f[ix][0][8]=D*fnew[ix][0][6];g[ix][0][8]=D*gnew[ix][0][6];
      //bounceback iy = Ly-1
      f[ix][Ly-1][1]=D*fnew[ix][Ly-1][3];g[ix][Ly-1][1]=D*gnew[ix][Ly-1][3];
      f[ix][Ly-1][2]=D*fnew[ix][Ly-1][4];g[ix][Ly-1][2]=D*gnew[ix][Ly-1][4];
      f[ix][Ly-1][3]=D*fnew[ix][Ly-1][1];g[ix][Ly-1][3]=D*gnew[ix][Ly-1][1];
      f[ix][Ly-1][4]=D*fnew[ix][Ly-1][2];g[ix][Ly-1][4]=D*gnew[ix][Ly-1][2];
      f[ix][Ly-1][5]=D*fnew[ix][Ly-1][7];g[ix][Ly-1][5]=D*gnew[ix][Ly-1][7];
      f[ix][Ly-1][6]=D*fnew[ix][Ly-1][8];g[ix][Ly-1][6]=D*gnew[ix][Ly-1][8];
      f[ix][Ly-1][7]=D*fnew[ix][Ly-1][5];g[ix][Ly-1][7]=D*gnew[ix][Ly-1][5];
      f[ix][Ly-1][8]=D*fnew[ix][Ly-1][6];g[ix][Ly-1][8]=D*gnew[ix][Ly-1][6];
      }
}
void LB::Init(double Ux0,double Uy0,double gx,double gy){
  double W = 5.0; double phi0,rho0,p0,gr_x0,gr_y0,x,y;
  //double Fsx0, Fsy0, Fx, Fy;
  double rhoRT;//,Ux0,Uy0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      x = (double)ix; y = (double)iy;
      phi0 = (phi_h+phi_l)*0.5+(phi_h-phi_l)*0.5*tanh(2*(y-Ly/2-Lx*0.1*cos(2*M_PI*x/Lx))/W);
      rho0 = (rho_h+rho_l)*0.5+(rho_h-rho_l)*0.5*tanh(2*(y-Ly/2-Lx*0.1*cos(2*M_PI*x/Lx))/W);
      gr_x0=Grad_Psi_rho_x(ix,iy,false); gr_y0=Grad_Psi_rho_y(ix,iy,false);
      //rhoRT = rho0*Cs2; Fsx0=Fsx(ix,iy,false); Fsy0=Fsy(ix,iy,false); Fx=gx*rho0; Fy=gy*rho0;
      //Ux0=Jx(Fsx0,Fx,ix,iy,false)/rhoRT;  Uy0=Jy(Fsy0,Fy,ix,iy,false)/rhoRT;
      p0=p(phi0,rho0,gr_x0,gr_y0,Ux0,Uy0,ix,iy,false);
      for(int i=0;i<Q;i++){
	      f[ix][iy][i]=feq(Ux0,Uy0,phi0,i);
	      g[ix][iy][i]=geq(p0,rho0,Ux0,Uy0,i);
        }
      }
}
void LB::ImposeFields(double gx, double gy){
  int i,ix,iy; double phi0,rho0,rhoRT,gr_x0,gr_y0,Ux0,Uy0;
  double Fsx0, Fsy0, Fx, Fy, p0;
  for(ix=0;ix<Lx;ix++){
      rhoRT=rho0*Cs2;
      gr_x0=0; gr_y0=0;
      //Fsx0=Fsx(ix,iy,true);  Fx=gx*rho0;  Fsy0=Fsy(ix,iy,true);  Fy=gy*rho0;
      Ux0=0;  Uy0=0;
	    for(i=0;i<Q;i++){
        fnew[ix][0][i]=feq(0,0,phi_l,i);
	      gnew[ix][0][i]=geq(p(phi_l,rho_l,0,0,0,0,ix,0,false),rho_l,0,0,i);
        fnew[ix][Ly-1][i]=feq(0,0,phi_h,i);
	      gnew[ix][Ly-1][i]=geq(p(phi_h,rho_h,0,0,0,0,ix,Ly-1,false),rho_h,0,0,i);
      }
    }
}
void LB::Print(const char * NombreArchivo,double gx,double gy){
  ofstream MiArchivo(NombreArchivo); 
  double phi0,rho0,rhoRT,Ux0,Uy0,gr_x0,gr_y0; 
  double Fsx0, Fsy0, Fx, Fy,  p0, tau0;
  int ix=0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      phi0=phi(ix,iy,true);  rho0=rho(ix,iy,true); rhoRT = rho0*Cs2;  tau0=tau(phi0);
      Fsx0=Fsx(ix,iy,true); Fsy0=Fsy(ix,iy,true);
      Fx=gx*rho0; Fy=gy*rho0;
      gr_x0=Grad_Psi_phi_x(ix,iy,true); gr_y0=Grad_Psi_phi_y(ix,iy,true);
      Ux0=Jx(Fsx0,Fx,ix,iy,true);  Uy0=Jy(Fsy0,Fy,ix,iy,true);
      p0=p(phi0,rho0,gr_x0,gr_y0,Ux0,Uy0,ix,iy,true);
      MiArchivo<<ix<<"\t"<<iy<<"\t"<<rho0<<"\t"<<phi0<<"\t"<<Uy0<<endl;
    }
  MiArchivo.close();
}


int main(void){
  LB HCZ_RT;
  int t,tmax=4000;
  double g=-6.4e-5;
  
  HCZ_RT.Init(0,0,0,g);
  
  for(t=0;t<tmax;t++){
    HCZ_RT.Collision(0,g);
    HCZ_RT.ImposeFields(g,0);
    HCZ_RT.Advection();
  }
  
  HCZ_RT.Print("rt_4000.dat",0,g);

  return 0;
}


