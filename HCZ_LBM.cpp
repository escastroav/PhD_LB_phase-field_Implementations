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
const double Cs2 = 1.0/U_Cs2; // Cs2 = RT
const double TresUmUsobre2tau=U_Cs2*(1-1.0/(2*tau));

class LB{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
  double g[Lx][Ly][Q], gnew[Lx][Ly][Q]; // f[ix][iy][i]
  //Multiphase parameters
  double kappa = 1.0;
  double phi_l = 1.0, phi_h = 0.0;
  double rho_l = 1.0, rho_h = 0.5;
  double a = 12*Cs2, b = 4.0;
  double tau_l = 0.8, tau_h = 1.1;
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
  double rho(double phi);
  double nu(double phi);
  double tau(double phi);
  //EOS and Psi computations
  double p_th(double phi);  //EOS Carnahan-Starling
  double d_p_th_d_phi(double phi);
  double Psi_phi(double phi, double p_th);
  double d_Psi_d_phi(double phi);
  double Grad_Psi_x(double dPs_dphi, int ix, int iy, bool UseNew);
  double Grad_Psi_y(double dPs_dphi, int ix, int iy, bool UseNew);
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
  void Init(double Ux0,double Uy0);
  void ImposeFields(void);
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
    jy=(Ly+iy+V[1][i])%Ly;
    sum+=w[i]*V[0][i]*phi(jx,jy,UseNew);
  }
  return sum*U_Cs2;
}
double LB::Grad_y_phi(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    jy=(Ly+iy+V[1][i])%Ly;
    sum+=w[i]*V[1][i]*phi(jx,jy,UseNew);
  }
  return sum*U_Cs2;
}
double LB::Laplacian_phi(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    jy=(Ly+iy+V[1][i])%Ly;
    sum+=w[i]*(phi(jx,jy,UseNew)-phi(ix,iy,UseNew));
  }
  return sum*2*U_Cs2; 
}
double LB::Fsx(int ix, int iy, bool UseNew){
 int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    jy=(Ly+iy+V[1][i])%Ly;
    sum+=w[i]*V[0][i]*Laplacian_phi(jx,jy,UseNew);
  } 
  return kappa*phi(ix,iy,UseNew)*sum;
}
double LB::Fsy(int ix, int iy, bool UseNew){
 int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    jy=(Ly+iy+V[1][i])%Ly;
    sum+=w[i]*V[1][i]*Laplacian_phi(jx,jy,UseNew);
  }
  return kappa*phi(ix,iy,UseNew)*sum;
}
double LB::rho(double phi){
  return rho_l+(phi-phi_l)*(rho_h-rho_l)/(phi_h-phi_l);  
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
    if(UseNew) suma+=fnew[ix][iy][i]*V[0][i]; else suma+=f[ix][iy][i]*V[0][i];
  return suma+0.5*Cs2*(Fsx+Fx);
}
double LB::Jy(double Fsy, double Fy, int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[1][i]; else suma+=f[ix][iy][i]*V[1][i];
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
double LB::Psi_phi(double phi, double p_th){
  return p_th - Cs2*phi;
}
double LB::d_Psi_d_phi(double phi){
  return d_p_th_d_phi(phi) - Cs2;
}
double LB::Grad_Psi_x(double phi, int ix, int iy, bool UseNew)
{
  return d_Psi_d_phi(phi) * Grad_x_phi(ix,iy,UseNew);
}
double LB::Grad_Psi_y(double phi, int ix, int iy, bool UseNew)
{
  return d_Psi_d_phi(phi) * Grad_y_phi(ix,iy,UseNew);
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
  return w[i]*p + rho * Gu_G0; 
}
void LB::Collision(double gx,double gy){
  int ix,iy,i; double phi0,gr_x0,gr_y0;
  double Fsx0, Fsy0, Fx, Fy;
  double rho0, rhoRT, nu0, tau0, Ux0,Uy0,p0;
  double UmUtau, Utau;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Compute macroscopic fields
      phi0=phi(ix,iy,false); gr_x0=Grad_Psi_x(phi0,ix,iy,false); gr_y0=Grad_Psi_y(phi0,ix,iy,false);
      rho0=rho(phi0);  nu0=nu(phi0);  tau0=tau(phi0); rhoRT = rho0*Cs2; UmUtau=1.0-1.0/tau0; Utau=1.0/tau0;
      Fsx0=Fsx(ix,iy,false); Fsy0=Fsy(ix,iy,false); Fx=gx*rho0; Fy=gy*rho0;
      Ux0=Jx(Fsx0,Fx,ix,iy,false)/rhoRT;  Uy0=Jy(Fsy0,Fy,ix,iy,false)/rhoRT;
      p0=p(phi0,rho0,gr_x0,gr_y0,Ux0,Uy0,ix,iy,false);
      for(i=0;i<Q;i++){
    	  fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(Ux0,Uy0,phi0,i)+Fi(tau0,Ux0,Uy0,gr_x0,gr_y0,i);
    	  gnew[ix][iy][i]=UmUtau*g[ix][iy][i]+Utau*geq(p0,rho0,Ux0,Uy0,i)+Gi(tau0,Ux0,Uy0,gr_x0,gr_y0,Fsx0,Fsy0,Fx,Fy,i);
      }
    }
}
void LB::Advection(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++){
	      f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
	      g[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=gnew[ix][iy][i];
      }
}
void LB::Init(double Ux0,double Uy0){
  int W = 3; double phi0, rho0, p0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      phi0 = (phi_l+phi_h)*0.5 + (phi_l-phi_h)*0.5*tanh(2*(ix-Lx/2)/W);
      rho0 = (rho_l+rho_h)*0.5 + (rho_l-rho_h)*0.5*tanh(2*(ix-Lx/2)/W);
      for(int i=0;i<Q;i++){
	      f[ix][iy][i]=feq(Ux0,Uy0,phi0,i);
	      g[ix][iy][i]=geq(0.0,rho0,Ux0,Uy0,i);
        }
      }
}
void LB::ImposeFields(void){
  int i,ix,iy; double phi0, rho0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      phi0=phi(ix,iy,false);
      //Walls
      if(iy==0 || iy==Ly-1) 
	for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(0,0,phi0,i);
    }
}
void LB::Print(const char * NombreArchivo,double gx,double gy){
  ofstream MiArchivo(NombreArchivo); double phi0,rho0,rhoRT,Ux0; 
  double Fsx, Fx;
  int ix=0;
    for(int iy=0;iy<Ly;iy++){
      phi0=phi(ix,iy,true);  rho0=rho(phi0); rhoRT = rho0*Cs2;
      //Fsx0=Fsx(ix,iy,true);  Fx=gx*rho0;
      //Ux0=Jx(Fsx,Fx,ix,iy,true)/rhoRT;
      MiArchivo<<iy<<" "<<rho0<<endl;
    }
  MiArchivo.close();
}


int main(void){
  LB Aire;
  int t,tmax=100;
  double g=0.01;
  
  Aire.Init(0,0);
  
  for(t=0;t<tmax;t++){
    Aire.Collision(g,0);
    Aire.ImposeFields();
    Aire.Advection();
  }
  
  Aire.Print("Aire.dat",g,0);

  return 0;
}


