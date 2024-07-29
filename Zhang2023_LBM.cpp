#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1;
const int Ly=64;

const int Q=9;
const double U_Cs2 = 3.0;
const double Cs2 = 1.0/U_Cs2;
//No parcial miscibility for now!
class LB{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i] fluid
  double h[Lx][Ly][Q], hnew[Lx][Ly][Q]; // h[ix][iy][i] molar fraction
  double g1[Lx][Ly][Q], g1new[Lx][Ly][Q]; // g[ix][iy][i] component concentration
  double g2[Lx][Ly][Q], g2new[Lx][Ly][Q]; // g[ix][iy][i] component concentration
  //Average densities
  double rho1_bar;
  //Adjustable parameters
  double W;               //interface thickness
  double rho_l,rho_g;     //rho_l liquid density; rho_g gas density;
  double sigma;           //surface tension
  double nu_l,nu_g;       //nu_l liquid viscosity; nu_g gas viscosity;
  double M;               //mobility
  double tau_phi=M*U_Cs2+0.5;
  //Auxiliary parameters
  double mu_l,  mu_g;
  double tau_l;
  double tau_g;
  double beta;
  //kappas
  double k_phi;
  double k_c1,k_c2;
  //cache variable phiu for dt
  double old_phiux[Lx][Ly];
  double old_phiuy[Lx][Ly];
  double old_H[Lx][Ly];
public:
  LB(double W,double rho_l0,double rho_g0,double nu_l0,double nu_g0,double M0,double sigma0);
  //Scalar fields
  double phi(int ix,int iy,bool UseNew);
  double c1(int ix,int iy,bool UseNew);
  double c2(int ix,int iy,bool UseNew);
  double lambda(double phi);
  double rho(double c10, double c20);
  double tau(double phi);
  double mu_phi(int ix, int iy, bool UseNew);
  //Gradient components
  double dphi_dx(int ix, int iy, bool UseNew);
  double dc1_dx(int ix, int iy, bool UseNew);
  double dphi_dy(int ix, int iy, bool UseNew);
  double dc2_dx(int ix, int iy, bool UseNew);
  double Laplacian_phi(int ix, int iy, bool UseNew);
  double Laplacian_c1(int ix, int iy, bool UseNew);
  double Laplacian_c2(int ix, int iy, bool UseNew);
  double norm_gr(double gr_x, double gr_y);
  //Free energy 
  double f0(double phi0, double c10, double c20); 
  double df0_dphi(double phi0);
  double df0_dc1(double phi0, double c10, double c20); 
  double df0_dc2(double phi0, double c10, double c20); 
  //Chemical potentials
  double mu_phi(int ix, int iy, bool UseNew);
  double mu_c1(int ix, int iy, bool UseNew);
  double mu_c2(int ix, int iy, bool UseNew);
  //Normal components
  double nx(double gr_x, double norm);
  double ny(double gr_y, double norm);
  //Forces
  double Fsx_phi(double mu, double gr_x);
  double Fsy_phi(double mu, double gr_y);
  double Jx(int ix,int iy,bool UseNew,double Fx);
  double Jy(int ix,int iy,bool UseNew,double Fy);
  double Fi(double tau,double Ux0,double Uy0,double gr_x,double gr_y,double Fsx,double Fsy,double Fx,double Fy,int i);
  double Gi(double nx0,double ny0,double dt_phiux,double dt_phiuy,double lmd,int i,double tau0);
  double Hi(double nx0,double ny0,double dt_phiux,double dt_phiuy,double H,double dH_dt,double lmd,int i,double tau0);
  double p(double Ux0,double Uy0,double gr_x,double gr_y,double rho0,int ix,int iy, bool UseNew);
  double Gamma(double Ux0, double Uy0, int i);
  double si(double Ux0, double Uy0, int i);
  double feq(double p0,double rho0,double Ux0,double Uy0,int i);
  double geq(double rho0,double p0,double Ux0,double Uy0,int i);
  double heq(double phi0,double Ux0,double Uy0,double gr_x,double gr_y,int i);
  void Collision(double gx,double gy);
  void Advection(void);
  void Init(double Ux0,double Uy0);
  void ImposeFields(void);
  void Print(const char * NombreArchivo,double gx,double gy);
};
LB::LB(double W0,double rho_l0,double rho_g0,double nu_l0,double nu_g0,double M0,double sigma0){
  //Setting parameteres up
  W=W0;                         //interface thickness
  rho_l=rho_l0, rho_g=rho_g0;   //rho_l liquid density; rho_g gas density;
  sigma=sigma0;                 //surface tension
  nu_l =nu_l0, nu_g =nu_g0;     //nu_l liquid viscosity; nu_g gas viscosity;
  M=M0;                         //Mobility
  //Average densities
  rho1_bar=1.0; rho2_bar=1.0;
  //Auxiliary parameters
  mu_l =nu_l*rho_l, mu_g =nu_g*rho_g;
  tau_l = U_Cs2*nu_l + 0.5;
  tau_g = U_Cs2*nu_g + 0.5;
  beta=12.0*sigma/W;
  //Kappas
  k_phi = 1.5*sigma*W;
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
double LB::phi(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=hnew[ix][iy][i]; else suma+=h[ix][iy][i];
  return suma;
}
double LB::c1(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=g1new[ix][iy][i]; else suma+=g1[ix][iy][i];
  return suma;
}
double LB::c2(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=g2new[ix][iy][i]; else suma+=g2[ix][iy][i];
  return suma;
}
double LB::rho(double c10,double c20){
  return rho1_bar*c10+rho2_bar*c20;
}
double LB::lambda(double phi){
  return 4*phi*(1.0 - phi)/W; 
}
double LB::tau(double phi){
  return phi*(tau_l - tau_g) + tau_g;
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
double LB::dc1_dx(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[0][i]*c1(jx,jy,UseNew);
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
double LB::dc2_dx(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[0][i]*c2(jx,jy,UseNew);
  }
  return sum*U_Cs2;
}
double LB::Laplacian_phi(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*(phi(jx,jy,UseNew)-phi(ix,iy,UseNew));
  }
  return sum*2*U_Cs2;
}
double LB::Laplacian_c1(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*(c1(jx,jy,UseNew)-c1(ix,iy,UseNew));
  }
  return sum*2*U_Cs2;
}
double LB::Laplacian_c2(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*(c2(jx,jy,UseNew)-c2(ix,iy,UseNew));
  }
  return sum*2*U_Cs2;
}
double LB::norm_gr(double gr_x, double gr_y){
  return (gr_x == 0 && gr_y == 0) ? 1 : sqrt(gr_x*gr_x + gr_y*gr_y);
}
double LB::f0(double phi0, double c10, double c20)
{
  double W_phi = phi0*phi0*(1.0 - phi0*phi0);
  double omega_phi = 1e-3;
  return omega_phi*W_phi;
} 
double LB::df0_dphi(double phi0)
{
  double dW_dphi = 2.0*phi0*(1.0 - 2.0*phi0*phi0);
  double omega_phi = 1e-3;
  return omega_phi*dW_phi;
} 
double LB::df0_dc1(double phi0, double c10, double c20){return 0.0;}  //Future implementation 
double LB::df0_dc2(double phi0, double c10, double c20){return 0.0;}  //Future implementation
double LB::mu_phi(int ix, int iy, bool UseNew){
  double phi0=phi(ix,iy,UseNew);
  return df0_dphi(phi0)-k_phi*Laplacian_phi(ix,iy,UseNew); 
}
double LB::mu_c1(int ix, int iy, bool UseNew){
  return -k_c1*Laplacian_c1(ix,iy,UseNew); 
}
double LB::mu_c2(int ix, int iy, bool UseNew){
  return -k_c2*Laplacian_c2(ix,iy,UseNew); 
}
double LB::nx(double gr_x, double norm){
  return gr_x / norm;
}
double LB::ny(double gr_y, double norm){
  return gr_y / norm;
}
double LB::Fsx(int ix, int iy, bool UseNew){
  int i; double sum_phi=0,sum_c1=0,sum_c2=0;
  int jx, jy;
  for(i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum_phi+=w[i]*V[0][i]*mu_phi(jx,jy,UseNew);
    sum_c1+=w[i]*V[0][i]*mu_c1(jx,jy,UseNew);
    sum_c2+=w[i]*V[0][i]*mu_c2(jx,jy,UseNew);
  }
  return -U_Cs2*(phi(ix,iy,UseNew)*sum_phi+c1(ix,iy,UseNew)*sum_c1+c2(ix,iy,UseNew)*sum_c2);
}
double LB::Fsy(double mu, double gr_y){
  int i; double sum_phi=0,sum_c1=0,sum_c2=0;
  int jx, jy;
  for(i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum_phi+=w[i]*V[1][i]*mu_phi(jx,jy,UseNew);
    sum_c1+=w[i]*V[0][i]*mu_c1(jx,jy,UseNew);
    sum_c2+=w[i]*V[0][i]*mu_c2(jx,jy,UseNew);
  }
  return -U_Cs2*(phi(ix,iy,UseNew)*sum_phi+c1(ix,iy,UseNew)*sum_c1+c2(ix,iy,UseNew)*sum_c2);
}
double LB::Jx(int ix,int iy,bool UseNew,double Fx){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=gnew[ix][iy][i]*V[0][i]; else suma+=g[ix][iy][i]*V[0][i];
  return suma+0.5*Cs2*Fx;
}
double LB::Jy(int ix,int iy,bool UseNew,double Fy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=gnew[ix][iy][i]*V[1][i]; else suma+=g[ix][iy][i]*V[1][i];
  return suma+0.5*Cs2*Fy;
}
double LB::Gamma(double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return w[i]*(1+U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2);
}
double LB::Fi(double tau,double Ux0,double Uy0,double gr_x, double gr_y,double Fsx,double Fsy,double Fx,double Fy,int i){
  double VU_x=(V[0][i]-Ux0),VU_y=(V[1][i]-Uy0);
  double Gu=Gamma(Ux0,Uy0,i);
  double Gu_G0=Gamma(Ux0,Uy0,i)-Gamma(0,0,i);
  double VUdotFG=VU_x*(Fsx+Fx)+VU_y*(Fsy+Fy);
  double VUdotGr=VU_x*gr_x + VU_y*gr_y;
  double UmU2tau = 1.0 - 1.0/(2.0*tau);
  return UmU2tau*(VUdotFG*Gu-VUdotGr*Gu_G0);
}
double LB::Gi(double nx0,double ny0,double dt_phiux,double dt_phiuy,double lmd,int i,double tau0){
  double DtdotVi=dt_phiux*V[0][i]+dt_phiuy*V[1][i];
  double UmU2tau = 1.0 - 1.0/(2.0*tau0);
  return U_Cs2*UmU2tau*w[i]*DtdotVi;
}
double LB::Hi(double nx0,double ny0,double dt_phiux,double dt_phiuy,double H,double dH_dt,double lmd,int i,double tau0){
  double DtdotVi=dt_phiux*V[0][i]+dt_phiuy*V[1][i];
  double UmU2tau = 1.0 - 1.0/(2.0*tau0);
  return UmU2tau*w[i]*(U_Cs2*DtdotVi+H+0.5*dH_dt);
}
double LB::p(double Ux0,double Uy0,double gr_x,double gr_y,double rho0,int ix,int iy, bool UseNew){
  int i; double sum;
  double UdotGr = Ux0*gr_x + Uy0*gr_y;
  for(sum=0,i=0;i<Q;i++)
    if(UseNew) sum+=fnew[ix][iy][i]; else sum+=f[ix][iy][i];
  return sumr+0.5*Cs2*UdotGr;
}
double LB::feq(double phi0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return phi0*w[i]*(1+U_Cs2*UdotVi);
}
double LB::geq(double rho0,double p0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  double rho_siu=rho0*w[i]*(U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2);
  if(i>0) return p0*w[i]*U_Cs2+rho_siu;
  else return p0*(w[i]-1.0)*U_Cs2+rho_siu;

}
void LB::Collision(double gx,double gy){
  int ix,iy,i; 
  double phi0,rho0,p0,lmd0,mu0; 
  double tau0,Utau,UmUtau,Utau_phi,UmUtau_phi;
  double gr_x,gr_y,norm0,nx0,ny0;
  double Fx,Fy;
  double Ux0,Uy0;
  double phiux0,phiuy0,dt_phiux0,dt_phiuy0;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macroscopicas
      phi0=phi(ix,iy,false);  gr_x=dphi_dx(ix,iy,false);  gr_y=dphi_dy(ix,iy,false);
      norm0=norm_gr(gr_x,gr_y); nx0=nx(gr_x,norm0); ny0=ny(gr_y,norm0);
      mu0=mu_phi(ix,iy,false);  rho0=rho(phi0); lmd0=lambda(phi0);  
      tau0=tau(phi0); Utau=1.0/tau0;  UmUtau=1.0-Utau;
      Utau_phi=1.0/tau_phi;  UmUtau_phi=1.0-Utau_phi;
      Fx=Fsx(mu0,gr_x)+gx*rho0;  Fy=Fsy(mu0,gr_y)+gy*rho0;
      Ux0=Jx(ix,iy,false,Fx)/rho0;  Uy0=Jy(ix,iy,false,Fy)/rho0;
      p0=p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,false);
      phiux0=Ux0*phi0; dt_phiux0=phiux0 - old_phiux[ix][iy];
      phiuy0=Uy0*phi0; dt_phiuy0=phiuy0 - old_phiuy[ix][iy];
      for(i=0;i<Q;i++){
	      fnew[ix][iy][i]=UmUtau_phi*f[ix][iy][i]+Utau_phi*feq(phi0,Ux0,Uy0,i)+Fi(nx0,ny0,dt_phiux0,dt_phiuy0,lmd0,i,tau_phi);
	      gnew[ix][iy][i]=UmUtau*g[ix][iy][i]+Utau*geq(rho0,p0,Ux0,Uy0,i)+Gi(Ux0,Uy0,Fx,Fy,gr_x,gr_y,i,tau0);
      }
      old_phiux[ix][iy]=phiux0; old_phiuy[ix][iy]=phiuy0;
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
void LB::Init(double Ux0,double Uy0){
  double phi0, rho0, gr_x, gr_y, p0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      phi0 =0.5-0.5*tanh(2*(double)(iy-Ly/2)/W);
      rho0=rho(phi0);
      gr_x=dphi_dx(ix,iy,false);  gr_y=dphi_dy(ix,iy,false);
      p0=p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,false);
      old_phiux[ix][iy] = phi0*Ux0; old_phiuy[ix][iy] = phi0*Uy0;
      for(int i=0;i<Q;i++){
	      f[ix][iy][i]=feq(phi0,Ux0,Uy0,i);
	      g[ix][iy][i]=geq(rho0,p0,Ux0,Uy0,i);
      }
  }
}
void LB::ImposeFields(void){
  int i,ix,iy; double phi0;
  for(ix=0;ix<Lx;ix++){
      //Walls
    	for(i=0;i<Q;i++){
        fnew[ix][0][i]=feq(1,0,0,i);
        fnew[ix][Ly-1][i]=feq(0,0,0,i);
        gnew[ix][0][i]=geq(rho_l,p(0,0,0,0,rho_l,ix,0,false),0,0,i);
        gnew[ix][Ly-1][i]=geq(rho_g,p(0,0,0,0,rho_g,ix,Ly-1,false),0,0,i);
      }
    }
}
void LB::Print(const char * NombreArchivo,double gx,double gy){
  ofstream MiArchivo(NombreArchivo); 
  double phi0,rho0,mu0,tau0,lmd0,Ux0,Uy0,p0; double Fx,Fy,gr_x,gr_y,gr_xloc,gr_yloc;
  double dt_phiux,dt_phiuy;
  int ix=Lx/2;
    for(int iy=0;iy<Ly;iy++){
      phi0=phi(ix,iy,true);
      gr_x=dphi_dx(ix,iy,true);  gr_y=dphi_dy(ix,iy,true);
      mu0=mu_phi(ix,iy,true);  rho0=rho(phi0);tau0=tau(phi0);lmd0=lambda(phi0);
      Fx=Fsx(mu0,gr_x)+gx*rho0;  Fy=Fsy(mu0,gr_y)+gy*rho0;
      Ux0=Jx(ix,iy,true,Fx)/rho0;  Uy0=Jy(ix,iy,true,Fy)/rho0;
      dt_phiux=phi0*Ux0 - old_phiux[ix][iy];
      dt_phiuy=phi0*Uy0 - old_phiuy[ix][iy];
      p0=p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,false);
      MiArchivo<<iy<<" "<<phi0<<" "<<rho0<<" "<<Ux0<<endl;
    }
  MiArchivo.close();
}


int main(void){
  double W=5;
  double rho_l=100, rho_g=1;  //rho_l liquid density; rho_g gas density;
  double sigma=1.0e-3;             //surface tension
  double nu_l =0.1, nu_g =0.01;
  double M=0.1;
  LB Liang(W,rho_l,rho_g,nu_l,nu_g,M,sigma);
  int t,tmax=100000;
  double Uc=1e-4,g=4*Uc*(rho_l*nu_l+rho_g*nu_g)/(Ly*Ly);
  cout << g << endl; 
  Liang.Init(g,0);
  
  for(t=0;t<tmax;t++){
    Liang.Collision(g,0);
    //Liang.ImposeFields();
    Liang.Advection();
  }
  
  Liang.Print("Liang.dat",g,0);

  return 0;
}


