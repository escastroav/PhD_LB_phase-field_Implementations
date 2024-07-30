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
  //Average densities per component
  double rho1_bar=1,rho2_bar=1;
  double eta1_bar=6.7e-2,eta2_bar=6.7e-2;
  //Adjustable parameters
  double W=4;                 //interface thickness
  double sigma=2.5e-3;             //surface tension
  double M_c=9.3e-2;
  //kappas
  double k_phi=1.5*sigma*W;
  double k_c1=k_phi,k_c2=k_phi;
  double M_phi=M_c/(W*W);
  double beta_phi=8*k_phi/(W*W);
  double beta_c=0.1*beta_phi;
  //Relaxation times
  double tau_g1=(M_c*U_Cs2/eta1_bar)+0.5;
  double tau_g2=(M_c*U_Cs2/eta2_bar)+0.5;
  double tau_h=(M_phi*k_phi*U_Cs2)+0.5;//B=0;
  double Utaug1=1.0/tau_g1;  UmUtaug1=1.0-Utaug1;
  double Utaug2=1.0/tau_g2;  UmUtaug2=1.0-Utaug2;
  double Utauh=1.0/tau_h;  UmUtauh=1.0-Utauh;
  //double tau_phi=M_phi*U_Cs2+0.5;
  //cache variable phiu for dt
  double old_c1ux[Lx][Ly];
  double old_c1uy[Lx][Ly];
  double old_c2ux[Lx][Ly];
  double old_c2uy[Lx][Ly];
  double old_phiux[Lx][Ly];
  double old_phiuy[Lx][Ly];
  double old_H[Lx][Ly];
public:
  LB(void);
  //Scalar fields
  double phi(int ix,int iy,bool UseNew);
  double c1(int ix,int iy,bool UseNew);
  double c2(int ix,int iy,bool UseNew);
  double rho(int ix, int iy, bool UseNew);
  double eta(double c10,double c20);
  double tau_f(double eta0,double rho0);
  //Derivative components
  double drho_dx(int ix, int iy, bool UseNew);
  double drho_dy(int ix, int iy, bool UseNew);
  double Laplacian_phi(int ix, int iy, bool UseNew);
  double Laplacian_c1(int ix, int iy, bool UseNew);
  double Laplacian_c2(int ix, int iy, bool UseNew);
  //Free energy 
  double f0(double phi0, double c10, double c20); 
  double df0_dphi(double phi0);
  double df0_dc1(double phi0, double c10, double c20); 
  double df0_dc2(double phi0, double c10, double c20); 
  //Chemical potentials
  double mu_phi(int ix, int iy, bool UseNew);
  double mu_c1(int ix, int iy, bool UseNew);
  double mu_c2(int ix, int iy, bool UseNew);
  //Forces
  double Fsx(int ix, int iy, bool UseNew);
  double Fsy(int ix, int iy, bool UseNew);
  double Jx(int ix,int iy,bool UseNew,double Fx);
  double Jy(int ix,int iy,bool UseNew,double Fy);
  //Auxiliary fields
  double H(int ix, int iy, bool UseNew);
  double Gamma(double Ux0, double Uy0, int i);
  double Fi(double tau,double Ux0,double Uy0,double gr_x,double gr_y,double Fx,double Fy,int i);
  double Gi(double tau,double dt_cjux,double dt_cjuy,int i);
  double Hi(double dt_phiux,double dt_phiuy,double H,double dH_dt,int i);
  double p(double Ux0,double Uy0,double gr_x,double gr_y,double rho0,int ix,int iy, bool UseNew);
  double feq(double Ux0,double Uy0,double p0,double rho0,int i);
  double geq(double Ux0,double Uy0,double cj0,double etaj0,double muj0,int i);
  double heq(double phi0,double Ux0,double Uy0,double gr_x,double gr_y,int i);
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
double eta(double c10,double c20){
  return exp(c10*log(eta1_bar)+c20*log(eta2_bar));//Gij = 0 
}
double LB::tau_f(double eta0,double rho0){
  return eta0*U_Cs2/rho0 + 0.5;
}
double LB::drho_dx(int ix, int iy, bool UseNew){
  int i; double sum_c1,sum_c2;
  int jx, jy;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum_c1+=w[i]*V[0][i]*c1(jx,jy,UseNew);
    sum_c2+=w[i]*V[0][i]*c2(jx,jy,UseNew);
  }
  return U_Cs2*(rho1_bar*sum_c1+rho2_bar*sum_c2);
}
double LB::drho_dy(int ix, int iy, bool UseNew){
  int i; double sum_c1,sum_c2;
  int jx, jy;
  for(sum=0,i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum_c1+=w[i]*V[1][i]*c1(jx,jy,UseNew);
    sum_c2+=w[i]*V[1][i]*c2(jx,jy,UseNew);
  }
  return U_Cs2*(rho1_bar*sum_c1+rho2_bar*sum_c2);
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
double LB::Fsx(int ix, int iy, bool UseNew){
  int i; double sum_phi=0,sum_c1=0;//,sum_c2=0;
  int jx, jy;
  for(i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum_phi+=w[i]*V[0][i]*mu_phi(jx,jy,UseNew);
    sum_c1+=w[i]*V[0][i]*mu_c1(jx,jy,UseNew);
    //sum_c2+=w[i]*V[0][i]*mu_c2(jx,jy,UseNew);
  }
  return -U_Cs2*(phi(ix,iy,UseNew)*sum_phi+c1(ix,iy,UseNew)*sum_c1);//+c2(ix,iy,UseNew)*sum_c2);
}
double LB::Fsy(int ix, int iy, bool UseNew){
  int i; double sum_phi=0,sum_c1=0;//,sum_c2=0;
  int jx, jy;
  for(i=1;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    if(iy==0||iy==Ly-1) jy=iy;
    sum_phi+=w[i]*V[1][i]*mu_phi(jx,jy,UseNew);
    sum_c1+=w[i]*V[0][i]*mu_c1(jx,jy,UseNew);
    //sum_c2+=w[i]*V[0][i]*mu_c2(jx,jy,UseNew);
  }
  return -U_Cs2*(phi(ix,iy,UseNew)*sum_phi+c1(ix,iy,UseNew)*sum_c1);//+c2(ix,iy,UseNew)*sum_c2);
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
double LB::H(int ix,int iy,bool UseNew){
  double mu_phi0 = mu_phi(ix,iy,UseNew);
  double lapl_phi= Laplacian_phi(ix,iy,UseNew);
  return -M_phi*(mu_phi0+k_phi*lapl_phi);
}
double LB::Gamma(double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return w[i]*(1+U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2);
}
double LB::Fi(double tau,double Ux0,double Uy0,double gr_x, double gr_y,double Fx,double Fy,int i){
  double VU_x=(V[0][i]-Ux0),VU_y=(V[1][i]-Uy0);
  double Gu=Gamma(Ux0,Uy0,i);
  double Gu_G0=Gamma(Ux0,Uy0,i)-Gamma(0,0,i);
  double VUdotFG=VU_x*Fx+VU_y*Fy;
  double VUdotGr=VU_x*gr_x + VU_y*gr_y;
  double UmU2tau = 1.0 - 1.0/(2.0*tau);
  return UmU2tau*(VUdotFG*Gu-VUdotGr*Gu_G0);
}
double LB::Gi(double tau,double dt_cjux,double dt_cjuy,int i)
  double DtdotVi=dt_cjux*V[0][i]+dt_cjuy*V[1][i];
  double UmU2tau = 1.0 - 1.0/(2.0*tau);
  return U_Cs2*UmU2tau*w[i]*DtdotVi;
}
double LB::Hi(double dt_phiux,double dt_phiuy,double H,double dH_dt,int i){
  double DtdotVi=dt_phiux*V[0][i]+dt_phiuy*V[1][i];
  double UmU2tau = 1.0 - 1.0/(2.0*tau_h);
  return UmU2tau*w[i]*(U_Cs2*DtdotVi+H+0.5*dH_dt);
}
double LB::p(double Ux0,double Uy0,double gr_x,double gr_y,double rho0,int ix,int iy, bool UseNew){
  int i; double sum;
  double UdotGr = Ux0*gr_x + Uy0*gr_y;
  for(sum=0,i=0;i<Q;i++)
    if(UseNew) sum+=fnew[ix][iy][i]; else sum+=f[ix][iy][i];
  return sumr+0.5*Cs2*UdotGr;
}
double LB::feq(double Ux0,double Uy0,double p0,double rho0,int i){
  double Gu=Gamma(Ux0,Uy0,i);
  double Gu_G0=Gamma(Ux0,Uy0,i)-Gamma(0,0,i);
  return w[i]*p0+rho0*Gu_G0*Cs2;
}
double LB::geq(double Ux0,double Uy0,double cj0,double etaj0,double muj0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i];
  double rho_siu=cj0*w[i]*U_Cs2*UdotVi;
  if(i>0) return w[i]*muj0*etaj0+rho_siu;
  else return etaj0*muj0*(w[i]-1.0)+cj0;
}
double heq(double phi0,double Ux0,double Uy0,double gr_x,double gr_y,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return phi0*w[i]*(1+U_Cs2*UdotVi);
}
void LB::Collision(double gx,double gy){
  int ix,iy,i; 
  double phi0,c10,c20,rho0,p0,eta0;
  double tau_phi0,Utau,UmUtau;
  double gr_x,gr_y,mu_c10,mu_c20;
  double Fx,Fy,Ux0,Uy0;
  double phiux0,phiuy0,dt_phiux0,dt_phiuy0;
  double c1ux0,c1uy0,dt_c1ux0,dt_c1uy0;
  double c2ux0,c2uy0,dt_c2ux0,dt_c2uy0;
  double H0,dH_dt;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macroscopicas
      phi0=phi(ix,iy,false);  c10=c1(ix,iy,false);  c20=c2(ix,iy,false);  rho0=rho(ix,iy,false);
      eta0=eta(c10,c20);  tau_phi0=tau_phi(eta0,rho0);  Utau=1.0/tau_phi0;  UmUtau=1.0-Utau;
      gr_x=drho_dx(ix,iy,false);  gr_y=drho_dy(ix,iy,false);
      mu_c10=mu_c1(ix,iy,false);  mu_c20=mu_c2(ix,iy,false);
      Fx=Fsx(ix,iy,false)+gx*rho0;  Fy=Fsy(ix,iy,false)+gy*rho0;
      Ux0=Jx(ix,iy,false,Fx)*U_Cs2/rho0;  Uy0=Jy(ix,iy,false,Fy)*U_Cs2/rho0;
      p0=p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,false);
      //Computing the time derivatives
      phiux0=Ux0*phi0; dt_phiux0=phiux0 - old_phiux[ix][iy];
      phiuy0=Uy0*phi0; dt_phiuy0=phiuy0 - old_phiuy[ix][iy];
      c1ux0=Ux0*c10; dt_c1ux0=c1ux0 - old_c1ux[ix][iy];
      c1uy0=Uy0*c10; dt_c1uy0=c1uy0 - old_c1uy[ix][iy];
      c2ux0=Ux0*c20; dt_c2ux0=c2ux0 - old_c2ux[ix][iy];
      c2uy0=Uy0*c20; dt_c2uy0=c2uy0 - old_c2uy[ix][iy];
      H0=H(ix,iy,false);  dH_dt=H0 - old_H[ix][iy];
      for(i=0;i<Q;i++){
	      fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(Ux0,Uy0,p0,rho0,i)+Fi(tau_phi0,Ux0,Uy0,gr_x,gr_y,Fx,Fy,i);
	      g1new[ix][iy][i]=UmUtaug1*g1[ix][iy][i]+Utaug1*geq(Ux0,Uy0,c10,eta1_bar,mu_c10,i)+Gi(tau_g1,dt_c1ux0,dt_c1uy0,i);
	      g2new[ix][iy][i]=UmUtaug2*g2[ix][iy][i]+Utaug2*geq(Ux0,Uy0,c20,eta2_bar,mu_c20,i)+Gi(tau_g2,dt_c2ux0,dt_c2uy0,i);
	      hnew[ix][iy][i]=UmUtauh*h[ix][iy][i]+Utauh*heq(phi0,Ux0,Uy0,gr_x,gr_y,i)+Hi(dt_phiux0,dt_phiuy0,H0,dH_dt,i);
      }
      old_phiux[ix][iy]=phiux0; old_phiuy[ix][iy]=phiuy0;
      old_c1ux[ix][iy]=c1ux0; old_c1uy[ix][iy]=c1uy0;
      old_c2ux[ix][iy]=c2ux0; old_c2uy[ix][iy]=c2uy0;
      old_H[ix][iy]=H0;
    }
}
void LB::Advection(void){
  double D = 1.0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++){
	      f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
	      g1[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=g1new[ix][iy][i];
	      g2[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=g2new[ix][iy][i];
	      h[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=hnew[ix][iy][i];
      }
}
void LB::Init(double Ux0,double Uy0){
  double cl0=0.1,cg0=0.8188;
  double phi0, c10, c20, rho0;
  double gr_x, gr_y, mu_c10, mu_c20, p0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      phi0 =0.5*tanh(2*(double)(iy-0.4*Ly)/W)-0.5*tanh(2*(double)(iy-0.6*Ly)/W);
      c10 =(cg0-cl0)*phi0+cl0;  c20 =(cl0-cg0)*phi0+cg0;  rho0=c10*rho1_bar+c20*rho2_bar;
      gr_x=drho_dx(ix,iy,false);  gr_y=drho_dy(ix,iy,false);
      mu_c10=mu_c1(ix,iy,false);  mu_c20=mu_c2(ix,iy,false);
      p0=p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,false);
      //Initialize the time derivatives!
      old_phiux[ix][iy] = phi0*Ux0; old_phiuy[ix][iy] = phi0*Uy0;
      old_c1ux[ix][iy]=c10*Ux0; old_c1uy[ix][iy]=c10*Uy0;
      old_c2ux[ix][iy]=c20*Ux0; old_c2uy[ix][iy]=c20*Uy0;
      old_H[ix][iy]=H(ix,iy,false);
      for(int i=0;i<Q;i++){
	      f[ix][iy][i]=feq(Ux0,Uy0,p0,rho0,i);
	      g1[ix][iy][i]=geq(Ux0,Uy0,c10,eta1_bar,mu_c10,i);
	      g2[ix][iy][i]=geq(Ux0,Uy0,c20,eta2_bar,mu_c20,i);
        h[ix][iy][i]=heq(phi0,Ux0,Uy0,gr_x,gr_y,i);
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
  double phi0,c10,c20,rho0;
  int ix=Lx/2;
    for(int iy=0;iy<Ly;iy++){
      phi0=phi(ix,iy,true); c10=c1(ix,iy,true); c20=c2(ix,iy,true); rho0=rho(ix,iy,true);
      MiArchivo<<iy<<" "<<phi0<<" "<<rho0<<" "<<c10<<" "<<c20<<endl;
    }
  MiArchivo.close();
}


int main(void){
  LB Liang();
  int t,tmax=1;
  Liang.Init(0,0);
  
  for(t=0;t<tmax;t++){
    Liang.Collision(0,0);
    //Liang.ImposeFields();
    Liang.Advection();
  }
  
  Liang.Print("Zhang.dat",g,0);

  return 0;
}


