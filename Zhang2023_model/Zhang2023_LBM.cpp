#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1;
const int Ly=320;

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
  //double g2[Lx][Ly][Q], g2new[Lx][Ly][Q]; // g[ix][iy][i] component concentration
  //Average densities per component
  double rho1_bar=1,rho2_bar=1;
  double eta1_bar=1.7e-1,eta2_bar=1.7e-1;
  double cl1=0.25,cg1=0.85;
  double cl2=1.0-cl1,cg2=1.0-cg1;
  //Adjustable parameters
  double W=4;                 //interface thickness
  double sigma=2.5e-3;             //surface tension
  double M_c=9.3e-2;
  double B=0.1;
  //kappas
  double k_phi=1.5*sigma*W;
  double k_c1=k_phi,k_c2=k_phi;
  double M_phi=M_c/(W*W);
  double beta_phi=8*k_phi/(W*W);
  double beta_c=0.1*beta_phi;
  //Relaxation times
  double tau_g1=(M_c*U_Cs2/eta1_bar)+0.5;
  double tau_g2=(M_c*U_Cs2/eta2_bar)+0.5;
  double tau_h=(M_phi*k_phi*U_Cs2)+0.5+B;//B=0;
  double Utaug1=1.0/tau_g1,  UmUtaug1=1.0-Utaug1;
  double Utaug2=1.0/tau_g2,  UmUtaug2=1.0-Utaug2;
  double Utauh=1.0/tau_h,  UmUtauh=1.0-Utauh;
  //omega for f0 (from ref [] w_phi = sigma/W but no clue about w_mix!)
  double omega_phi = beta_phi;
  double omega_mix = beta_c;
  double lmd_21_lmd_11 = -7.0;
  //cache variable phiu for dt
  double old_c1ux[Lx][Ly];
  double old_c1uy[Lx][Ly];
  //double old_c2ux[Lx][Ly];
  //double old_c2uy[Lx][Ly];
  double old_phiux[Lx][Ly];
  double old_phiuy[Lx][Ly];
  double old_H[Lx][Ly];
  double old_p[Lx][Ly];
public:
  LB(void);
  //Scalar fields
  double phi(int ix,int iy,bool UseNew);
  double c1(int ix,int iy,bool UseNew);
  //double c2(int ix,int iy,bool UseNew);
  double rho(double c10,double c20);
  double eta(double c10,double c20);
  double tau_phi(double eta0,double rho0);
  //Derivative components
  double dphi_dx(int ix, int iy, bool UseNew);
  double dphi_dy(int ix, int iy, bool UseNew);
  double drho_dx(int ix, int iy, bool UseNew);
  double drho_dy(int ix, int iy, bool UseNew);
  double Laplacian_phi(int ix, int iy, bool UseNew);
  double Laplacian_c1(int ix, int iy, bool UseNew);
  double Laplacian_c2(int ix, int iy, bool UseNew);
  //Free energy 
  double f0(double phi0,double c1l,double c2l,double c1g,double c2g); 
  double fl(double c1l,double c2l);
  double fg(double c1g,double c2g);
  double dfl_dc1(double c1l,double c2l);
  double dfg_dc1(double c1g,double c2g);
  double dfl_dc2(double c1l,double c2l);
  double dfg_dc2(double c1g,double c2g);
  double df0_dphi(double phi0,double c1l,double c2l,double c1g,double c2g);
  double df0_dc1(double phi0,double c1l,double c2l,double c1g,double c2g);
  double df0_dc2(double phi0,double c1l,double c2l,double c1g,double c2g);
  //Chemical potentials
  double mu_phi(int ix, int iy, bool UseNew);
  double mu_c1(int ix, int iy, bool UseNew);
  double mu_c2(int ix, int iy, bool UseNew);
  double mutilde_cj(double mu_cj,double p0);
  //Forces and vector fields
  double Fsx(int ix, int iy, bool UseNew);
  double Fsy(int ix, int iy, bool UseNew);
  double Jx(int ix,int iy,bool UseNew,double Fx);
  double Jy(int ix,int iy,bool UseNew,double Fy);
  //Auxiliary fields
  double H(int ix, int iy, bool UseNew);
  double Gamma(double Ux0, double Uy0, int i);
  //Forcing terms
  double Fi(double tau,double Ux0,double Uy0,double gr_x,double gr_y,double Fx,double Fy,int i);
  double Gi(double tau,double dt_cjux,double dt_cjuy,int i);
  double Hi(double dt_phiux,double dt_phiuy,double H,double dH_dt,int i);
  //Pressure
  double p(double Ux0,double Uy0,double gr_x,double gr_y,double rho0,int ix,int iy, bool UseNew);
  //Equilibrium functions
  double feq(double Ux0,double Uy0,double p0,double rho0,int i);
  double geq(double Ux0,double Uy0,double cj0,double etaj0,double muj0,int i);
  double heq(double phi0,double Ux0,double Uy0,double gr_x,double gr_y,int i);
  //LB steps
  void Collision(int t,double gx,double gy);
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
//Scalar fields
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
/**double LB::c2(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=g2new[ix][iy][i]; else suma+=g2[ix][iy][i];
  return suma;
}**/
double LB::rho(double c10,double c20){
  return rho1_bar*c10+rho2_bar*c20;
}
double LB::eta(double c10,double c20){
  return exp(c10*log(eta1_bar)+c20*log(eta2_bar));//Gij = 0 
}
double LB::tau_phi(double eta0,double rho0){
  return eta0*U_Cs2/rho0 + 0.5;
}
//Derivative components
double LB::dphi_dx(int ix, int iy, bool UseNew){
  int i; double sum=0;
  int jx, jy;
  for(i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[0][i]*phi(jx,jy,UseNew);
  }
  return U_Cs2*sum;
}
double LB::dphi_dy(int ix, int iy, bool UseNew){
  int i; double sum=0;
  int jx, jy;
  for(i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*V[1][i]*phi(jx,jy,UseNew);
  }
  return U_Cs2*sum;
}
double LB::drho_dx(int ix, int iy, bool UseNew){
  int i; double sum_c1=0,sum_c2=0;
  int jx, jy;
  for(i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum_c1+=w[i]*V[0][i]*c1(jx,jy,UseNew);
    sum_c2+=w[i]*V[0][i]*(1.0-c1(jx,jy,UseNew));
  }
  return U_Cs2*(rho1_bar*sum_c1+rho2_bar*sum_c2);
}
double LB::drho_dy(int ix, int iy, bool UseNew){
  int i; double sum_c1=0,sum_c2=0;
  int jx, jy;
  for(i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum_c1+=w[i]*V[1][i]*c1(jx,jy,UseNew);
    sum_c2+=w[i]*V[1][i]*(1.0-c1(jx,jy,UseNew));
  }
  return U_Cs2*(rho1_bar*sum_c1+rho2_bar*sum_c2);
}
double LB::Laplacian_phi(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*(phi(jx,jy,UseNew)-phi(ix,iy,UseNew));
  }
  return sum*2*U_Cs2;
}
double LB::Laplacian_c1(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*(c1(jx,jy,UseNew)-c1(ix,iy,UseNew));
  }
  return sum*2*U_Cs2;
}
double LB::Laplacian_c2(int ix, int iy, bool UseNew){
  int i; double sum;
  int jx, jy;
  for(sum=0,i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum+=w[i]*(c1(ix,iy,UseNew)-c1(jx,jy,UseNew));
  }
  return sum*2*U_Cs2;
}
//Free Energy
double LB::f0(double phi0,double c1l,double c2l,double c1g,double c2g)
{
  double W_phi = phi0*phi0*(1.0 - phi0*phi0)*(1.0 - phi0*phi0);
  double g=phi0*phi0*(3.0 - 2.0*phi0);
  return omega_phi*W_phi+omega_mix*(fl(c1l,c2l)*(1-g)+fg(c1g,c2g)*g);
} 
double LB::fl(double c1l,double c2l){
  double lmd = lmd_21_lmd_11;
  double v1v2= rho2_bar/rho1_bar;
  double v2v1= rho1_bar/rho2_bar;
  double term1=c1l*log(c1l)+c2l*log(c2l);
  double A_12 =v2v1*exp(-lmd), A_21=v1v2*exp(lmd);
  double term2=c1l*log(c1l+A_12*c2l)+c1l*log(c1l*A_21+c2l);
  return term1-term2;
}
double LB::fg(double c1g,double c2g){
  double lmd = lmd_21_lmd_11;
  double v1v2= rho2_bar/rho1_bar;
  double v2v1= rho1_bar/rho2_bar;
  double term1=c1g*log(c1g)+c2g*log(c2g);
  double A_12 =v2v1*exp(lmd), A_21=v1v2*exp(-lmd);
  double term2=c1g*log(c1g+A_12*c2g)+c2g*log(c1g*A_21+c2g);
  return term1-term2;
}
double LB::dfl_dc1(double c1l,double c2l){
  double lmd = lmd_21_lmd_11;
  double v1v2= rho2_bar/rho1_bar;
  double v2v1= rho1_bar/rho2_bar;
  double term1=log(c1l)+1;
  double A_12 =v2v1*exp(-lmd), A_21=v1v2*exp(lmd);
  double term2=log(c1l+A_12*c2l);
  double term3=c1l/(c1l+A_12*c2l);
  double term4=c2l*A_21/(c1l*A_21+c2l);
  return term1-term2-term3-term4;
}
double LB::dfg_dc1(double c1g,double c2g){
  double lmd = lmd_21_lmd_11;
  double v1v2= rho2_bar/rho1_bar;
  double v2v1= rho1_bar/rho2_bar;
  double term1=log(c1g)+1;
  double A_12 =v2v1*exp(lmd), A_21=v1v2*exp(-lmd);
  double term2=log(c1g+A_12*c2g);
  double term3=c1g/(c1g+A_12*c2g);
  double term4=c2g*A_21/(c1g*A_21+c2g);
  return term1-term2-term3-term4;
}
double LB::dfl_dc2(double c1l,double c2l){
  double lmd = lmd_21_lmd_11;
  double v1v2= rho2_bar/rho1_bar;
  double v2v1= rho1_bar/rho2_bar;
  double term1=log(c2l)+1;
  double A_12 =v2v1*exp(-lmd), A_21=v1v2*exp(lmd);
  double term2=log(c1l*A_21+c2l);
  double term3=c2l/(c1l*A_21+c2l);
  double term4=c1l*A_12/(c1l+A_12*c2l);
  return term1-term2-term3-term4;
}
double LB::dfg_dc2(double c1g,double c2g){
  double lmd = lmd_21_lmd_11;
  double v1v2= rho2_bar/rho1_bar;
  double v2v1= rho1_bar/rho2_bar;
  double term1=log(c2g)+1;
  double A_12 =v2v1*exp(lmd), A_21=v1v2*exp(-lmd);
  double term2=log(c1g*A_21+c2g);
  double term3=c2g/(c1g*A_21+c2g);
  double term4=c1g*A_12/(c1g+A_12*c2g);
  return term1-term2-term3-term4;
}
double LB::df0_dphi(double phi0,double c1l,double c2l,double c1g,double c2g){
  double dW_dphi = 2.0*phi0*(1.0 - 2.0*phi0)*(1.0 - phi0);
  double dg_dphi = 6.0*phi0*(1.0 - phi0);
  double fg_fl = fg(c1l,c2l)-fl(c1g,c2g);
  return omega_phi*dW_dphi+omega_mix*fg_fl*dg_dphi;
} 
double LB::df0_dc1(double phi0,double c1l,double c2l,double c1g,double c2g){
  double g=phi0*phi0*(3.0 - 2.0*phi0);
  double dfl = dfl_dc1(c1l,c2l);
  double dfg = dfg_dc1(c1g,c2g);
  return omega_mix*(dfl*(1.0-g)-dfg*g);
}
double LB::df0_dc2(double phi0,double c1l,double c2l,double c1g,double c2g){
  double g=phi0*phi0*(3.0 - 2.0*phi0);
  double dfl = dfl_dc2(c1l,c2l);
  double dfg = dfg_dc2(c1g,c2g);
  return omega_mix*(dfl*(1.0-g)-dfg*g);
}
//Chemical potentials
double LB::mu_phi(int ix,int iy,bool UseNew){
  double phi0=phi(ix,iy,UseNew);
  double  c10=c1(ix,iy,UseNew), c20=1.0-c10;
  //double c1l=c1(ix,0,UseNew), c1g=c1(ix,Ly/2,UseNew); //at y=0 should be always liquid and y=Ly/2 always gas.
  //double c2l=c2(ix,0,UseNew), c2g=c2(ix,Ly/2,UseNew); //at y=0 should be always liquid and y=Ly/2 always gas.
  //c1l=(c1l==0)?cl1:c1l; c1g=(c1g==0)?cg1:c1g;
  c10=(c10==0)?cl1:c10;
  //c2l=(c2l==0)?cl2:c2l; c2g=(c2g==0)?cg2:c2g;
  c20=(c20==0)?cl2:c20;
  //cout << c2g << endl;
  return df0_dphi(phi0,c10,c20,c10,c20)-k_phi*Laplacian_phi(ix,iy,UseNew); 
}
double LB::mu_c1(int ix,int iy,bool UseNew){
  double phi0=phi(ix,iy,UseNew);
  double  c10=c1(ix,iy,UseNew), c20=1.0-c10;
  //double c1l=c1(ix,0,UseNew), c1g=c1(ix,Ly/2,UseNew); //at y=0 should be always liquid and y=Ly/2 always gas.
  //double c2l=c2(ix,0,UseNew), c2g=c2(ix,Ly/2,UseNew); //at y=0 should be always liquid and y=Ly/2 always gas.
  //c1l=(c1l==0)?cl1:c1l; c1g=(c1g==0)?cg1:c1g;
  c10=(c10==0)?cl1:c10;
  //c2l=(c2l==0)?cl2:c2l; c2g=(c2g==0)?cg2:c2g;
  c20=(c20==0)?cl2:c20;
  return df0_dc1(phi0,c10,c20,c10,c20)-k_c1*Laplacian_c1(ix,iy,UseNew); 
}
double LB::mutilde_cj(double mu_cj,double p0){
  double g_i = (rho2_bar - rho1_bar)/rho2_bar;
  return mu_cj + g_i * p0;
}
double LB::mu_c2(int ix,int iy,bool UseNew){
  double phi0=phi(ix,iy,UseNew);
  double  c10=c1(ix,iy,UseNew), c20=1.0-c10;
  //double c1l=c1(ix,0,UseNew), c1g=c1(ix,Ly/2,UseNew); //at y=0 should be always liquid and y=Ly/2 always gas.
  //double c2l=c2(ix,0,UseNew), c2g=c2(ix,Ly/2,UseNew); //at y=0 should be always liquid and y=Ly/2 always gas.
  //c1l=(c1l==0)?cl1:c1l; c1g=(c1g==0)?cg1:c1g;
  c10=(c10==0)?cl1:c10;
  //c2l=(c2l==0)?cl2:c2l; c2g=(c2g==0)?cg2:c2g;
  c20=(c20==0)?cl2:c20;
  return df0_dc2(phi0,c10,c20,c10,c20)+k_c2*Laplacian_c1(ix,iy,UseNew); 
}
//Forces and vector fields
double LB::Fsx(int ix, int iy, bool UseNew){
  int i; double sum_phi=0,sum_c1=0;
  int jx, jy;
  for(i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum_phi+=w[i]*V[0][i]*mu_phi(jx,jy,UseNew);
    sum_c1+=w[i]*V[0][i]*mu_c1(jx,jy,UseNew);
  }
  return -U_Cs2*(phi(ix,iy,UseNew)*sum_phi+c1(ix,iy,UseNew)*sum_c1);
}
double LB::Fsy(int ix, int iy, bool UseNew){
  int i; double sum_phi=0,sum_c1=0;
  int jx, jy;
  for(i=0;i<Q;i++){
    jx=(Lx+ix+V[0][i])%Lx;
    //if(ix==0||ix==Lx-1) jx=ix;
    jy=(Ly+iy+V[1][i])%Ly;
    //if(iy==0||iy==Ly-1) jy=iy;
    sum_phi+=w[i]*V[1][i]*mu_phi(jx,jy,UseNew);
    sum_c1+=w[i]*V[1][i]*mu_c1(jx,jy,UseNew);
  }
  return -U_Cs2*(phi(ix,iy,UseNew)*sum_phi+c1(ix,iy,UseNew)*sum_c1);
}
double LB::Jx(int ix,int iy,bool UseNew,double Fx){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[0][i]; else suma+=f[ix][iy][i]*V[0][i];
  return suma+0.5*Cs2*Fx;
}
double LB::Jy(int ix,int iy,bool UseNew,double Fy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[1][i]; else suma+=f[ix][iy][i]*V[1][i];
  return suma+0.5*Cs2*Fy;
}
//Auxiliary fields
double LB::H(int ix,int iy,bool UseNew){
  double mu_phi0 = mu_phi(ix,iy,UseNew);
  double lapl_phi= Laplacian_phi(ix,iy,UseNew);
  return -M_phi*(mu_phi0+k_phi*lapl_phi);
}
double LB::Gamma(double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return w[i]*(1+U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2);
}
//Forcing terms
double LB::Fi(double tau,double Ux0,double Uy0,double gr_x, double gr_y,double Fx,double Fy,int i){
  double VU_x=(V[0][i]-Ux0),VU_y=(V[1][i]-Uy0);
  double Gu=Gamma(Ux0,Uy0,i);
  double Gu_G0=Gamma(Ux0,Uy0,i)-Gamma(0,0,i);
  double VUdotFG=VU_x*Fx+VU_y*Fy;
  double VUdotGr=VU_x*gr_x + VU_y*gr_y;
  double UmU2tau = 1.0 - 1.0/(2.0*tau);
  return w[i]*(UmU2tau*VUdotFG*Gu+Cs2*VUdotGr*Gu_G0);
}
double LB::Gi(double tau,double dt_cjux,double dt_cjuy,int i){
  double DtdotVi=dt_cjux*V[0][i]+dt_cjuy*V[1][i];
  double UmU2tau = 1.0 - 1.0/(2.0*tau);
  return U_Cs2*UmU2tau*w[i]*DtdotVi;
}
double LB::Hi(double dt_phiux,double dt_phiuy,double H,double dH_dt,int i){
  double DtdotVi=dt_phiux*V[0][i]+dt_phiuy*V[1][i];
  double UmU2tau = 1.0 - 1.0/(2.0*tau_h);
  return w[i]*(UmU2tau*U_Cs2*DtdotVi+H+0.5*dH_dt);
}
//Pressure
double LB::p(double Ux0,double Uy0,double gr_x,double gr_y,double rho0,int ix,int iy, bool UseNew){
  int i; double sum;
  double UdotGr = Ux0*gr_x + Uy0*gr_y;
  for(sum=0,i=0;i<Q;i++)
    if(UseNew) sum+=fnew[ix][iy][i]; else sum+=f[ix][iy][i];
  return sum+0.5*Cs2*UdotGr;
}
//Equilibrium functions
double LB::feq(double Ux0,double Uy0,double p0,double rho0,int i){
  double Gu=Gamma(Ux0,Uy0,i);
  double Gu_G0=Gamma(Ux0,Uy0,i)-Gamma(0,0,i);
  return w[i]*p0+rho0*Gu_G0*Cs2;
}
double LB::geq(double Ux0,double Uy0,double cj0,double etaj0,double muj0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i];
  double rho_siu=cj0*w[i]*U_Cs2*UdotVi;
  if(i>0) return w[i]*etaj0*muj0+rho_siu;
  else return etaj0*muj0*(w[0]-1.0)+cj0;
}
double LB::heq(double phi0,double Ux0,double Uy0,double gr_x,double gr_y,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i];
  double GrdotVi=gr_x*V[0][i]+gr_y*V[1][i];
  return phi0*w[i]*(1+U_Cs2*UdotVi)+w[i]*B*GrdotVi;
}
//LB steps
void LB::Collision(int t,double gx,double gy){
  int ix,iy,i; 
  double phi0,c10,c20,rho0,p0,eta0;
  double tau_phi0,Utau,UmUtau;
  double gr_x,gr_y,grphi_x,grphi_y,mu_c10,mu_c20;
  double Fx,Fy,Ux0,Uy0;
  double phiux0,phiuy0,dt_phiux0,dt_phiuy0;
  double c1ux0,c1uy0,dt_c1ux0,dt_c1uy0;
  double c2ux0,c2uy0,dt_c2ux0,dt_c2uy0;
  double H0,dH_dt;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macroscopicas
      phi0=phi(ix,iy,false);  c10=c1(ix,iy,false);  c20=1.0-c10;  rho0=rho(c10,c20);
      eta0=eta(c10,c20);  tau_phi0=tau_phi(eta0,rho0);  Utau=1.0/tau_phi0;  UmUtau=1.0-Utau;
      gr_x=drho_dx(ix,iy,false);  gr_y=drho_dy(ix,iy,false);  grphi_x=dphi_dx(ix,iy,false);  grphi_y=dphi_dy(ix,iy,false);
      mu_c10=mu_c1(ix,iy,false);  mu_c20=mu_c2(ix,iy,false);
      Fx=Fsx(ix,iy,false)+gx*rho0;  Fy=Fsy(ix,iy,false)+gy*rho0;
      Ux0=Jx(ix,iy,false,Fx)*U_Cs2/rho0;  Uy0=Jy(ix,iy,false,Fy)*U_Cs2/rho0;
      p0=p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,false);
      //Computing the time derivatives
      phiux0=Ux0*phi0; dt_phiux0=phiux0 - old_phiux[ix][iy];
      phiuy0=Uy0*phi0; dt_phiuy0=phiuy0 - old_phiuy[ix][iy];
      c1ux0=Ux0*c10; dt_c1ux0=c1ux0 - old_c1ux[ix][iy];
      c1uy0=Uy0*c10; dt_c1uy0=c1uy0 - old_c1uy[ix][iy];
      //c2ux0=Ux0*c20; dt_c2ux0=c2ux0 - old_c2ux[ix][iy];
      //c2uy0=Uy0*c20; dt_c2uy0=c2uy0 - old_c2uy[ix][iy];
      H0=H(ix,iy,false);  dH_dt=H0 - old_H[ix][iy];
      //if(t==49)cout<<gr_y<<endl;
      for(i=0;i<Q;i++){
	      fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(Ux0,Uy0,p0,rho0,i)+Fi(tau_phi0,Ux0,Uy0,gr_x,gr_y,Fx,Fy,i);
	      g1new[ix][iy][i]=UmUtaug1*g1[ix][iy][i]+Utaug1*geq(Ux0,Uy0,c10,eta1_bar,mu_c10,i)+Gi(tau_g1,dt_c1ux0,dt_c1uy0,i);
	      //g2new[ix][iy][i]=UmUtaug2*g2[ix][iy][i]+Utaug2*geq(Ux0,Uy0,c20,eta2_bar,mu_c20,i)+Gi(tau_g2,dt_c2ux0,dt_c2uy0,i);
	      hnew[ix][iy][i]=UmUtauh*h[ix][iy][i]+Utauh*heq(phi0,Ux0,Uy0,grphi_x,grphi_y,i)+Hi(dt_phiux0,dt_phiuy0,H0,dH_dt,i);
      }
      old_phiux[ix][iy]=phiux0; old_phiuy[ix][iy]=phiuy0;
      old_c1ux[ix][iy]=c1ux0; old_c1uy[ix][iy]=c1uy0;
      //old_c2ux[ix][iy]=c2ux0; old_c2uy[ix][iy]=c2uy0;
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
	      //g2[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=g2new[ix][iy][i];
	      h[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=hnew[ix][iy][i];
      }
}
void LB::Init(double Ux0,double Uy0){
  double phi0, c10, c20, rho0;
  double gr_x, gr_y,grphi_x,grphi_y, mu_c10, mu_c20, p0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      phi0 =0.5*tanh(2*(double)(iy-0.4*Ly)/W)-0.5*tanh(2*(double)(iy-0.6*Ly)/W);
      c10 =(cg1-cl1)*phi0+cl1;  c20 =(cg2-cl2)*phi0+cl2;  rho0=rho(c10,c20);
      gr_x=drho_dx(ix,iy,false);  gr_y=drho_dy(ix,iy,false);  grphi_x=dphi_dx(ix,iy,false);  grphi_y=dphi_dy(ix,iy,false);
      mu_c10=mu_c1(ix,iy,false);  mu_c20=mu_c2(ix,iy,false);
      p0=0; old_p[ix][iy]=p0;//p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,false);
      //Initialize the time derivatives!
      old_phiux[ix][iy] = phi0*Ux0; old_phiuy[ix][iy] = phi0*Uy0;
      old_c1ux[ix][iy]=c10*Ux0; old_c1uy[ix][iy]=c10*Uy0;
      //old_c2ux[ix][iy]=c20*Ux0; old_c2uy[ix][iy]=c20*Uy0;
      old_H[ix][iy]=H(ix,iy,false);
      for(int i=0;i<Q;i++){
	      f[ix][iy][i]=feq(Ux0,Uy0,p0,rho0,i);
	      g1[ix][iy][i]=geq(Ux0,Uy0,c10,eta1_bar,mu_c10,i);
	      //g2[ix][iy][i]=geq(Ux0,Uy0,c20,eta2_bar,mu_c20,i);
        h[ix][iy][i]=heq(phi0,Ux0,Uy0,grphi_x,grphi_y,i);
    }
  }
}
void LB::ImposeFields(void){
  int i,ix,iy; double phi0;
  for(ix=0;ix<Lx;ix++){
      //Walls
    	for(i=0;i<Q;i++){
        fnew[ix][0][i]=feq(0,0,0,rho1_bar*cl1+rho2_bar*cl2,i);
        fnew[ix][Ly-1][i]=feq(0,0,0,rho1_bar*cl1+rho2_bar*cl2,i);
        g1new[ix][0][i]=geq(0,0,cl1,eta1_bar,0,i);
        g1new[ix][Ly-1][i]=geq(0,0,cl1,eta1_bar,0,i);
      }
    }
}
void LB::Print(const char * NombreArchivo,double gx,double gy){
  ofstream MiArchivo(NombreArchivo); 
  double phi0,c10,c20,rho0,gr_x,gr_y,grphi_x,grphi_y,p0,Fx,Fy,Ux0,Uy0;
  int ix=Lx/2;
    for(int iy=0;iy<Ly;iy++){
      phi0=phi(ix,iy,true); c10=c1(ix,iy,true); c20=1.0-c10; rho0=rho(c10,c20);
      gr_x=drho_dx(ix,iy,true);  gr_y=drho_dy(ix,iy,true);
      grphi_x=dphi_dx(ix,iy,true);  grphi_y=dphi_dy(ix,iy,true);
      Fx=Fsx(ix,iy,true)+gx*rho0;  Fy=Fsy(ix,iy,true)+gy*rho0;
      Ux0=Jx(ix,iy,true,Fx)*U_Cs2/rho0;  Uy0=Jy(ix,iy,true,Fy)*U_Cs2/rho0;
      p0=p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,true);
      //MiArchivo<<iy<<"\t"<<fnew[ix][iy][0]-feq(Ux0,Uy0,p0,rho0,0)
      //             <<"\t"<<g1new[ix][iy][0]-geq(Ux0,Uy0,c10,eta1_bar,mu_c1(ix,iy,true),0)
      //             <<"\t"<<hnew[ix][iy][0]-heq(phi0,Ux0,Uy0,grphi_x,grphi_y,0)<<"\t"<<endl;
      MiArchivo<<iy<<"\t"<<phi0
                   <<"\t"<<c10
                   <<"\t"<<p0
                   <<"\t"<<Uy0
                   <<"\t"<<mu_phi(ix,iy,true)
                   <<"\t"<<mu_c1(ix,iy,true)<<endl;
    }
  MiArchivo.close();
}


int main(void){
  LB Zhang;
  int t,tmax=50000;
  Zhang.Init(0,0);
  
  for(t=0;t<tmax;t++){
    Zhang.Collision(t,0,0);
    Zhang.Advection();
  }
  
  Zhang.Print("zhang.dat",0,0);

  return 0;
}


