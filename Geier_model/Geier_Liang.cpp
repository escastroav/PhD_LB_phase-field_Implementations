#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1;
const int Ly=100;

const int Q=9;

const double U_Cs2 = 3.0;
const double Cs2 = 1.0/U_Cs2;

class LB{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
  double g[Lx][Ly][Q], gnew[Lx][Ly][Q]; // g[ix][iy][i]
  //Adjustable parameters
  double W;               //interface thickness
  double rho_l,rho_g;     //rho_l liquid density; rho_g gas density;
  double sigma;           //surface tension
  double nu_l,nu_g;       //nu_l liquid viscosity; nu_g gas viscosity;
  double M;               //mobility
  double tau_phi;
  //Auxiliary parameters
  double mu_l,  mu_g;
  double tau_l;
  double tau_g;
  double beta;
  double k;
  //cache variable phiu for dt
  double old_phiux[Lx][Ly];
  double old_phiuy[Lx][Ly];
public:
  LB(double W,double rho_l0,double rho_g0,double nu_l0,double nu_g0,double M0,double sigma0);
  double phi(int ix,int iy,bool UseNew);
  double mu_phi(int ix, int iy, bool UseNew);
  double lambda(double phi);
  double rho(double phi);
  double tau(double phi);
  //Gradient components
  double dphi_dx(int ix, int iy, bool UseNew);
  double dphi_dy(int ix, int iy, bool UseNew);
  double Laplacian_phi(int ix, int iy, bool UseNew);
  double norm_gr(double gr_x, double gr_y);
  //Alternative derivative approach
  //double dphi_dx_local(double phi0,double Ux0,double Uy0,double lmd,double dt_phiux,double dt_phiuy,int ix,int iy,bool UseNew);
  //double dphi_dy_local(double phi0,double Ux0,double Uy0,double lmd,double dt_phiux,double dt_phiuy,int ix,int iy,bool UseNew);
  //Normal components
  double nx(double gr_x, double norm);
  double ny(double gr_y, double norm);
  double dphix_moments(double phi0,double Ux0,double Uy0,int ix,int iy,bool UseNew);
  double dphiy_moments(double phi0,double Ux0,double Uy0,int ix,int iy,bool UseNew);
  double nx_moments(double Ux0,double Uy0,int ix,int iy,bool UseNew);
  double ny_moments(double Ux0,double Uy0,int ix,int iy,bool UseNew);
  //Forces
  double Fsx(double mu, double gr_x);
  double Fsy(double mu, double gr_y);
  double Jx(int ix,int iy,bool UseNew,double Fx);
  double Jy(int ix,int iy,bool UseNew,double Fy);
  double Fi(double nx0,double ny0,double dt_phiux,double dt_phiuy,double lmd,int i,double tau0);
  double Gi(double Ux0,double Uy0,double Fx0,double Fy0,double gr_x,double gr_y,int i,double tau0);
  double p(double Ux0,double Uy0,double gr_x,double gr_y,double rho0,int ix,int iy, bool UseNew);
  double feq(double nx0,double ny0,double phi0,double Ux0,double Uy0,int i);
  double geq(double rho0,double p0,double Ux0,double Uy0,int i);
  void Collision(double gx,double gy);
  void Advection(void);
  void Init(double Ux0,double Uy0);
  void ImposeFields(void);
  void Print(const char * NombreArchivo,double mu_l,double mu_g,double gx,double gy);
};
LB::LB(double W0,double rho_l0,double rho_g0,double nu_l0,double nu_g0,double M0,double sigma0){
  //Setting parameteres up
  W=W0;                         //interface thickness
  rho_l=rho_l0, rho_g=rho_g0;   //rho_l liquid density; rho_g gas density;
  sigma=sigma0;                 //surface tension
  nu_l =nu_l0, nu_g =nu_g0;     //nu_l liquid viscosity; nu_g gas viscosity;
  M=M0;                         //Mobility
  tau_phi=M*U_Cs2+0.5;
  //Auxiliary parameters
  mu_l =nu_l*rho_l, mu_g =nu_g*rho_g;
  tau_l = U_Cs2*nu_l + 0.5;
  tau_g = U_Cs2*nu_g + 0.5;
  beta=12.0*sigma/W;
  k = 1.5*sigma*W;
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
    if(UseNew) suma+=fnew[ix][iy][i]; else suma+=f[ix][iy][i];
  return suma;
}
double LB::lambda(double phi){
  return (1.0 - 4*phi*phi)/W; 
}
double LB::rho(double phi){
  return (phi+0.5)*(rho_l-rho_g) + rho_g;
}
double LB::tau(double phi){
  //return phi*(tau_l - tau_g) + tau_g;
  return (phi<0) ? tau_g : tau_l;
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
double LB::dphix_moments(double phi0,double Ux0,double Uy0,int ix,int iy,bool UseNew){
  double k10=0,k01=0,knorm;
  for(int i=0;i<Q;i++){
    k10+=UseNew ? fnew[ix][iy][i]*(V[0][i]-Ux0) : f[ix][iy][i]*(V[0][i]-Ux0);
    k01+=UseNew ? fnew[ix][iy][i]*(V[1][i]-Uy0) : f[ix][iy][i]*(V[1][i]-Uy0);
  }
  knorm =(k10==0&&k01==0)?1:sqrt(k10*k10+k01*k01);
  return 3*(k10+k10*M*lambda(phi0)/knorm)/tau_phi;
}
double LB::dphiy_moments(double phi0,double Ux0,double Uy0,int ix,int iy,bool UseNew){
  double k10=0,k01=0,knorm;
  for(int i=0;i<Q;i++){
    k10+=UseNew ? fnew[ix][iy][i]*(V[0][i]-Ux0) : f[ix][iy][i]*(V[0][i]-Ux0);
    k01+=UseNew ? fnew[ix][iy][i]*(V[1][i]-Uy0) : f[ix][iy][i]*(V[1][i]-Uy0);
  }
  knorm =(k10==0&&k01==0)?1:sqrt(k10*k10+k01*k01);
  return 3*(k01+k01*M*lambda(phi0)/knorm)/tau_phi;
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
double LB::norm_gr(double gr_x, double gr_y){
  return sqrt(gr_x*gr_x + gr_y*gr_y)+1e-20;//(gr_x == 0 && gr_y == 0) ? 1 : sqrt(gr_x*gr_x + gr_y*gr_y);
}
double LB::nx(double gr_x, double norm){
  return gr_x / norm;
}
double LB::ny(double gr_y, double norm){
  return gr_y / norm;
}
double LB::nx_moments(double Ux0,double Uy0,int ix,int iy,bool UseNew){
  double k10=0,k01=0,knorm;
  for(int i=0;i<Q;i++){
    k10+=UseNew ? fnew[ix][iy][i]*(V[0][i]-Ux0) : f[ix][iy][i]*(V[0][i]-Ux0);
    k01+=UseNew ? fnew[ix][iy][i]*(V[1][i]-Uy0) : f[ix][iy][i]*(V[1][i]-Uy0);
  }
  knorm=(k10==0&&k01==0)?1:sqrt(k10*k10+k01*k01);
  return k10/knorm;
}
double LB::ny_moments(double Ux0,double Uy0,int ix,int iy,bool UseNew){
  double k10=0,k01=0,knorm;
  for(int i=0;i<Q;i++){
    k10+=UseNew ? fnew[ix][iy][i]*(V[0][i]-Ux0) : f[ix][iy][i]*(V[0][i]-Ux0);
    k01+=UseNew ? fnew[ix][iy][i]*(V[1][i]-Uy0) : f[ix][iy][i]*(V[1][i]-Uy0);
  }
  knorm =(k10==0&&k01==0)?1:sqrt(k10*k10+k01*k01);
  return k01/knorm;
}
double LB::mu_phi(int ix, int iy, bool UseNew){
  double phi0=phi(ix,iy,UseNew);
  return 4.0*beta*(phi0+0.5)*(phi0-0.5)*phi0-k*Laplacian_phi(ix,iy,UseNew);
}
double LB::Fsx(double mu, double gr_x){
  return mu*gr_x;
}
double LB::Fsy(double mu, double gr_y){
  return mu*gr_y;
}
double LB::Jx(int ix,int iy,bool UseNew,double Fx){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=gnew[ix][iy][i]*V[0][i]; else suma+=g[ix][iy][i]*V[0][i];
  return suma+0.5*Fx;
}
double LB::Jy(int ix,int iy,bool UseNew,double Fy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=gnew[ix][iy][i]*V[1][i]; else suma+=g[ix][iy][i]*V[1][i];
  return suma+0.5*Fy;
}
double LB::Fi(double nx0,double ny0,double dt_phiux,double dt_phiuy,double lmd,int i,double tau0){
  double DtdotVi=dt_phiux*V[0][i]+dt_phiuy*V[1][i], NdotVi=nx0*V[0][i]+ny0*V[1][i];
  return U_Cs2*(1.0 - 1.0/(2*tau0))*w[i]*(DtdotVi+Cs2*lmd*NdotVi);
}
double LB::Gi(double Ux0,double Uy0,double Fx0,double Fy0,double gr_x,double gr_y,int i,double tau0){
  int vx=V[0][i],  vy=V[1][i];
  double FdotVi=Fx0*vx+Fy0*vy,  UdotVi=Ux0*vx+Uy0*vy, GrdotVi=gr_x*vx+gr_y*vy;
  return U_Cs2*(1.0 - 1.0/(2*tau0))*w[i]*(FdotVi+(rho_l-rho_g)*UdotVi*GrdotVi);
}
double LB::p(double Ux0,double Uy0,double gr_x,double gr_y,double rho0,int ix,int iy, bool UseNew){
  int i; double suma;
  double UdotGr=Ux0*gr_x+Uy0*gr_y;
  double U2=Ux0*Ux0+Uy0*Uy0;
  double cs2_UmW0 = 0.6;//Cs2 / (1.0 - w[0]);
  for(suma=0,i=1;i<Q;i++)
    if(UseNew) suma+=gnew[ix][iy][i]; else suma+=g[ix][iy][i];
  return cs2_UmW0*(suma+0.5*(rho_l-rho_g)*UdotGr-rho0*w[0]*U2*0.5*U_Cs2);
}
double LB::feq(double nx0,double ny0,double phi0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0, NdotVi=nx0*V[0][i]+ny0*V[1][i];
  return phi0*w[i]*(1+U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2)+w[i]*M*U_Cs2*lambda(phi0)*NdotVi;
}
double LB::geq(double rho0,double p0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  double rho_siu=rho0*w[i]*(U_Cs2*UdotVi+U_Cs2*U_Cs2*0.5*UdotVi*UdotVi-U_Cs2*0.5*U2);
  if(i>0) return p0*w[i]*U_Cs2+rho_siu;
  else return p0*(w[0]-1.0)*U_Cs2+rho_siu;

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
      Fx=Fsx(mu0,gr_x)+gx;  Fy=Fsy(mu0,gr_y)+gy;
      Ux0=Jx(ix,iy,false,Fx)/rho0;  Uy0=Jy(ix,iy,false,Fy)/rho0;
      p0=p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,false);
      phiux0=Ux0*phi0; dt_phiux0=phiux0 - old_phiux[ix][iy];
      phiuy0=Uy0*phi0; dt_phiuy0=phiuy0 - old_phiuy[ix][iy];
      for(i=0;i<Q;i++){
	      fnew[ix][iy][i]=UmUtau_phi*f[ix][iy][i]+Utau_phi*feq(nx0,ny0,phi0,Ux0,Uy0,i);//+Fi(nx0,ny0,dt_phiux0,dt_phiuy0,lmd0,i,tau_phi);
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
  double phi0, rho0, gr_x, gr_y, p0,norm0,nx0,ny0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      phi0 =-0.5*tanh(2*(double)(iy-Ly/2)/W);
      rho0=rho(phi0);
      //gr_x=dphi_dx(ix,iy,false);  gr_y=dphi_dy(ix,iy,false);
      gr_x=dphix_moments(phi0,Ux0,Uy0,ix,iy,false);  gr_y=dphiy_moments(phi0,Ux0,Uy0,ix,iy,false);
      //norm0=norm_gr(gr_x,gr_y); nx0=nx(gr_x,norm0); ny0=ny(gr_y,norm0);
      nx0=nx_moments(phi0,Ux0,Uy0,ix,iy,false);
      ny0=ny_moments(phi0,Ux0,Uy0,ix,iy,false);
      p0=p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,false);
      old_phiux[ix][iy] = phi0*Ux0; old_phiuy[ix][iy] = phi0*Uy0;
      for(int i=0;i<Q;i++){
	      f[ix][iy][i]=feq(nx0,ny0,phi0,Ux0,Uy0,i);
	      g[ix][iy][i]=geq(rho0,p0,Ux0,Uy0,i);
      }
  }
}
void LB::ImposeFields(void){
  int i,ix,iy; double phi0;
  for(ix=0;ix<Lx;ix++){
      //Walls
    	for(i=0;i<Q;i++){
        fnew[ix][0][i]=feq(0,0,1,0,0,i);
        fnew[ix][Ly-1][i]=feq(0,0,0,0,0,i);
        gnew[ix][0][i]=geq(rho_l,p(0,0,0,0,rho_l,ix,0,false),0,0,i);
        gnew[ix][Ly-1][i]=geq(rho_g,p(0,0,0,0,rho_g,ix,Ly-1,false),0,0,i);
      }
    }
}
void LB::Print(const char * NombreArchivo,double mu_l,double mu_g,double gx,double gy){
  ofstream MiArchivo(NombreArchivo); 
  double phi0,rho0,mu0,tau0,lmd0,Ux0,Uy0,p0; double Fx,Fy,gr_x,gr_y,gr_xloc,gr_yloc;
  double norm0,ny0,ny0_m;
  double U_l0 = gx*Ly*Ly/(8*mu_l);
  double a = (mu_g-mu_l)/(mu_g+mu_l);
  double b = 2/(mu_g+mu_l);
  double U_g0 = gx*Ly*Ly/(8*mu_g);
  double dt_phiux,dt_phiuy;
  double Ly0=(double)Ly-1.0;
  double Ux_t;
  int ix=Lx/2;
    for(int iy=0;iy<Ly;iy++){
      phi0=phi(ix,iy,true);
      gr_x=dphi_dx(ix,iy,true);  gr_y=dphi_dy(ix,iy,true);
      norm0=norm_gr(gr_x,gr_y); ny0=ny(gr_y,norm0);
      mu0=mu_phi(ix,iy,true);  rho0=rho(phi0);tau0=tau(phi0);lmd0=lambda(phi0);
      Fx=Fsx(mu0,gr_x)+gx;  Fy=Fsy(mu0,gr_y)+gy;
      Ux0=Jx(ix,iy,true,Fx)/rho0;  Uy0=Jy(ix,iy,true,Fy)/rho0;
      ny0_m=ny_moments(Ux0,Uy0,ix,iy,true);
      dt_phiux=phi0*Ux0 - old_phiux[ix][iy];
      dt_phiuy=phi0*Uy0 - old_phiuy[ix][iy];
      Ux_t=(iy<Ly/2) ? U_l0*(-(iy-Ly0/2)*(iy-Ly0/2)*4/(Ly0*Ly0)-(iy-Ly0/2)*2*a/Ly0+mu_l*b) : U_g0*(-(iy-Ly0/2)*(iy-Ly0/2)*4/(Ly0*Ly0)-(iy-Ly0/2)*2*a/Ly0+mu_g*b);
      p0=p(Ux0,Uy0,gr_x,gr_y,rho0,ix,iy,false);
      MiArchivo<<iy<<" "
              <<phi0<<" "
              <<Ux_t<<" "
              <<Ux0<<" "
              <<ny0<<" "
              <<ny0_m<<endl;
    }
  MiArchivo.close();
}


int main(void){
  double W=5;
  double sigma=1.0e-3;             //surface tension
  double rho_l=10, rho_g=1;  //rho_l liquid density; rho_g gas density;
  double nu_l=0.1, nu_g=0.1;
  double mu_l =nu_l*rho_l, mu_g =nu_g*rho_g;
  double M=0.001;
  LB Geier(W,rho_l,rho_g,nu_l,nu_g,M,sigma);
  int t,tmax=100000;
  double Uc=1e-4,g=4*Uc*(mu_l+mu_g)/(Ly*Ly);
  cout << g << endl; 
  Geier.Init(g,0);
  
  for(t=0;t<tmax;t++){
    Geier.Collision(g,0);
    //Geier.ImposeFields();
    Geier.Advection();
  }
  
  Geier.Print("Geier.dat",mu_l,mu_g,g,0);

  return 0;
}


