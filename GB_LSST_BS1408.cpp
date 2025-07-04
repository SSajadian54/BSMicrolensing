#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include "VBBinaryLensingLibrary.h"
using namespace std;


const int    Num=5000;
const double MaxD=20.0;///kpc
const double step=MaxD/(double)Num/1.0;///step in kpc
const double RA=180.0/M_PI;
const double pi=M_PI; 
const double binary_fraction=double(2.0/3.0);
const double velocity=299792458.0;//velosity of light
const double Msun=1.98892*pow(10.,30); //in [kg].
const double Rsun=6.957*pow(10.0,8.0); ///solar radius [meter]
const double KP=  3.08568025*pow(10.,19); // in meter.
const double G=   6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double AU=  1.4960*pow(10.0,11.0);
const double vro_sun=226.0;
const double VSunR =11.1;
const double VSunT =vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ =7.25;
const double year  =365.2425;//days 
const double tE_max=300.0;  
const double Avks=double(8.20922);

///============================ Besancon constant ===========================///

const double Dsun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419,6.2/3.52112,4.0/2.27237,5.8/3.29525,4.9/2.78402, 6.6/3.74991, 3.96/2.24994};

const int M=6;      //filter ugrizy
const int Nl=20000; //filter_LSST.txt
const int No=1500;  // points for stellar orbits
const int N1=29786, N2=28809, N3=4650, N4=134;///CMD_BESANCON, thinD, bulge, thickD, halo 10/08/1402
const int NB=1000;  //Bessel calculation  
const double thre=0.001;//threshold for bessell  


///============================== LSST constant =============================///
const double Tobs=   7.0*30.0;///LSST observational time 10 years
const double cadence=3.0;//in days 
const double sea=double(7.0/12.0);    
const double texp=   30.0;///in second
const double Delta2= 0.005;///systematic errors
const double gama[M]={0.037,0.038,0.039,0.039,0.040,0.040};///a function of sky brightness and airmass in different wavelengths.
const double seeing[M]={0.77,0.73,0.70,0.67,0.65,0.63};///seeing
const double msky[M]={22.9,22.3,21.2,20.5,19.6,18.6};
const double Cm[M]=  {22.92,24.29,24.33,24.20,24.07,23.69};
const double Dci[M]= {0.67,0.21,0.11,0.08,0.05,0.04};
const double km[M]=  {0.451,0.163,0.087,0.065,0.043,0.138};///sky extinction
const double wav[M]= {0.3671, 0.4827, 0.6223, 0.7546, 0.8691, 0.9712};///in[micrometer] u g r i z y
const double sigma[M]={0.022,0.02,0.017,0.017,0.027,0.027};
const double detect[M]={23.4,24.6,24.3,23.6,22.9,21.7};///depth of single visit in ugriz
const double satu[M]=  {15.2,16.3,16.0,15.3,14.6,13.4};///saturation limit of single visit in ugriz
const double M50[M]=   {23.68,24.89,24.43,24.00,24.45,22.60};
const double FWHM[M]=  {1.22087,1.10136,0.993103,0.967076,0.951766,0.936578};//LSST [arcsec]
const double AlAv[M]=  {1.55214, 1.17507, 0.870652, 0.665363, 0.509012, 0.423462};// u g r i z y 

//double x, y, a, b, AlAv; 
//for (int i=0; i<M; ++i) {
//x=double(1.0/wav[i]); 
//y=double(x-1.82);
//a=1.0+0.17699*y-0.50447*y*y-0.02427*y*y*y+0.72085*y*y*y*y+ 0.01979*y*y*y*y*y-0.7753*y*y*y*y*y*y+0.32999*y*y*y*y*y*y*y;
///b=1.41338*y+2.28305*y*y+1.07233*y*y*y-5.38434*y*y*y*y-0.62251*y*y*y*y*y+5.3026*y*y*y*y*y*y-2.09002*y*y*y*y*y*y*y;   
///if(x<=1.1){
//a= 0.574*pow(x, 1.61); 
//b=-0.527*pow(x, 1.61);}
//AlAv=double(a+b/3.1);
//cout<<"wav[i]:   "<<wav[i]<<"\t x: "<<x<<"\t  AlAv:  "<<AlAv<<endl;} 
   
///=============================================================================

struct source{
    int    nums,struc,cl1, cl2;
    double Ds,TET,FI, lat, lon;
    double od_disk, od_ThD, od_bulge, od_halo, opt;///optical  depth
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs, Romins,nstart, nstarti;
    double blend[M], Fluxb[M], magb[M], Ai[M], nsbl[M]; 
    double vs, semi, q, period, ecen, inc, tet;
    double radi1,radi2, mass1, mass2, tef1, tef2, typ1, typ2, age;
    double Mab2[M],Map2[M], Mab1[M], Map1[M]; 
    double xv, yv, zv; 
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
    double deltao;  
};

struct lens{
    double Ml,Dl,vl,Vt, xls, u0, magni;
    double rhomaxl, tE, RE, ksi;
    double ro_star1, ro_star2; 
    int    numl,struc;
};

struct detection{
   int    peak;
   double phi,ksi,weight; 
   double tmin, tmax,t0,tp, dto;
   double x1, y1, x2, y2;  
   double uxs1, uys1, uxs2, uys2; 
   double Astar1, Astar2, Astar0; 
   double ulx, uly, f21;   
   double chi1, chi2, chib, chit, chif, dchiL;
};
struct CMD{  
    double tef_d[N1], age_d[N1],mass_d[N1], Mab_d[M][N1],rs_d[N1],typ_d[N1];  int cl_d[N1];///thin disk
    double tef_b[N2], age_b[N2],mass_b[N2], Mab_b[M][N2],rs_b[N2],typ_b[N2];  int cl_b[N2];/// bulge
    double tef_t[N3], age_t[N3],mass_t[N3], Mab_t[M][N3],rs_t[N3],typ_t[N3];  int cl_t[N3];///thick disk
    double tef_h[N4], age_h[N4],mass_h[N4], Mab_h[M][N4],rs_h[N4],typ_h[N4];  int cl_h[N4];/// halo
};
struct lsst{
    double airm[Nl];
    int  filter[Nl];         
};
struct extinc{
    double dis[100];///distance
    double Extks[100];///ks-band extinction
    double Aks;
    int flag;
};
///===================== FUNCTION ==============================================

int    Extinction(extinc & ex,source & s);
void   read_cmd(CMD & cm);
void   func_source(source & s, CMD & cm, extinc & ex);
void   func_lens(  lens & l, source & s, detection & d);
void   vrel(source & s, lens & l);
void   Disk_model(source & s, int);
void   optical_depth(source & s);
void   Secondary(source & s, CMD & cm);
double ThirdKepler(source & s);  
double Interpol(double ds, extinc & ex);
double errlsst2(double, int, double); 
double Kepler(double , double); 
double Bessel(int n,double x); 
double RandN(double , double);
double RandR(double , double);

///=============================================================================

    time_t _timeNow;
      unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;

///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//
int main()
{

   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   
   VBBinaryLensing vbb;
   vbb.Tol=1.e-4;
   vbb.a1=0.0;  
   vbb.LoadESPLTable("./ESPL.tbl");
  
   CMD cm;  
   source s;
   lens l;
   detection d;
   extinc ex; 
   lsst ls; 
   read_cmd(cm); 
   cout<<"***** read_cmd was done *******"<<endl; 
   
    
    
   char filenam1[40], filenam2[40];  ;
   FILE* fil1; 
   FILE* fil2; 
   fil1=fopen("./files/light/lcF/files/resGBLsstA.dat", "a+");
   fil2=fopen("./files/light/lcF/files/resGBLsstB.dat", "a+");   
   fclose(fil1);  
   fclose(fil2);   
   FILE *model;
   FILE *data1;  

     
     
     
   FILE *lsstf; 
   lsstf=fopen("./files/filterLSST.txt","r");
   if(!lsstf){cout<<"cannot read filterLSST.txt: "<<endl; int yye;  cin>>yye;}
   for(int i=0; i<Nl; ++i){
   fscanf(lsstf,"%lf   %d\n",&ls.airm[i],&ls.filter[i]);  
   if(ls.airm[i]<0.0 or ls.filter[i]>=6 or ls.filter[i]<0){ 
   cout<<"ERROR airmass: "<<ls.airm[i]<<"\t filter: "<<ls.filter[i]<<endl;  int yye; cin>>yye; }}
   cout<<"**** File filterLSST.txt was read ****"<<endl;    
   fclose(lsstf);
   
  
                            
   int    Ndet=0, flag=0, rd, fi, ndata, fl1, fl2, fl3, threef;
   double u1, u2, u0, x0, y0, lonn;
   double prob1,prob2, magnio, errg, erra, tim, timc;
   double delta, magt, magf,  Af, Atot; 
   double A1, A2, Atoto, Atotf, fbl; 
   
///===================== Monte Carlo Simulation ================================
   for(int nsim=3020;  nsim<30000; ++nsim){
   s.lat= RandR(-6.5, -1.5);  
   lonn=  RandR(-10.0,-1.0); 
   if(lonn<=0.0) s.lon=360.0+lonn;
   else          s.lon=lonn;
   s.TET=(360.0-s.lon)/RA;///radian 
   s.FI=  s.lat/RA; 
   Disk_model(s,nsim); 
   ex.flag=-1.0;  
   ex.flag=Extinction(ex,s); 



   do{
   func_source(s,cm,ex);
   func_lens(l,s,d);
   fbl=RandR(0.0 , 1.0);  
   }while(l.tE<=0.5 or l.tE>tE_max or fbl>float(s.blend[2]) ); 
   optical_depth(s);
   
   
   
  
   flag=0; 
   if(nsim>=0){
   flag=1; 
   sprintf(filenam1,"./files/light/lcF/files/%c%d.dat",'M',nsim);
   sprintf(filenam2,"./files/light/lcF/files/%c%d.dat",'L',nsim);
   model=fopen(filenam1,"w"); 
   data1=fopen(filenam2,"w");}
    
    
    
   d.peak=0;      timc=0.0;  d.dchiL=0.0;  ndata=0;     threef=0; 
   d.chi2=0.0;  d.chi1=0.0;  d.chib=0.0;   d.chit=0.0;  d.chif=0.0;   
   rd=int(RandR(0.0,Nl-1));
   fl1=0;  fl2=0;   fl3=0; 
   
 
   for(tim=d.tmin;  tim<d.tmax; tim+=d.dto){ 
   d.phi=double(tim-d.tp)*2.0*pi/s.period;   
   if(s.ecen<0.01) d.ksi=d.phi;//circular
   else            d.ksi=Kepler(d.phi, s.ecen);//elliptic
   x0=s.semi*(cos(d.ksi)-s.ecen);//[m]
   y0=s.semi*sin(d.ksi)*sqrt(1.0-s.ecen*s.ecen);//[m]
   d.x1= -sin(s.inc)*( -y0*sin(s.tet) + x0*cos(s.tet) )*s.mass2/(s.mass1+s.mass2);//[m]  
   d.x2= +sin(s.inc)*( -y0*sin(s.tet) + x0*cos(s.tet) )*s.mass1/(s.mass1+s.mass2);//[m]
   d.y1= +(y0*cos(s.tet)+x0*sin(s.tet))*s.mass2/(s.mass1+s.mass2);//[m] 
   d.y2= -(y0*cos(s.tet)+x0*sin(s.tet))*s.mass1/(s.mass1+s.mass2);//[m] 
   d.ulx=-l.u0*sin(l.ksi) +(tim-d.t0)/l.tE *cos(l.ksi);  
   d.uly=+l.u0*cos(l.ksi) +(tim-d.t0)/l.tE *sin(l.ksi); 

   d.x1= d.x1*l.xls/l.RE; 
   d.y1= d.y1*l.xls/l.RE; 
   d.x2= d.x2*l.xls/l.RE; 
   d.y2= d.y2*l.xls/l.RE;
   d.uxs1=d.ulx - d.x1;  
   d.uys1=d.uly - d.y1;  
   d.uxs2=d.ulx - d.x2;  
   d.uys2=d.uly - d.y2;
     
   u0=sqrt(d.ulx*d.ulx   + d.uly*d.uly  );
   u1=sqrt(d.uxs1*d.uxs1 + d.uys1*d.uys1); 
   u2=sqrt(d.uxs2*d.uxs2 + d.uys2*d.uys2); 
   
   d.Astar0=vbb.ESPLMag2(u0, 0.5*(l.ro_star1+l.ro_star2) );  
   d.Astar1=vbb.ESPLMag2(u1, l.ro_star1);  
   d.Astar2=vbb.ESPLMag2(u2, l.ro_star2);
   
   magt=-2.5*log10(pow(10.0,-0.4*s.Map1[2])*(d.Astar1-1.0)+ pow(10.0,-0.4*s.Map2[2])*(d.Astar2-1.0)+s.Fluxb[2]);//r-band
   Atot=double(pow(10.0,-0.4*s.Map1[2])*(d.Astar1-1.0)+pow(10.0,-0.4*s.Map2[2])*(d.Astar2-1.0)+s.Fluxb[2])/s.Fluxb[2];//r-band
   A1=  double(pow(10.0,-0.4*s.Map1[2])*(d.Astar1-1.0)+s.Fluxb[2])/s.Fluxb[2];  
   A2=  double(pow(10.0,-0.4*s.Map2[2])*(d.Astar2-1.0)+s.Fluxb[2])/s.Fluxb[2];  
   Af=  double(pow(10.0,-0.4*s.Map1[2])*(d.Astar0-1.0)+pow(10.0,-0.4*s.Map2[2])*(d.Astar0-1.0)+s.Fluxb[2])/s.Fluxb[2];
   
   if(flag>0) 
   fprintf(model,"%.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.6lf  %.6lf  %.6lf\n", 
   tim,d.ulx,d.uly, d.uxs1,d.uys1, d.uxs2,d.uys2, A1, A2, Atot);//10
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    if(tim>0.0 and tim<Tobs){  
    timc+=d.dto; 
    if(timc>cadence){
    timc=timc-cadence; 
    prob1=RandR(0.0,100.0);///probability of bad weather
    prob2=RandR(0.0,100.0);///probability of uniform observation
    if(fabs(tim-d.t0)<double(3.0*cadence)) d.peak=1; 
    if(prob1>17.0 and prob2>=10.0){
    
    rd+=1; 
    if(rd>=(Nl-1)) rd=0;
    fi=int(ls.filter[rd]);
    Atotf=double(pow(10.0,-0.4*s.Map1[fi])*(d.Astar1-1.0)+pow(10.0,-0.4*s.Map2[fi])*(d.Astar2-1.0)+s.Fluxb[fi])/s.Fluxb[fi]; 
    magt =s.magb[fi]-2.5*log10(Atotf);  
    if(magt>=satu[fi] and magt<=detect[fi]){
    errg=errlsst2(magt, fi, double(ls.airm[rd]));
    erra= Atotf*fabs( pow(10.0 , -0.4*errg)-1.0 ); 
    
    
    Atoto=Atot+RandN(erra,2.5); 
    d.chib +=(Atoto-  1.)*(Atoto-  1.)/(erra*erra);//baseline
    d.chi1 +=(Atoto-  A1)*(Atoto-  A1)/(erra*erra);//magnification of the first star
    d.chi2 +=(Atoto-  A2)*(Atoto-  A2)/(erra*erra);//magnification of the second star
    d.chit +=(Atoto-Atot)*(Atoto-Atot)/(erra*erra);//lensing from binary source star(real_model)  
    d.chif +=(Atoto-  Af)*(Atoto-  Af)/(erra*erra);//simple light curve due to CM of binary stars  
    ndata+=1;
    fl1=fl2;
    fl2=fl3;
    if(double(Atoto-1.0)>double(3.0*erra)) fl3=1;
    else fl3=0; 
    if(ndata>=3 and double(fl1+fl2+fl3*1.0)>2.0)  threef=1; 
    if(flag>0) fprintf(data1,"%.4lf  %.7lf  %.7lf  %d\n",tim, Atoto, erra, fi);//4 
    }}}}}
    if(flag>0){fclose(data1);  fclose(model);} 
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH    
   d.dchiL=fabs(d.chib-d.chit);
   if(d.peak>0 and d.dchiL>300 and threef>0)  Ndet+=1;  
   fil1=fopen("./files/light/lcF/files/resGBLsstA.dat", "a+");
   fil2=fopen("./files/light/lcF/files/resGBLsstB.dat", "a+");   
   fprintf(fil1,
   "%d  %d  %.3lf  %.3lf "//4
   "%d  %.4lf  %.4lf  %.4lf   %.4lf  %.4lf  %.3lf  %.3lf  %.4lf  %.4lf  %d  %d  %.2lf  %.2lf  "//18
   "%.8lf  %.5lf  %.5lf  %.5lf  %.5lf   %.5lf  %.5lf  " //25
   "%d  %.4lf  %.4lf  %.5lf  %.5lf  %.5lf  %.6lf  %.9lf  %.9lf   %.6lf  %.5lf  "//36
   "%.5lf  %.5lf  %.5lf  %.5lf  %.5lf   %d  %.2lf  %e  %.1lf  %.1lf  %.1lf %.1lf  %.1lf  %.6lf %d  %d\n", //52
    nsim, Ndet, s.lon, s.lat,//4 
    s.struc, s.Ds, s.vs, s.mass1, s.mass2, s.age, s.tef1, s.tef2, s.radi1, s.radi2, s.cl1, s.cl2, s.typ1, s.typ2,//18
    s.semi/Rsun, s.period, s.ecen,  s.q, s.inc*RA, s.tet*RA, s.opt*1.0e6,//25
    l.struc, l.Dl, l.vl, l.Ml, l.Vt, l.RE/AU, l.tE, l.ro_star1, l.ro_star2, l.u0, l.ksi*RA,//36
    d.t0, d.tmin, d.tmax, d.tp, d.dto, d.peak, d.dchiL, d.f21, d.chib, d.chi1,d.chi2,d.chit,d.chif,s.blend[2],ndata,threef);//52
    for(int i=0; i<M; ++i){
    fprintf(fil2,"%d  %d  %.4lf  %.1lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.3lf\n",
    nsim,i,s.magb[i],s.nsbl[i],s.blend[i],s.Mab1[i],s.Map1[i],s.Mab2[i],s.Map2[i],s.Ai[i]);}
    fclose(fil1);  
    fclose(fil2);   
    
    cout<<"============================================================="<<endl;
    cout<<"nsim:  "<<nsim<<"\t Ndet:  "<<Ndet<<endl;
    cout<<"Galactic longitude:  "<<s.lon<<"\t Galactic Latitude:  "<<s.lat<<"\t weight:  "<<d.weight<<endl;
    cout<<"mass1:  "<<s.mass1<<"\t mass2:  "<<s.mass2<<"\t age:  "<<s.age<<endl;
    cout<<"cl1:  "<<s.cl1<<"\t cl2:  "<<s.cl2<<"\t typ1:  "<<s.typ1<<"\t tp2:   "<<s.typ2<<endl;
    cout<<"tef1:  "<<s.tef1<<"\t tef2:  "<<s.tef2<<"\t radi1:  "<<s.radi1<<"\t radi2:  "<<s.radi2<<endl;
    cout<<"semi[Rsun]:  "<<s.semi/Rsun<<"\t period(days): "<<s.period<<"\t ecen:  "<<s.ecen<<endl;
    cout<<"inclination:  "<<s.inc*RA<<"\t tet:  "<<s.tet*RA<<"\t Ds:  "<<s.Ds<<endl;
    cout<<"DL:   "<<l.Dl<<"\t tE:  "<<l.tE<<"\t u0:  "<<l.u0<<endl;
    cout<<"f21:  "<<d.f21<<"\t peak:  "<<d.peak<<"\t Dchi:  "<<d.dchiL<<endl; 
    cout<<"============================================================="<<endl;
    }
    fclose(_randStream);
    return(0);
}
////==========================================================================
double RandR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return( p*(up-down)+down );
}
///=======================================================================//
///                                                                       //
///                  Photometry error LSST                                //
///                                                                       //
///=======================================================================//
double errlsst2(double maga, int fi, double airm) 
{  
    double seei, mfive, x, Delta1;    
    seei=double(seeing[fi]+RandN(0.2,1));
    if(seei<=0.4) seei=0.4; 
    mfive=Cm[fi]+0.0+0.5*(msky[fi]-21.0)+2.5*log10(0.7/seei)+1.25*log10(texp/30.0)-km[fi]*(airm-1.0) ;
    x=pow(10.0,0.4*(maga-mfive));
    Delta1=sqrt(fabs((0.04-gama[fi])*x+gama[fi]*x*x)); 
    if(Delta1<0.0001)   Delta1=0.0001; 
    if(mfive<0.0 or mfive>40.0  or Delta1<0.00001 or Delta1>0.9) {
    cout<<"Error(LSST_phot_acuracy):  mfive:    "<<mfive<<"\t Delta1:  "<<Delta1<<"\t maga:  "<<maga<<endl;
    cout<<"msky[fi]:  "<<msky[fi]<<"\t seei:  "<<seei<<"\t airm:   "<<airm<<"\t  fi: "<<fi<<endl;
    int uue; cin>>uue;}   
    return( sqrt(Delta2*Delta2+Delta1*Delta1) ); 
}
///==============================================================//6
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*Msun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opt=0.0;
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC; }
    s.opt= fabs(s.od_disk+s.od_ThD+s.od_bulge+s.od_halo );///total
    //cout<<"total_opticalD: "<<s.opt<<"\t od_disk: "<<s.od_disk<<endl;
    //cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///==============================================================//
///                                                              //
///                  Kepler Equation                             //
///                                                              //
///==============================================================//
double Kepler(double phi, double ecen)
{
    double ksi=0;   
    double term, term0;  
    phi=double(phi*RA); 
    while(phi>360.0) phi=phi-360.0; 
    while(phi<0.0)   phi=phi+360.0;      
    if(phi>180)      phi=double(phi-360.0);
    if(phi<-181.0 or phi>181.0){ 
    cout<<"Error :  Phi:  "<<phi<<"\t ecent:  "<<ecen<<endl;   int yye;  cin>>yye;}
    phi=double(phi/RA);
    ksi=phi; 
    for(int i=1; i<NB; ++i){
    term= Bessel(i,i*ecen)*sin(i*phi)*2.0/i;  
    ksi+=term; 
    if(i==1)  term0=fabs(term); 
    if(fabs(term)<double(thre*term0) and i>5)  break;}        
    return(ksi); 
}       
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
double Bessel(int n,double x)
{
    double j1=0.00000001,tet;
    int kmax=10000;
    for(int k=0; k<kmax; ++k){
    tet=double(k*M_PI/kmax);
    j1+=double(M_PI/kmax/1.0)*cos(n*tet-x*sin(tet)); }
    return(j1/M_PI);
}    
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
///==============================================================//
///                                                              //
///                  func_source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD  &  cm, extinc & ex)
{
    int    struc, nums, num;
    double rho,rf, Alv, Avk, Akv;
    double Ds,Ai[M],Av;
    double Map[M],Mab[M];        
    double test,frand;
    double maxnb=0.0;
    double emax; 
    
    
    
    for(int i=0; i<M; ++i){
    s.Fluxb[i]=0.0;
    s.nsbl[i]= fabs(s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0));
    s.nsbl[i]+=RandN(sqrt(s.nsbl[i]),2.0);
    if(s.nsbl[i]<=2.0)   s.nsbl[i]=2.0;
    if(s.nsbl[i]>maxnb)  maxnb=s.nsbl[i];}
    
   
  
    for(int k=1; k<=int(maxnb+0.000000034756346); ++k){  
    do{
    nums=int(RandR(5.0 , Num-5.0));
    rho=RandR(s.Romins , s.Romaxs); 
    Ds=double(nums*step);
    }while(rho>s.Rostari[nums] or Ds<0.0 or Ds>MaxD);///distance larger than 50.0
    //if(Ds>MaxD or Ds<0.0){
    //cout<<"ERROR (1): Ds: "<<Ds<<"\t MaxD: "<<MaxD<<"\t step: "<<step<<"\t num: "<<nums<<endl; int yye; cin>>yye;}
    //cout<<"Ds:  "<<Ds<<endl;
    
    
    rf=RandR(0.0,s.Rostar0[nums]); 
     if (rf<=  s.rho_disk[nums])                    struc=0;///thin disk
else if (rf<=( s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
else if (rf<=( s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
else if (rf<=(s .rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums])) struc=3;///halo


    if(struc==0){///thin disk
    num=int(RandR(0.0 , N1-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_d[i][num];}
    if(k==1){ 
    s.mass1=cm.mass_d[num]; 
    s.radi1=  cm.rs_d[num]; 
    s.age =  cm.age_d[num]; 
    s.tef1 = cm.tef_d[num];
    s.cl1 =   cm.cl_d[num];
    s.typ1=  cm.typ_d[num];}}


    if(struc==1){///bulge
    num=int(RandR(0.0,N2-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_b[i][num];}
    if(k==1){ 
    s.mass1= cm.mass_b[num];
    s.radi1=   cm.rs_b[num]; 
    s.age=    cm.age_b[num]; 
    s.tef1=   cm.tef_b[num];
    s.cl1=     cm.cl_b[num];
    s.typ1=   cm.typ_b[num];}}


    if(struc==2){///thick disk
    num=int(RandR(0.0,N3-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_t[i][num];}
    if(k==1){ 
    s.mass1=cm.mass_t[num];
    s.radi1=  cm.rs_t[num]; 
    s.age=   cm.age_t[num]; 
    s.tef1=  cm.tef_t[num];
    s.cl1=    cm.cl_t[num];
    s.typ1=  cm.typ_t[num];}}


    if(struc==3){///stellar halo
    num=int(RandR(0.0,N4-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_h[i][num];}
    if(k==1){ 
    s.mass1=cm.mass_h[num];
    s.radi1=  cm.rs_h[num]; 
    s.age=   cm.age_h[num]; 
    s.tef1=  cm.tef_h[num];
    s.cl1=    cm.cl_h[num]; 
    s.typ1=  cm.typ_h[num];}}
    
    
    if(ex.flag>0)  Av=double(Interpol(Ds,ex)*Avks);
    else           Av=double(0.7*Ds); 
    if(Av<0.0)     Av=0.0;
    for(int i=0; i<M; ++i){
    Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i],1.5);
    if(Ai[i]<0.0) Ai[i]=0.0;   
    Map[i]=Mab[i]+5.0*log10(Ds*100.0)+Ai[i];
    if(s.nsbl[i]>=k) s.Fluxb[i]+=pow(10.0,-0.4*Map[i]);}
 
    
    
    if(k==1){
    s.struc=struc;   
    s.Ds=Ds;  
    s.nums=nums;


    do{///https://articles.adsabs.harvard.edu/pdf/1991A%26A...248..485D  
    s.q= RandR(0.0 , 1.0);
    test=RandR(0.0 , 1.0);  
    frand=exp(-(s.q-0.23)*(s.q-0.23)/(2.0*0.42*0.42));
    }while(test>frand);
    
      
    s.mass2=double(s.q*s.mass1); 
    Secondary(s, cm); 
    s.q=s.mass2/s.mass1; 
    
    s.semi=RandR(log10(2.0) , log10(10000.0) );//in the radius of the primary lens log-uniform
    s.semi=pow(10.0,s.semi); 
    s.semi=double(s.semi*s.radi1*Rsun);//[m]   
    
    s.period = ThirdKepler(s);//[days]
    emax=double(0.8-8.0*exp(-pow(6.0*s.period,0.35)));//Fig 7
    if(emax<0.0)  emax=0.0; 
    if(emax>1.0)  emax=1.0;
    if(emax<0.0 or emax>1.0){cout<<"Error emax:  "<<emax<<"\t period:  "<<s.period<<endl; int yyw; cin>> yyw; } 
    s.ecen=RandR(0.0, emax);
    s.inc=RandR(0.1,89.9)/RA; 
    s.tet=RandR(0.1, 359.9)/RA;
    
    for(int i=0; i<M; ++i){
    s.Ai[i]=Ai[i];   
    s.Map1[i]=Map[i];   
    s.Mab1[i]=Mab[i];
    s.Map2[i]=s.Mab2[i]+5.0*log10(s.Ds*100.0)+Ai[i];
    s.Fluxb[i]+=pow(10.0,-0.4*s.Map2[i]);} 
    k+=1;}    
    }///loop over the stars 


    for(int i=0; i<M; ++i){
    s.magb[i]=-2.5*log10(s.Fluxb[i]);
    s.blend[i]=double(pow(10.0,-0.4*s.Map1[i])+pow(10.0,-0.4*s.Map2[i]))/s.Fluxb[i];
    if(s.Fluxb[i]<=0.0 or int(s.nsbl[i])<0.0 or s.nsbl[i]==0.0 or s.blend[i]>1.00001 or s.blend[i]<=0.0 or 
    (s.nsbl[i]==2.0 and s.blend[i]<1.0) or fabs(s.Map1[0]-s.Mab1[0]-s.Ai[0]-s.Map1[1]+s.Mab1[1]+s.Ai[1])>0.1 or 
    s.Ds<0.0 or s.Ds>MaxD or s.ecen>1.0 or s.ecen<0.0){ 
    cout<<"nsbl[i]: "<<s.nsbl[M]<<"\t Ds:  "<<s.Ds<<endl; 
    cout<<"BIg ERROR nsbl: "<<s.nsbl[i]<<"\t nsbl[i]: "<<s.nsbl[i]<<"\t blend: "<<s.blend[i]<<endl; 
    int uue; cin>>uue;}}
    //cout<<"********** End of Func_source *************"<<endl;
}
///==============================================================//
///                                                              //
///                  Secondary source stars                      //
///                                                              //
///==============================================================//
void Secondary(source & s, CMD & cm){
   int count=0; 
   int kj=-1; 
   
   
     
   if(s.struc==0){///thin disk
   if(s.mass2<=cm.mass_d[0]) kj=0; 
   else if(s.mass2>=cm.mass_d[N1-1]) kj=int(N1-1); 
   else{ 
   for(int i=0; i<int(N1-1); ++i){
   if( float((s.mass2-cm.mass_d[i])*(s.mass2-cm.mass_d[i+1]))<=0.0){kj=i;  i=N1;}}}
   if(kj<0){cout<<"Error(1) kj<0: "<<kj<<"\t mass2:  "<<s.mass2<<endl;  int yye;  cin>>yye; }
   count=kj; 
   
   if(s.mass2<0.0 or fabs(s.mass2-cm.mass_d[count])>0.2 or kj<0.0 or s.struc<0 or s.struc!=0){
   cout<<"mass2:  "<<s.mass2<<"\t mass_d:  "<<cm.mass_d[count]<<endl; 
   cout<<"strucs:  "<<s.struc<<"\t kj:  "<<kj<<"\t count:  "<<count<<endl;   int uue; cin>>uue;}
   s.mass2=   cm.mass_d[count];  
   s.tef2=     cm.tef_d[count];
   s.radi2=     cm.rs_d[count];  
   s.cl2=       cm.cl_d[count];
   s.typ2=     cm.typ_d[count];
   for(int i=0; i<M; ++i) s.Mab2[i]=cm.Mab_d[i][count];}
   
   
   
  
   kj=-1; 
   if(s.struc==1){///Bulge
   if(s.mass2<=cm.mass_b[0]) kj=0; 
   else if(s.mass2>=cm.mass_b[N2-1]) kj=int(N2-1); 
   else{ 
   for(int i=0; i<int(N2-1); ++i){
   if(float((s.mass2-cm.mass_b[i])*(s.mass2-cm.mass_b[i+1]))<=0.0) {kj=i;  i=N2;}}}
   if(kj<0){cout<<"Error(2) kj<0: "<<kj<<"\t mass2:  "<<s.mass2<<endl;  int yye;  cin>>yye;}
   count=kj; 
   if(s.mass2<0.0 or fabs(s.mass2-cm.mass_b[count])>0.2 or kj<0.0 or s.struc!=1){
   cout<<"mass2:  "<<s.mass2<<"\t mass_b:  "<<cm.mass_b[count]<<endl; 
   cout<<"strucs:  "<<s.struc<<"\t kj:  "<<kj<<"\t count:  "<<count<<endl; int uue; cin>>uue;}
   s.mass2=   cm.mass_b[count];  
   s.radi2=     cm.rs_b[count];
   s.tef2=     cm.tef_b[count];  
   s.cl2=       cm.cl_b[count];
   s.typ2=     cm.typ_b[count];
   for(int i=0; i<M; ++i)  s.Mab2[i]=cm.Mab_b[i][count];}
   
   
   kj=-1; 
   if(s.struc==2){//Thick disk
   if(s.mass2<=cm.mass_t[0])         kj=0; 
   else if(s.mass2>=cm.mass_t[N3-1]) kj=int(N3-1); 
   else{ 
   for(int i=0; i<int(N3-1); ++i){
   if(float((s.mass2-cm.mass_t[i])*(s.mass2-cm.mass_t[i+1]))<=0.0){kj=i;  i=N3;}}}
   if(kj<0){cout<<"Error(3) kj<0: "<<kj<<"\t mass2:  "<<s.mass2<<endl;  int yye;  cin>>yye;}
   count=kj; 
   if(s.mass2<0.0 or fabs(s.mass2-cm.mass_t[count])>0.2 or kj<0.0 or s.struc!=2){
   cout<<"mass2:  "<<s.mass2<<"\t mass_t:  "<<cm.mass_t[count]<<endl; 
   cout<<"strucs:  "<<s.struc<<"\t kj:  "<<kj<<"\t count:  "<<count<<endl;   int uue; cin>>uue;  }
   s.mass2=   cm.mass_t[count];  
   s.radi2=     cm.rs_t[count];
   s.tef2=     cm.tef_t[count];  
   s.cl2=       cm.cl_t[count];
   s.typ2=     cm.typ_t[count];
   for(int i=0; i<M; ++i) s.Mab2[i]=cm.Mab_t[i][count];}
   
   
   
   
   kj=-1;  
   if(s.struc==3){///Stellar halo
   if(s.mass2<=cm.mass_h[0]) kj=0; 
   else if(s.mass2>=cm.mass_h[N4-1]) kj=int(N4-1); 
   else{ 
   for(int i=0; i<int(N4-1); ++i){
   if(float((s.mass2-cm.mass_h[i])*(s.mass2-cm.mass_h[i+1]))<=0.0){kj=i;  i=N4;}}}
   if(kj<0){cout<<"Error(4) kj<0: "<<kj<<"\t mass2:  "<<s.mass2<<endl;  int yye;  cin>>yye; }
   count=kj; 
   if(s.mass2<0.0 or fabs(s.mass2-cm.mass_h[count])>0.2 or kj<0 or s.struc!=3){
   cout<<"mass2:  "<<s.mass2<<"\t mass_h:  "<<cm.mass_h[count]<<endl; 
   cout<<"strucs:  "<<s.struc<<"\t kj:  "<<kj<<"\t count:  "<<count<<endl;   int uue; cin>>uue;  }
   s.mass2=   cm.mass_h[count];  
   s.radi2=     cm.rs_h[count];
   s.tef2=     cm.tef_h[count];  
   s.cl2=       cm.cl_h[count];
   s.typ2=     cm.typ_h[count];
   for(int i=0; i<M; ++i) s.Mab2[i]=cm.Mab_h[i][count];}
   //cout<<"mass2:  "<<s.mass2<<"\t radi2:  "<<s.radi2<<"\t tef2:  "<<s.tef2<<endl;  
   //cout<<"cl2:  "<<s.cl2<<"\t type2:  "<<s.typ2<<endl;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double ThirdKepler(source & s){
    return( sqrt(4.0*pi*pi*s.semi*s.semi*s.semi/(G*Msun*(s.mass1+s.mass2)) )/(24.0*3600.0));///days  
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts              //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s, detection & d)
{
    double f,test, tt;
    double rholens[s.nums+2]={0.0};
    l.rhomaxl=0.0;
    for(int k=1; k<int(s.nums-1); ++k){
    rholens[k]=0.0;
    l.Dl= double(k*step);
    l.xls=double(l.Dl/s.Ds);
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl)  l.rhomaxl=rholens[k];}



    do{
    l.numl =int(RandR(1.0,s.nums-2.0));
    test =RandR(0.0,l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){
    cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;  int ue; cin>>ue;}
    }while(test>rholens[l.numl]);
    l.Dl=double(l.numl*step);
    
   
 double randflag=RandR(0.0,s.Rostar0[l.numl]);
     if (randflag<=fabs(s.rho_disk[l.numl]) )  l.struc=0;//disk
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]))  l.struc=1;//bulge
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]))  l.struc=2;//thick
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]))  l.struc=3;//halo
else {cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}
     //cout<<"Dl:  "<<l.Dl<<"\t struc_lens :  "<<l.struc<<endl;


  if(l.struc==0){///thin disk
  do{
  l.Ml=RandR(0.08,4.5);
  test=RandR(0.0,57.0);
  if(l.Ml<=1.0) f=pow(l.Ml,-1.6);
  if(l.Ml>1.0)  f=pow(l.Ml,-3.0);
  }while(test>f);}




  if(l.struc==1){///Galactic bulge
  do{
  l.Ml=RandR(0.08,1.4);
  test=RandR(0.3,378.3);
  f=pow(l.Ml,-2.35);
  }while(test>f);}


  if(l.struc==2){///thick disk
  do{
  l.Ml=RandR(0.08,1.4);
  test=RandR(0.8,3.8);
  f=pow(l.Ml,-0.5);
  }while(test>f);}


  if(l.struc==3){//stellar halo
  do{
  l.Ml=RandR(0.08,0.8);
  test=RandR(0.8, 3.8);
  f=pow(l.Ml,-0.5);
  }while(test>f);}


  l.ksi=RandR(0.0,2.0*pi);  
  l.xls=l.Dl/s.Ds;
  l.RE=sqrt(4.0*G*l.Ml*Msun*s.Ds*KP)/velocity;
  l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
  vrel(s , l);
  l.tE=fabs(l.RE/(l.Vt*1000.0*3600.0*24.0));///in day
  l.ro_star1=s.radi1*Rsun*l.xls/l.RE; 
  l.ro_star2=s.radi2*Rsun*l.xls/l.RE; 
  l.u0=RandR(0.0,1.0);
  if(l.u0==0.0) l.u0=1.0e-50; 
  //cout<<"RE(AU):  "<<l.RE/AU<<"\t tE(days):  "<<l.tE<<"\t Vt:  "<<l.Vt<<endl;
  //cout<<"ro_*(1):  "<<l.ro_star1<<"\t ro_*(2):  "<<l.ro_star2<<"\t l.u0:  "<<l.u0<<endl;
 
  if(l.tE<0.0 or l.tE==0.0 or l.RE<0.0 or  l.Dl>s.Ds or l.xls>=1.0 or l.Dl<0.0 or s.Ds<0.0 or s.Ds>MaxD or 
  l.Ml<0.06 or l.Ml>5.0){ 
  cout<<"Dl:  "<<l.Dl<<"\t Ds:   "<<s.Ds<<"\t xls:  "<<l.xls<<"\t Ml:   "<<l.Ml<<endl;
  cout<<"BIG ERROR te: "<<l.tE<<"\t RE: "<<l.RE<<"\t V_rel: "<<l.Vt<<endl; int iie; cin>>iie;} 
 
 
 
  d.t0=RandR(10.0,Tobs-10.0);
  d.tmin=d.t0-3.5*l.tE;
  d.tmax=d.t0+3.5*l.tE;
  d.dto=(d.tmax-d.tmin)/No;
  if(d.dto>cadence){d.dto=0.1;}
  //cout<<"Error cadence"<<cadence<<"\t dto:  "<<d.dto<<"\t tE:  "<<l.tE<<endl;   int uus; cin>>uus;  } 
  d.tp=d.t0+RandR(-s.period*0.5 , s.period*0.5);     
  d.f21=pow(10.0,-0.4*(s.Map2[2]-s.Map1[2])); ///r_LSST
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm)
{
    double metal;  
    FILE *fil1; 
    
    fil1=fopen("./files/CMDBESLSST09081402/Disk_Besan_LSSTs.dat","r"); //sorted by mass 8/4/1404
    for(int i=0; i<N1; ++i){
    fscanf(fil1,"%lf %lf  %lf  %lf %lf %d  %lf  %lf %lf %lf  %lf  %lf %lf\n",
    &cm.mass_d[i],   &cm.age_d[i],    &metal,         &cm.tef_d[i],   &cm.rs_d[i],    &cm.cl_d[i], &cm.typ_d[i], 
    &cm.Mab_d[0][i], &cm.Mab_d[1][i], &cm.Mab_d[2][i],&cm.Mab_d[3][i],&cm.Mab_d[4][i],&cm.Mab_d[5][i]);}
    fclose(fil1);  
    
    
    fil1=fopen("./files/CMDBESLSST09081402/Bulge_Besan_LSSTs.dat","r"); 
    for(int i=0; i<N2; ++i){
    fscanf(fil1,"%lf %lf  %lf  %lf %lf %d  %lf  %lf %lf %lf  %lf  %lf %lf\n",
    &cm.mass_b[i], &cm.age_b[i], &metal, &cm.tef_b[i], &cm.rs_b[i], &cm.cl_b[i], &cm.typ_b[i], 
    &cm.Mab_b[0][i], &cm.Mab_b[1][i], &cm.Mab_b[2][i], &cm.Mab_b[3][i], &cm.Mab_b[4][i], &cm.Mab_b[5][i]);}
    fclose(fil1);  
   
    
    fil1=fopen("./files/CMDBESLSST09081402/DiskTh_Besan_LSSTs.dat","r"); 
    for(int i=0; i<N3; ++i){
    fscanf(fil1,"%lf %lf  %lf  %lf %lf %d  %lf  %lf %lf %lf  %lf  %lf %lf\n",
    &cm.mass_t[i], &cm.age_t[i], &metal, &cm.tef_t[i], &cm.rs_t[i], &cm.cl_t[i], &cm.typ_t[i], 
    &cm.Mab_t[0][i], &cm.Mab_t[1][i], &cm.Mab_t[2][i], &cm.Mab_t[3][i], &cm.Mab_t[4][i], &cm.Mab_t[5][i]);}
    fclose(fil1);  
    
 
    fil1=fopen("./files/CMDBESLSST09081402/Halo_Besan_LSSTs.dat","r"); 
    for(int i=0; i<N4; ++i){
    fscanf(fil1,"%lf %lf  %lf  %lf %lf %d  %lf  %lf %lf %lf  %lf  %lf %lf\n",
    &cm.mass_h[i], &cm.age_h[i], &metal, &cm.tef_h[i], &cm.rs_h[i], &cm.cl_h[i], &cm.typ_h[i], 
    &cm.Mab_h[0][i], &cm.Mab_h[1][i], &cm.Mab_h[2][i], &cm.Mab_h[3][i], &cm.Mab_h[4][i], &cm.Mab_h[5][i]);}
    fclose(fil1);  
    
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s, int numt)
{
   int flagf=0;
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nnf=0.4/0.8;
   double mBarre;//stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   double fd=1.0;//see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;//0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;//No limitation 
   double Rdd=2.17;//2.53;///2.17;
   double Rhh=1.33;//1.32;//1.33;
   double frac=0.05;//fraction of halo in the form of compact objects
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   s.Romins=10000000000.0;
   char filename[40];
   FILE *fill;
   
   
   
   if(numt<10){
   flagf=1; 
   sprintf(filename,"./files/density/%c%.2lf%c%.2lf.dat",'D',s.lat,'_',s.lon);
   fill=fopen(filename,"w");
   if(!fill){cout<<"cannot open file longtitude : "<<s.lon<<"\t latitude: "<<s.lat<<endl;  exit(0);}}


   for(int i=1;i<Num;++i){


   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;
   x=double(i*step);
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = Dsun- x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);
   double rsun= sqrt(zb*zb+ yb*yb+ xb*xb); 


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///Msun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nnf)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*exp(nnf)*exp(-fabs(zb)/0.8)/(1.0+0.5*nnf);///Msun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=frac*(0.932*0.00001/867.067)*pow(0.5/Dsun,-2.44);
   else            s.rho_halo[i]=frac*(0.932*0.00001/867.067)*pow(rdi/Dsun,-2.44);///Msun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///Msun/pc^3

///==================================================================


s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]);///[Msun/pc^3]
s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Msun/deg^2]
s.Nstari[i]= binary_fraction*(s.rho_disk[i]*fd/0.403445+s.rho_ThD[i]*fh/0.4542+s.rho_halo[i]*fh/0.4542+s.rho_bulge[i]*fb/0.308571);////[Nt/pc^3] 
s.Nstari[i]=s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]
s.Nstart  +=  s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
if(s.Rostari[i]<s.Romins) s.Romins=s.Rostari[i];///source selection

  if(flagf>0)
  fprintf(fill,"%e %e %e %e  %e  %e %e\n",x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.Rostar0[i],s.Nstari[i]);}
  if(flagf>0) fclose(fill);    
   cout<<"Nstart [Nt/deg^2]: "<<s.Nstart<<"\t Ro_star [Mass/deg^2]: "<<s.Rostart<<endl;
   cout<<">>>>>>>>>>>>>>>>>>>>>>>>> END OF DISK MODLE <<<<<<<<<<<<<<<<<<<<"<<endl;    
}
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
///===========================================================================//
double RandN(double sigma, double nnd){
   double rr,f,frand;
   do{
   rr=RandR(-sigma*nnd, sigma*nnd);//[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/sigma/sigma);
   frand=RandR(0.0,1.0);
   }while(frand>f);
   return rr;
}
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex)
{
  double F=-1.0;
  if(ds<ex.dis[0])       F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{ 
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] && ds<ex.dis[i+1]){
  F = ex.Extks[i]+(double)(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0||F<0.0){cout<<"ERROR big Ext(ds): "<<F<<"\t ds: "<<ds<<endl; int uut; cin>>uut;}
  return(F);
}
///#############################################################################
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
   double sig,Lon,Lat;
   int uue, flag=0;
   if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
   else sig=1.0;
   double delt=fabs(s.lon)-floor(fabs(s.lon));

     
   if(delt>1.0 or delt<0.0){cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
   else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
   else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
   else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
   else               Lon=(floor(fabs(s.lon))+0.75)*sig;
   if(fabs(s.lon)<0.24999999)      Lon=360.00;
   if(fabs(s.lon-360.0)<0.2499999) Lon=360.00;
   //cout<<"s.lat:  "<<s.lon<<"\t s.lat:  "<<s.lat<<endl;
   //cout<<"Lon:    "<<Lon<<"\t     Lat:  "<<Lat<<endl;


   if(s.lat<0.0) sig=-1.0;
   else sig=1.0;
   delt=fabs(s.lat)-floor(fabs(s.lat));
   if(delt>1.0 or delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
   else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
   else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
   else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
   else                Lat=(floor(fabs(s.lat))+0.75)*sig;
    
   if(fabs(Lon)<0.2499999)  Lon=360.00; 
   if(Lat==-0.00)  Lat=0.00;
   if(fabs(Lon)<0.24999999) Lon=360.00;
   cout<<"Lon:    "<<Lon<<"\t     Lat:  "<<Lat<<endl;
     
     
     
     
   char filename[40];
   FILE *fpd;
   sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',float(Lat),'_',float(Lon) );
   fpd=fopen(filename,"r");
   if(!fpd){
   cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
   //FILE *SD;
   //SD=fopen("./files/Ext/saved_direction.txt","r");
   //for(int i=0; i<64881; ++i) {
   //fscanf(SD,"%lf %lf \n",&latit,&lonti);
   //if(fabs(Lat-latit)<0.1 and fabs(Lon-lonti)<0.1){
   //cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
   //cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
   flag=-1;}
   else{
   flag=1;
   for(int i=0; i<100; ++i){
   fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
   if(ex.dis[i]<0.2  or ex.dis[i]>50.0 or ex.Extks[i]<0.0){
   cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
   cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0;}}
   fclose(fpd);}
   //cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<endl
   return(flag);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
  if (l.Dl==0.0)  l.Dl=0.00034735;
  double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*l.Dl*cos(s.TET)*cos(s.FI));
  double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*s.Ds*cos(s.TET)*cos(s.FI));
  if(Rlc==0.0) Rlc=0.0000000000034346123;
  if(Rsc==0.0) Rsc=0.000000000004762654134; 
 
  double LVx, SVx;
  double SVT, SVR, SVZ, LVT, LVR, LVZ;
  double fv, testfv, test, age;
  double  VSunx, vls2, vls1;
  double betal, betas, deltal, deltas, tetd ;



  double NN=2.5;
  double sigma_R_Disk, sigma_T_Disk, sigma_Z_Disk;
  double sigma_R_DiskL,  sigma_T_DiskL, sigma_Z_DiskL;
  double sigma_R_DiskS,  sigma_T_DiskS, sigma_Z_DiskS;
  double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
  double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
  double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
  double Rho[8]={00.0}; double maxr=0.0;
  for(int i=0; i<8; ++i){ Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}

 
for (int i=0;i<2; ++i){
     test= RandR(0.0,maxr);
     if(test<=Rho[0])                                     {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0; age= 0.075;}
else if(test<=(Rho[0]+Rho[1]))                            {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0; age=0.575; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]))                     {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;age=1.5;  }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]))              {sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2; age=2.5; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))       {sigma_R_Disk=36.7;sigma_T_Disk=23.7; sigma_Z_Disk=15.8; age=4.0; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5])){sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4; age=6.0; }
else if(test<=maxr)                                       {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5; age=8.5; }
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}


    if(i==0) {
    sigma_R_DiskS= sigma_R_Disk;
    sigma_T_DiskS= sigma_T_Disk;
    sigma_Z_DiskS= sigma_Z_Disk;}

    if(i==1){
    sigma_R_DiskL= sigma_R_Disk;
    sigma_T_DiskL= sigma_T_Disk;
    sigma_Z_DiskL= sigma_Z_Disk; }}
    
    
    

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_DiskS, NN);
    SVT= RandN(sigma_T_DiskS, NN);
    SVZ= RandN(sigma_Z_DiskS, NN); }

    else if(s.struc==1){///Galactic bulge
    SVR= RandN(sigma_R_Bulge, NN);
    SVT= RandN(sigma_T_Bulge, NN);
    SVZ= RandN(sigma_Z_Bulge, NN); }

    else if(s.struc==2){///thick disk
    SVR= RandN(sigma_R_TDisk, NN);
    SVT= RandN(sigma_T_TDisk, NN);
    SVZ= RandN(sigma_Z_TDisk, NN); }

    else if(s.struc==3){///stellar halo
    SVR= RandN(sigma_R_halo, NN);
    SVT= RandN(sigma_T_halo, NN);
    SVZ= RandN(sigma_Z_halo, NN); }
    
   
        
    if(s.struc==0 or s.struc==2)  SVT =SVT+ vro_sun*(1.00762*pow(Rsc/Dsun,0.0394) + 0.00712);
    
    s.vs=sqrt( SVR*SVR + SVT*SVT + SVZ*SVZ );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   if(l.struc==0){///Galactic disk
   LVR= RandN(sigma_R_DiskL, NN) ;
   LVT= RandN(sigma_T_DiskL, NN) ;
   LVZ= RandN(sigma_Z_DiskL, NN) ; }

   else if(l.struc==1){///Galactic bulge
   LVR= RandN(sigma_R_Bulge, NN) ;
   LVT= RandN(sigma_T_Bulge, NN) ;
   LVZ= RandN(sigma_Z_Bulge, NN) ; }

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN) ;
   LVT= RandN(sigma_T_TDisk, NN) ;
   LVZ= RandN(sigma_Z_TDisk, NN) ; }

   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN);
   LVT= RandN(sigma_T_halo, NN);
   LVZ= RandN(sigma_Z_halo, NN); }
   
  
   
   
   if(l.struc==0 or l.struc==2)  LVT = LVT+vro_sun *(1.00762*pow(Rlc/Dsun,0.0394) + 0.00712);
   l.vl=sqrt( LVT*LVT + LVZ*LVZ + LVR*LVR );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  BETA  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    betal=betas=0.0;
    tetd= s.TET;
    test= double(l.Dl*cos(s.FI)*sin(tetd)/Rlc);
    if(fabs(test-1.0)<0.01)       betal= pi/2.0;
    else if(fabs(test+1.0)<0.01)  betal=-pi/2.0;
    else                          betal=asin(test);
    
    test= double(s.Ds*cos(s.FI)*sin(tetd)/Rsc); 
    if( fabs(test-1.0)<0.01)     betas=pi/2.0;
    else if(fabs(test+1.0)<0.01) betas=-pi/2.0;
    else                         betas=asin(test);
    
    if(Dsun < fabs(l.Dl*cos(s.FI)*cos(tetd))) betal= pi-betal; 
    if(Dsun < fabs(s.Ds*cos(s.FI)*cos(tetd))) betas= pi-betas; 

   if(fabs(l.Dl*cos(s.FI)*sin(tetd))>Rlc or fabs(test)>1.0){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<tetd<<"\t betal:  "<<betal<<endl;
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<"\t betas:   "<<betas<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(tetd)/Rlc<<"\t sin(s): "<<test<<endl;
   int ew; cin>>ew;}
///HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.TET>pi)  tetd=s.TET-2.0*pi; 
    deltal= pi - fabs(tetd) -fabs(betal);
    deltas= pi - fabs(tetd) -fabs(betas);  
    if(betal<0.0)  deltal= -1.0*deltal;
    if(betas<0.0)  deltas= -1.0*deltas;
    s.deltao= pi-fabs(tetd);
    if(tetd<0.0)  s.deltao=-1.0*s.deltao; 

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    s.SV_n1 =+SVR * sin(deltas)- SVT * cos(deltas);
    s.LV_n1 =+LVR * sin(deltal)- LVT * cos(deltal);
    s.VSun_n1=+VSunR*sin(s.deltao)-VSunT*cos(s.deltao);
    
    SVx= -SVR*cos(deltas)- SVT*sin(deltas);
    LVx= -LVR*cos(deltal)- LVT*sin(deltal);
    VSunx= -VSunR*cos(s.deltao) -VSunT*sin(s.deltao);
    
    s.SV_n2=-sin(s.FI)*(SVx) + cos(s.FI)*SVZ;
    s.LV_n2=-sin(s.FI)*(LVx) + cos(s.FI)*LVZ;
    s.VSun_n2=-sin(s.FI)*(VSunx)+cos(s.FI)*(VSunZ);
 
    
    vls1= l.xls*s.SV_n1 - s.LV_n1 +(1.0-l.xls)*s.VSun_n1;  ///Source - lens 
    vls2= l.xls*s.SV_n2 - s.LV_n2 +(1.0-l.xls)*s.VSun_n2;  /// Source -lens
    l.Vt=sqrt(fabs( vls1*vls1 + vls2*vls2 ) );
    

    //if(vls1>0.0 and vls2>=0.0)         s.xi=atan(fabs(vls2)/fabs(vls1));//in radian
    //else if(vls1<=0.0 and vls2>0.0)    s.xi=atan(fabs(vls1)/fabs(vls2)) + M_PI/2.0;
    //else if(vls1<0.0 and vls2<=0.0)    s.xi=atan(fabs(vls2)/fabs(vls1)) + M_PI;
    //else if(vls1>=0.0 and vls2<0.0)    s.xi=atan(fabs(vls1)/fabs(vls2)) + 3.0*M_PI/2.0;
    //else if(s.xi > 2.0*M_PI)           s.xi-= 2.0*M_PI; 
    
    
    
    if(vls1==0.0 and vls2==0.0){
    cout<<"Big error both zeros:  "<<vls1<<"\t vls_b:  "<<vls2<<endl;
    int iee;  cin>>iee;}
    
    
    if (l.Vt<0.0 or l.Vt>1.0e6){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<s.vs<<endl;
    cout<<"source type:  "<<s.cl1<<"\t s.struc:  "<<s.struc<<endl;
    int yee; cin>>yee;}
    //cout<<"(vrel_func):   s.deltao:  "<<s.deltao<<"FI: "<<s.FI<<endl;   
///*****************************************************
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
