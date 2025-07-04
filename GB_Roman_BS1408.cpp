#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include "VBBinaryLensingLibrary.h"
using namespace std;

const int    Num =5000;
const double MaxD=20.0;///kpc
const double step=MaxD/(double)Num/1.0;///step in kpc
const double RA  =180.0/M_PI;
const double pi  =M_PI; 
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
const double tE_max= 60.0;  
const double Avks  =double(8.20922);

const double Dsun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419,6.2/3.52112,4.0/2.27237,5.8/3.29525,4.9/2.78402, 6.6/3.74991, 3.96/2.24994};

const int M=5;       ///number of filters  VIKH,W149
const double satu[M]={12.0, 12.0, 13.0, 13.0, 14.8}; //it should be changed
const double dete[M]={20.0, 21.0, 21.0, 21.0, 26.0};
const double FWHM[M]={0.33, 0.33, 0.33, 0.33 , 0.33};//3*pixel_size (0.11") of WFIRST, VIKH W149 Z087
const double AlAv[M]={1.009,0.600,0.118,0.184,0.225};///From besancon model[VIKH W149]
const double sigma[M]={0.022,0.022,0.02,0.025,0.025};//MOAاستفاده از مقاله کاردلی
const double cadence=double(15.16/60.0/24.0);//W149_cadence in  day
const double Tobs=double(62.0); 

const int Nw=int(123);///sigma_WFIRST.dat
const int YZ=3615;///number of rows in file yzma
const int N1=36224, N2=25000, N3=3818, N4=3500;///CMD_BESANCON, thinD, bulge, thickD, halo
const int NB=1000;  //Bessel calculation  
const double thre=0.001;//threshold for bessell  

///=============================================================================

struct source{
    int    nums,struc,cl1, cl2;
    double Ds, TET, FI, lat, lon;
    double od_disk, od_ThD, od_bulge, od_halo, opt;///optical  depth
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart, Rostart, Romaxs, Romins, nstart, nstarti;
    double blend[M], Fluxb[M], magb[M], Ai[M], nsbl[M]; 
    double vs, semi, q, period, ecen, inc, tet;
    double radi1,radi2, mass1, mass2, tef1, tef2, typ1, typ2, age, logl1;
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
   double phi, ksi; 
   double tmin, tmax,t0, tp, dto;
   double x1, y1, x2, y2;  
   double uxs1, uys1, uxs2, uys2; 
   double Astar1, Astar2, Astar0; 
   double ulx, uly, f21;   
   double chi1, chi2, chib, chit, chif, dchiL;
};
struct CMD{  
    double tef_d[N1],age_d[N1],logl_d[N1],Mab_d[5][N1],Rs_d[N1],mass_d[N1],type_d[N1]; int cl_d[N1];  ///thin disk
    double tef_b[N2],age_b[N2],logl_b[N2],Mab_b[5][N2],Rs_b[N2],mass_b[N2],type_b[N2]; int cl_b[N2];  /// bulge
    double tef_t[N3],age_t[N3],logl_t[N3],Mab_t[5][N3],Rs_t[N3],mass_t[N3],type_t[N3]; int cl_t[N3];  ///thick disk
    double tef_h[N4],age_h[N4],logl_h[N4],Mab_h[5][N4],Rs_h[N4],mass_h[N4],type_h[N4]; int cl_h[N4];  /// halo
};
struct extinc{
    double dis[100];///distance
    double Extks[100];///ks-band extinction
    double Aks;
    int    flag;
};
struct cmdw{
    double MJ, MH; 
    double Age1[YZ]; double B1[YZ];  double M1[YZ];   double mm1[YZ];
    double Age2[YZ]; double B2[YZ];  double M2[YZ];   double mm2[YZ];
    double Metal[70];  
    int    number[70]; int count[70];   
};
struct roman{
     double magw[Nw], errw[Nw];
};
///===================== FUNCTION ==============================================
int    Extinction(extinc & ex,source & s);
void   read_cmd(CMD & cm , cmdw & w);
void   func_source(source & s, CMD & cm, extinc & ex);
void   func_lens(lens & l, source & s, detection & d);
void   RomanFilter(cmdw & w,double,double,double,double);
void   vrel(source & s, lens & l);
void   Disk_model(source & s, int);
void   optical_depth(source & s);
void   Secondary(source & s, CMD & cm);
double errwfirst(roman & ro, double ghadr, int flag);  
double Interpol(double ds, extinc & ex);
double ThirdKepler(source & s);  
double Kepler(double , double); 
double Bessel(int n,double x); 
double RandN(double , double);
double RandR(double , double);
///============================================================================
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
   lens l;
   source s;
   detection d;
   extinc ex;
   roman ro; 
   cmdw  w;   
   
   char filenam1[40], filenam2[40];  ;
   FILE* fil1; 
   FILE* fil2; 
   fil1=fopen("./files/light/lcF/files/resGBRomanA.dat", "a+");
   fil2=fopen("./files/light/lcF/files/resGBRomanB.dat", "a+");   
   fclose(fil1);  
   fclose(fil2);   
   FILE *model;
   FILE *data1;  

     
     
     
   FILE * roman;
   roman=fopen("./files/sigma_WFIRST.txt","r");
   if(!roman){cout<<"cannot read sigma_WFIRST.txt: "<<endl; int yye;  cin>>yye;}
   for(int i=0; i<Nw;++i){
   fscanf(roman,"%lf  %lf\n",&ro.magw[i],&ro.errw[i]);}
   fclose(roman);
   cout<<"*****  sigma_wfirst.txt was read *****  "<<endl;
    
    
    
   FILE *meta;
   meta=fopen("./files/CMD_WFIRST/metal.txt","r");
   if(!meta){cout<<"cannot read metal.txt:    "<<endl; int yye;  cin>>yye;}
   for(int i=0; i<70; ++i){
   fscanf(meta,"%lf   %d    %d\n",&w.Metal[i],&w.count[i],&w.number[i]);}
   fclose(meta);
   cout<<"*****  metal.txt was read *****  "<<endl; 
   
    
   FILE *hks;
   hks=fopen("./files/CMD_WFIRST/HKS.txt", "r");
   if(!hks){cout<<"cannot read HKS.txt:    "<<endl; int yye;  cin>>yye;}
   for(int i=0; i<YZ; ++i){
   fscanf(hks,"%lf  %lf  %lf  %lf\n",&w.Age1[i],&w.mm1[i],&w.B1[i],&w.M1[i]);}
   fclose(hks);
   cout<<"***** HKS.txt was read *****  "<<endl; 
    
    
   FILE *ji;
   ji=fopen("./files/CMD_WFIRST/JI.txt", "r");
   if(!ji) {cout<<"cannot read JI.txt:    "<<endl; int yye;  cin>>yye;}
   for(int i=0; i<YZ; ++i){
   fscanf(ji,"%lf   %lf   %lf  %lf\n",&w.Age2[i], &w.mm2[i],&w.B2[i],&w.M2[i]);}
   fclose(ji);





   read_cmd(cm,w);     
   double obs1[6][2]={0.0};
   obs1[0][0]=0.0*year+0.0;    obs1[0][1]=0.0*year+62.0;
   obs1[1][0]=0.0*year+182.0;  obs1[1][1]=0.0*year+182.0+62.0;
   obs1[2][0]=1.0*year+0.0;    obs1[2][1]=1.0*year+62.0;
   obs1[3][0]=3.0*year+182.0;  obs1[3][1]=3.0*year+182.0+62.0;
   obs1[4][0]=4.0*year+0.0;    obs1[4][1]=4.0*year+62.0;
   obs1[5][0]=4.0*year+182.0;  obs1[5][1]=4.0*year+182.0+62.0;
    
   int    Ndet=0, flag=0, fi, ndata, fl1, fl2, fl3, threef, li, flago;
   double u1, u2, u0, x0, y0, lonn;
   double prob2, magnio, errg, erra, tim, timc, lonn0, lat0;
   double delta, mag1, mag2, magt, magb, magf, Atot; 
   double nfrac1=0.0 , nfrac2=0.0;  
   double Atoto, A1, A2, Af;  
///===================== Monte Carlo Simulation ================================
   for(int nsim=1450;  nsim<20000; ++nsim){
   
   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   
   
   
   li=int(RandR(1.01,7.9));
   if(li==1){lonn0=1.3;  lat0=-0.875;}
   if(li==2){lonn0=0.9;  lat0=-0.875;}
   if(li==3){lonn0=1.3;  lat0=-1.625;}
   if(li==4){lonn0=0.9;  lat0=-1.675;}
   if(li==5){lonn0=0.5;  lat0=-1.675;}
   if(li==6){lonn0=0.1;  lat0=-1.675;}
   if(li==7){lonn0=-0.3; lat0=-1.675;}
   s.lat= lat0+ double(RandR(-0.375,0.375));
   lonn =lonn0+ double(RandR(-0.200,0.200));
   if(lonn<=0.0) s.lon=360.0+lonn;
   else          s.lon=lonn;
   s.TET=double(360.0-s.lon)/RA;///radian 
   s.FI= double(s.lat/RA); 
   
   Disk_model(s,nsim); 
   ex.flag=-1.0;  
   ex.flag=Extinction(ex,s); 

   do{
   func_source(s,cm,ex);
   func_lens(  l,s,d);
   }while(l.tE<=0.01 or l.tE>tE_max); 
   optical_depth(s);
   

   
   flag=0; 
   if(nsim>=0){
   flag=1; 
   sprintf(filenam1,"./files/light/lcF/files/%c%d.dat",'m',nsim);
   sprintf(filenam2,"./files/light/lcF/files/%c%d.dat",'l',nsim);
   model=fopen(filenam1,"w"); 
   data1=fopen(filenam2,"w");}
    
    
    
   d.peak=0;      timc=0.0;  d.dchiL=0.0;  ndata=0;     threef=0; 
   d.chi2=0.0;  d.chi1=0.0;  d.chib=0.0;   d.chit=0.0;  d.chif=0.0;   
   fl1=0;  fl2=0;   fl3=0; 
   for(tim=d.tmin;  tim<d.tmax; tim+=d.dto){
   
   d.phi=double(tim-d.tp)*2.0*pi/s.period;   
   if(s.ecen<0.01) d.ksi=d.phi;//circular
   else            d.ksi=Kepler(d.phi,s.ecen);//elliptic
   x0=s.semi*(cos(d.ksi)-s.ecen);//_//[m]
   y0=s.semi*sin(d.ksi)*sqrt(1.0-s.ecen*s.ecen);//[m]
   d.x1= -sin(s.inc)*(-y0*sin(s.tet) + x0*cos(s.tet) )*s.mass2/(s.mass1+s.mass2);//[m]  
   d.x2= +sin(s.inc)*(-y0*sin(s.tet) + x0*cos(s.tet) )*s.mass1/(s.mass1+s.mass2);//[m]
   d.y1= +(y0*cos(s.tet)+x0*sin(s.tet))*s.mass2/(s.mass1+s.mass2);//[m] 
   d.y2= -(y0*cos(s.tet)+x0*sin(s.tet))*s.mass1/(s.mass1+s.mass2);//[m] 
   d.ulx=-l.u0*sin(l.ksi)+(tim-d.t0)/l.tE*cos(l.ksi);  
   d.uly=+l.u0*cos(l.ksi)+(tim-d.t0)/l.tE*sin(l.ksi); 

   d.x1=  d.x1*l.xls/l.RE; 
   d.y1=  d.y1*l.xls/l.RE; 
   d.x2=  d.x2*l.xls/l.RE; 
   d.y2=  d.y2*l.xls/l.RE;
   d.uxs1=d.ulx - d.x1;  
   d.uys1=d.uly - d.y1;  
   d.uxs2=d.ulx - d.x2;  
   d.uys2=d.uly - d.y2;
     
   u0=sqrt(d.ulx*d.ulx   + d.uly*d.uly  );//Lens-CM
   u1=sqrt(d.uxs1*d.uxs1 + d.uys1*d.uys1);//Lens-Source1 
   u2=sqrt(d.uxs2*d.uxs2 + d.uys2*d.uys2);//Lens-Source2 
   d.Astar0=vbb.ESPLMag2(u0, 0.5*(l.ro_star1+l.ro_star2) );  
   d.Astar1=vbb.ESPLMag2(u1, l.ro_star1);  
   d.Astar2=vbb.ESPLMag2(u2, l.ro_star2);
 
   mag1= -2.5*log10(pow(10.0,-0.4*s.Map1[4])*(d.Astar1-1.0)+ s.Fluxb[4] );//W149
   mag2= -2.5*log10(pow(10.0,-0.4*s.Map2[4])*(d.Astar2-1.0)+ s.Fluxb[4] );//W149
   magt= -2.5*log10(pow(10.0,-0.4*s.Map1[4])*(d.Astar1-1.0)+ pow(10.0,-0.4*s.Map2[4])*(d.Astar2-1.0)+s.Fluxb[4]);//W149
   magf= -2.5*log10(pow(10.0,-0.4*s.Map1[4])*(d.Astar0-1.0)+ pow(10.0,-0.4*s.Map2[4])*(d.Astar0-1.0)+s.Fluxb[4]);//W149
   
   A1=  double(pow(10.0,-0.4*s.Map1[4])*(d.Astar1-1.0)+ s.Fluxb[4])/s.Fluxb[4]; 
   A2=  double(pow(10.0,-0.4*s.Map2[4])*(d.Astar2-1.0)+ s.Fluxb[4])/s.Fluxb[4];
   Atot=double(pow(10.0,-0.4*magt)/s.Fluxb[4]);
   Af=  double(pow(10.0,-0.4*magf)/s.Fluxb[4]);
   
   if(flag>0) 
   fprintf(model,"%.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf\n", 
   tim,d.ulx,d.uly,d.uxs1, d.uys1, d.uxs2, d.uys2, A1, A2, Atot);//10
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(tim>0.0 and tim<Tobs){  
    flago=0; 
    for(int kl=0; kl<6; kl+=1){
    if((obs1[kl][0]-tim)*(obs1[kl][1]-tim)<=0.0){flago=1; break;}}
    if(flago>0){
    timc+=d.dto; 
    if(timc>cadence){
    timc=timc-cadence; 
    prob2=RandR(0.0,100.0);///probability of uniform observation
    if(fabs(tim-d.t0)<double(5.0*cadence)) d.peak=1; 
    if(prob2>=5.0){
    if(magt>=satu[4] and magt<=dete[4]){
    
    errg=errwfirst(ro, magt,0);//photometric error in magnitude
    erra= Atot*fabs(pow(10.0,-0.4*errg)-1.0);///filter
    delta= RandN(errg,2.5);
    magnio=magt+delta; 
    Atoto=Atot+RandN(erra,2.5);
    
    d.chib +=(Atoto- 1.0)*(Atoto- 1.0)/(erra*erra);//baseline
    d.chi1 +=(Atoto-  A1)*(Atoto-  A1)/(erra*erra);//magnification of the first star
    d.chi2 +=(Atoto-  A2)*(Atoto-  A2)/(erra*erra);//magnification of the second star
    d.chit +=(Atoto-Atot)*(Atoto-Atot)/(erra*erra);//lensing from binary source star(real_model)  
    d.chif +=(Atoto-  Af)*(Atoto-  Af)/(erra*erra);//simple light curve due to CM of binary stars  
    ndata+=1;
    fl1=fl2;
    fl2=fl3;
    if(double(Atot+RandN(erra,2.5)-1.0)>double(3.0*erra)) fl3=1;
    else fl3=0; 
    if(ndata>=3 and double(fl1+fl2+fl3*1.0)>2.0) threef=1;  
    if(flag>0) fprintf(data1,"%.4lf   %.4lf   %.5lf\n",tim, Atoto,erra);  
    }}}}}}
    if(flag>0){fclose(data1); fclose(model);} 
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
  
   d.dchiL=fabs(d.chib-d.chit);
   if(d.peak>0 and d.dchiL>400 and threef>0){
   Ndet+=1;  
   if( fabs(s.Map2[4]-s.Map1[4])<2.0) nfrac1+=1.0;  
   if(double(s.semi*sin(s.inc)/AU/s.Ds)<float(0.1*1000.0))  nfrac2+=1.0;}
   //if(fabs(s.Map2[4]-s.Map1[4])<2.0 and double(s.semi*sin(s.inc)/AU/s.Ds)<float(0.1*1000.0)){
   fil1=fopen("./files/light/lcF/files/resGBRomanA.dat", "a+");
   fil2=fopen("./files/light/lcF/files/resGBRomanB.dat", "a+");   
   fprintf(fil1,
   "%d  %d  %.3lf  %.3lf "//4
   "%d  %.4lf  %.4lf  %.4lf   %.4lf  %.4lf  %.3lf  %.3lf  %.4lf  %.4lf  %d  %d  %.2lf  %.2lf  "//18
   "%.8lf  %.5lf  %.5lf  %.5lf  %.5lf   %.5lf  %.5lf  " //25
   "%d  %.4lf  %.4lf  %.5lf  %.5lf  %.5lf  %.6lf  %.9lf  %.9lf   %.6lf  %.5lf  "//36
   "%.5lf  %.5lf  %.5lf  %.5lf  %.5lf   %d  %.2lf  %e  %.1lf  %.1lf  %.1lf %.1lf  %.1lf  %.6lf %.5lf  %d\n", //52
    nsim, Ndet, s.lon, s.lat,//4 
    s.struc, s.Ds, s.vs, s.mass1, s.mass2, s.age, s.tef1, s.tef2, s.radi1, s.radi2, s.cl1, s.cl2, s.typ1, s.typ2,//18
    s.semi/Rsun, s.period, s.ecen,  s.q, s.inc*RA, s.tet*RA, s.opt*1.0e6,//25
    l.struc, l.Dl, l.vl, l.Ml, l.Vt, l.RE/AU, l.tE, l.ro_star1, l.ro_star2, l.u0, l.ksi*RA,//36
    d.t0, d.tmin, d.tmax, d.tp, d.dto, d.peak, d.dchiL,d.f21,d.chib, d.chi1,d.chi2,d.chit,d.chif,s.blend[4],s.magb[4],threef);//52
    for(int i=0; i<M; ++i){
    fprintf(fil2,"%d  %d  %.4lf  %.1lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.3lf\n",
    nsim,i,s.magb[i],s.nsbl[i],s.blend[i],s.Mab1[i],s.Map1[i],s.Mab2[i],s.Map2[i],s.Ai[i]);}
    fclose(fil1);  
    fclose(fil2); 
    cout<<"============================================================="<<endl;
    cout<<"nsim:  "<<nsim<<"\t Ndet:  "<<Ndet<<"\t ndata: "<<ndata<<endl;
    cout<<"Galactic longitude:  "<<s.lon<<"\t Galactic Latitude:  "<<s.lat<<endl;
    cout<<"mass1:  "<<s.mass1<<"\t mass2:  "<<s.mass2<<"\t age:  "<<s.age<<endl;
    cout<<"cl1:  "<<s.cl1<<"\t cl2:  "<<s.cl2<<"\t typ1:  "<<s.typ1<<"\t tp2:   "<<s.typ2<<endl;
    cout<<"tef1:  "<<s.tef1<<"\t tef2:  "<<s.tef2<<"\t radi1:  "<<s.radi1<<"\t radi2:  "<<s.radi2<<endl;
    cout<<"semi[Rsun]:  "<<s.semi/Rsun<<"\t period(days): "<<s.period<<"\t ecen:  "<<s.ecen<<endl;
    cout<<"inclination:  "<<s.inc*RA<<"\t tet:  "<<s.tet*RA<<"\t Ds:  "<<s.Ds<<endl;
    cout<<"DL:   "<<l.Dl<<"\t tE:  "<<l.tE<<"\t u0:  "<<l.u0<<endl;
    cout<<"f21:  "<<d.f21<<"\t peak:  "<<d.peak<<"\t Dchi:  "<<d.dchiL<<endl; 
    cout<<"threef: "<<threef<<endl;
    cout<<"nfrac1:  "<<nfrac1<<"\t nfrac2:  "<<nfrac2<<endl;
    cout<<"============================================================="<<endl; }
    fclose(_randStream);
    return(0);
}
////==========================================================================
double RandR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return( p*(up-down)+down );
}
///===========================================================================
double errwfirst(roman & ro, double ghadr, int flag){
     double error=-10.0, shib=0.0; 
     
     if(flag==0){///calculating photometry errors ROMAN
     if(ghadr<ro.magw[0] or  ghadr==ro.magw[0])       error= ro.errw[0];
     
     else if(ghadr>ro.magw[Nw-1] or ghadr==ro.magw[Nw-1]){
     shib=(ro.errw[Nw-1]-ro.errw[Nw-2])/(ro.magw[Nw-1]-ro.magw[Nw-2]); 
     error=ro.errw[Nw-1]+shib*(ghadr-ro.magw[Nw-1]);     }
     
     else{
     for(int i=1; i<Nw; ++i){
     if(double((ghadr-ro.magw[i])*(ghadr-ro.magw[i-1]))<0.0 or  ghadr==ro.magw[i-1]){
     shib=(ro.errw[i]-ro.errw[i-1])/(ro.magw[i]-ro.magw[i-1]); 
     error=ro.errw[i-1]+shib*(ghadr-ro.magw[i-1]); 
     break;}}}}
     return(error);     
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
void func_source(source & s, CMD & cm, extinc & ex)
{
    int    struc, nums, num, yye;
    double rho,rf, cosi;
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
    nums=int(RandR(5.0,Num-5.0));
    rho=RandR(s.Romins , s.Romaxs); 
    Ds=double(nums*step);
    }while(rho>s.Rostari[nums] or Ds<0.0 or Ds>MaxD);
    //if(Ds>MaxD or Ds<0.0){
    //cout<<"ERROR (1): Ds: "<<Ds<<"\t MaxD: "<<MaxD<<"\t step: "<<step<<"\t num: "<<nums<<endl; int yye; cin>>yye;}
    //cout<<"Ds:  "<<Ds<<endl;
    
    
    rf=RandR(0.0, s.Rostar0[nums]); 
     if (rf<= s.rho_disk[nums])                    struc=0;///thin disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums])) struc=3;///halo




    if(struc==0){//thin disk
    num=int(RandR(0.0,N1-2.0) );
    for(int i=0; i<5; ++i){Mab[i]=cm.Mab_d[i][num]; }
    if(k==1){
    s.age=  cm.age_d[num];
    s.radi1=cm.Rs_d[num];
    s.typ1= cm.type_d[num];
    s.mass1=cm.mass_d[num];
    s.tef1= cm.tef_d[num];
    s.logl1=cm.logl_d[num];
    s.cl1=  cm.cl_d[num]; }
    if(cm.mass_d[num]<0.0 or int(cm.cl_d[num])>5 or cm.tef_d[num]<0.0 or cm.type_d[num]>8.0 or cm.mass_d[num]>1000.0){
    cout<<"Error(thin disk) temp: "<<cm.tef_d[num]<<"\t mass: "<<cm.mass_d[num]<<"\t counter: "<<num<<endl; cin>>yye; }}


    if(struc==1){//bulge
    num=int(RandR(0.0,N2-2.0));
    for(int i=0; i<5; ++i){ Mab[i]=cm.Mab_b[i][num]; }
    if(k==1){
    s.age=  cm.age_b[num];
    s.radi1=cm.Rs_b[num];
    s.typ1= cm.type_b[num];
    s.mass1=cm.mass_b[num];
    s.tef1= cm.tef_b[num];
    s.logl1=cm.logl_b[num];
    s.cl1=  cm.cl_b[num];}
    if(cm.mass_b[num]<0.0 or int(cm.cl_b[num])>5 or cm.tef_b[num]<0.0 or cm.type_b[num]>8.0 or cm.mass_b[num]>10000.0){
    cout<<"Error(bulge) temp: "<<cm.tef_b[num]<<"\t mass: "<<cm.mass_b[num]<<"\t counter: "<<num<<endl;   cin>>yye; }}



    if(struc==2){//thick disk
    num=int(RandR(0.0,N3-2.0));
    for(int i=0; i<5; ++i){ Mab[i]=cm.Mab_t[i][num]; }
    if(k==1){
    s.age=  cm.age_t[num];
    s.radi1=cm.Rs_t[num];
    s.typ1= cm.type_t[num];
    s.mass1=cm.mass_t[num];
    s.tef1= cm.tef_t[num];
    s.logl1=cm.logl_t[num];
    s.cl1=  cm.cl_t[num]; }
    if( cm.mass_t[num]<0.0 or int(cm.cl_t[num])>5 or cm.tef_t[num]<0.0 or cm.type_t[num]>8.0 or cm.mass_t[num]>100000.0){
    cout<<"Error(thick disk) temp: "<<cm.tef_t[num]<<"\t mass: "<<cm.mass_t[num]<<"\t counter: "<<num<<endl;  cin>>yye;}}




    if(struc==3){///stellar halo
    num=int(RandR(0.0,N4-2.0));
    for(int i=0; i<5; ++i){ Mab[i]=cm.Mab_h[i][num]; }
    if(k==1){
    s.age=  cm.age_h[num];
    s.radi1=cm.Rs_h[num];
    s.typ1= cm.type_h[num];
    s.mass1=cm.mass_h[num];
    s.tef1= cm.tef_h[num];
    s.logl1=cm.logl_h[num];
    s.cl1=  cm.cl_h[num];}
    if(cm.mass_h[num]<0.0 or int(cm.cl_h[num])>5 or cm.tef_h[num]<0.0 or cm.type_h[num]>8.0 or cm.mass_h[num]>10000.0){
    cout<<"Error(Galactic halo) temp: "<<cm.tef_h[num]<<"\t mass: "<<cm.mass_h[num]<<"\t counter: "<<num<<endl;  cin>>yye;}}




    Mab[4]=(Mab[2]+Mab[3]+Mab[4])/3.0;//absolute magnitude in W149: (K+H+J)/3    
    if(ex.flag>0) Av=double(Interpol(Ds,ex)*Avks);
    else          Av=double(0.7*Ds); 
    if(Av<0.0)    Av=0.0;
    for(int i=0; i<M; ++i){
    Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i],1.5);
    if(Ai[i]<0.0) Ai[i]=0.0;   
    Map[i]=Mab[i]+5.0*log10(Ds*100.0)+Ai[i];
    if(k==1){
    s.Ai[i]=Ai[i];   
    s.Map1[i]=Map[i];   
    s.Mab1[i]=Mab[i];}
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
    s.period=ThirdKepler(s);//[days]
    emax=double(0.8-8.0*exp(-pow(6.0*s.period,0.35)));//Fig 7
    if(emax<0.0)  emax=0.0; 
    if(emax>1.0)  emax=1.0;
    if(emax<0.0 or emax>1.0){cout<<"Error emax:  "<<emax<<"\t period:  "<<s.period<<endl; cin>>yye; } 
    s.ecen=RandR(0.0, emax);
    
    cosi=RandR(0.0,0.9999999);
    s.inc=acos(cosi);//radian 
    
    //s.inc=RandR(0.1,89.9)/RA; 
    s.tet=RandR(0.1, 359.9)/RA;
    for(int i=0; i<M; ++i){
    s.Map2[i]=s.Mab2[i]+5.0*log10(s.Ds*100.0)+s.Ai[i];
    s.Fluxb[i]+=pow(10.0,-0.4*s.Map2[i]);} 
    k+=1;}
        
    }///loop over the stars
 


    for(int i=0; i<M; ++i){
    s.magb[i]=-2.5*log10(s.Fluxb[i]);
    s.blend[i]=double(pow(10.0,-0.4*s.Map1[i])+ pow(10.0,-0.4*s.Map2[i]))/s.Fluxb[i];
    if(s.Fluxb[i]<=0.0 or int(s.nsbl[i])<0.0 or s.nsbl[i]==0.0 or s.blend[i]>1.00001 or s.blend[i]<=0.0 or 
    (s.nsbl[i]==2.0 and s.blend[i]<1.0) or fabs(s.Map1[0]-s.Mab1[0]-s.Ai[0]-s.Map1[1]+s.Mab1[1]+s.Ai[1])>0.1 or 
    s.Ds<0.0 or s.Ds>MaxD or s.ecen>1.0 or s.ecen<0.0){ 
    cout<<"nsbl[i]: "<<s.nsbl[M]<<"\t Ds:  "<<s.Ds<<endl; 
    cout<<"BIg ERROR nsbl: "<<s.nsbl[i]<<"\t nsbl[i]: "<<s.nsbl[i]<<"\t blend: "<<s.blend[i]<<endl; 
    cin>>yye;}}
    //cout<<"********** End of Func_source *************"<<endl;
}
///==============================================================//
///                                                              //
///                  Secondary source stars                      //
///                                                              //
///==============================================================//
void Secondary(source & s, CMD & cm){
   double dism=10000000.0,dis;
   int count=0; 
   int kj=-1;  
 
 
 
   if(s.struc==0){///thin disk
   if(s.mass2<=cm.mass_d[0]) kj=0; 
   else if(s.mass2>=cm.mass_d[N1-1]) kj=int(N1-1); 
   else{ 
   for(int i=0; i<int(N1-1); ++i){
   if(float((s.mass2-cm.mass_d[i])*(s.mass2-cm.mass_d[i+1]))<=0.0) {kj=i;  i=N1;}}}
   if(kj<0){cout<<"Error(1) kj<0: "<<kj<<"\t mass2:  "<<s.mass2<<endl;  int yye;  cin>>yye; }
   count=kj; 
   s.mass2=cm.mass_d[count];  
   s.tef2=  cm.tef_d[count];
   s.radi2=  cm.Rs_d[count];  
   s.cl2=    cm.cl_d[count];
   s.typ2= cm.type_d[count];
   for(int i=0; i<M; ++i) s.Mab2[i]=cm.Mab_d[i][count];}
   
   
   kj=-1; 
   if(s.struc==1){///Bulge
   if(s.mass2<=cm.mass_b[0]) kj=0; 
   else if(s.mass2>=cm.mass_b[N2-1]) kj=int(N2-1); 
   else{ 
   for(int i=0; i<int(N2-1); ++i){
   if(float((s.mass2-cm.mass_b[i])*(s.mass2-cm.mass_b[i+1]))<=0.0) {kj=i;  i=N2;}}}
   if(kj<0){cout<<"Error(2) kj<0: "<<kj<<"\t mass2:  "<<s.mass2<<endl;  int yye;  cin>>yye; }
   count=kj; 
   s.mass2=   cm.mass_b[count];  
   s.radi2=     cm.Rs_b[count];
   s.tef2=     cm.tef_b[count];  
   s.cl2=       cm.cl_b[count];
   s.typ2=     cm.type_b[count];
   for(int i=0; i<M; ++i)  s.Mab2[i]=cm.Mab_b[i][count];}
   
   
   kj=-1;
   if(s.struc==2){//Thick disk
   if(s.mass2<=cm.mass_t[0]) kj=0; 
   else if(s.mass2>=cm.mass_t[N3-1]) kj=int(N3-1); 
   else{ 
   for(int i=0; i<int(N3-1); ++i){
   if(float((s.mass2-cm.mass_t[i])*(s.mass2-cm.mass_t[i+1]))<=0.0){kj=i;  i=N3;}}}
   if(kj<0){cout<<"Error(3) kj<0: "<<kj<<"\t mass2:  "<<s.mass2<<endl;  int yye;  cin>>yye; }
   count=kj; 
   s.mass2=   cm.mass_t[count];  
   s.radi2=     cm.Rs_t[count];
   s.tef2=     cm.tef_t[count];  
   s.cl2=       cm.cl_t[count];
   s.typ2=     cm.type_t[count];
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
   s.mass2=   cm.mass_h[count];  
   s.radi2=     cm.Rs_h[count];
   s.tef2=     cm.tef_h[count];  
   s.cl2=       cm.cl_h[count];
   s.typ2=     cm.type_h[count];
   for(int i=0; i<M; ++i) s.Mab2[i]=cm.Mab_h[i][count];}
   s.Mab2[4]=double(s.Mab2[2]+s.Mab2[3]+s.Mab2[4])/3.0;//absolute magnitude in W149: (K+H+J)/3
   
   
   //cout<<"mass2:  "<<s.mass2<<"\t radi2:  "<<s.radi2<<"\t tef2:  "<<s.tef2<<endl;  
   //cout<<"cl2:  "<<s.cl2<<"\t type2:  "<<s.typ2<<endl;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double ThirdKepler(source & s){
    return( sqrt(4.0*pi*pi*s.semi*s.semi*s.semi/(G*Msun*(s.mass1+s.mass2)))/(24.0*3600.0));///days  
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
  l.Ml=RandR(0.08,1.5);
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
  d.tmin=d.t0-1.5*l.tE;
  d.tmax=d.t0+1.5*l.tE;
  d.dto=cadence*1.0;//(d.tmax-d.tmin)/No;
  //if(d.dto>cadence or l.tE>500.0){
  //cout<<"Error cadence"<<cadence<<"\t dto:  "<<d.dto<<"\t tE:  "<<l.tE<<endl;   int uus; cin>>uus;  } 
  d.tp=d.t0+RandR(-s.period*0.5 , s.period*0.5);     
  d.f21=pow(10.0,-0.4*(s.Map2[4]-s.Map1[4]));//W149
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
void RomanFilter(cmdw & w, double metal,  double age, double MK, double MI)
{
    int h, g, k1, k2, uui;       
    h=-1;
    if(metal<w.Metal[0] || metal==w.Metal[0])         h=0;
    else if(metal>w.Metal[69] ||  metal==w.Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if(float( (metal-w.Metal[i-1])*(metal-w.Metal[i]))<0.0 or metal==w.Metal[i-1]){ h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    
   
    k1=int(w.count[h]);   
    k2=int(w.count[h]+w.number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or w.mm1[k1]!=w.mm1[k2-1] or w.number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<w.count[h]<<"\t number[h]: "<<w.number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<w.Age1[k1]<<"\t Age1[k2-1]: "<<w.Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<w.mm1[k1]<<"\t metal[k2-1]: "<<w.mm1[k2-1]<<endl;  cin>>uui;}
    
    g=-1;  
    if(age<w.Age1[k1]  or age==w.Age1[k1])     g=k1;
    else if(age>w.Age1[k2-1] or age==w.Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    
    if(w.Age1[k-1]>w.Age1[k] or w.mm1[k-1]!=w.mm1[k]){    
    cout<<"Bad error: Age1[k-1]: "<<w.Age1[k-1]<<"\t Age1[k]: "<<w.Age1[k]<<endl;
    cout<<"mm[k-1]:"<<w.mm1[k-1]<<"\t mm[k]"<<w.mm1[k]<<endl;    cin>>uui;}

    if((age-w.Age1[k-1])*(age-w.Age1[k])<0.0 or  age==w.Age1[k-1]) {g=int(k-1);   break;}}}
    
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){
    cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    
     w.MH=double(w.B1[g]+w.M1[g]*MK*1.05263157894737); ///H-band  Mks/Mk=1.05263157894737
     w.MJ=double(w.B2[g]+w.M2[g]*MI);   ///J-band
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm, cmdw & w)
{
    //mass, Teff, Age, logL,  log(g),  Z,  Rs,  MB, MV, MI, MK, Cl, type (13)
    int yye;
    double metal, gravity, MB;
    char filename[40];
    FILE *fp2;


////=================================== THIN DISK ==============================
    int j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c%c.dat",'C','M','D','T','i','W','s');///sorted by the mass
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTiW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.tef_d[j],&cm.age_d[j],&cm.logl_d[j],&gravity,&metal,&cm.Rs_d[j],&MB,
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.cl_d[j],&cm.type_d[j]);
    ///*******************************************************
    RomanFilter(w, metal, cm.age_d[j] , double(cm.Mab_d[2][j]) , double(cm.Mab_d[1][j]) );
    cm.Mab_d[3][j]=w.MH; 
    cm.Mab_d[4][j]=w.MJ;
    
    if(fabs(cm.Mab_d[3][j]-cm.Mab_d[2][j])>1.5 or fabs(cm.Mab_d[4][j]-cm.Mab_d[1][j])>1.5 or cm.mass_d[j]<0.0 
    or cm.mass_d[j]==0.0 or cm.tef_d[j]<0.0 or metal>0.1 or cm.age_d[j]>10 or int(cm.cl_d[j])==6 or float(cm.type_d[j])>8.0 or cm.type_d[j]<1.0){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<endl;
    cout<<"ERROR: Mab_d(y-band): "<<cm.Mab_d[4][j]<<"\t Mab_d(z-band): "<<cm.Mab_d[3][j]<<"\t metal: "<<metal<<endl;
    cin>>yye;}
    
    j++;} 
    fclose(fp2);
    if(j!=N1){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;





////=================================== BULGE ==================================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','b','W','s');///sorted by mass
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDbW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.tef_b[j],&cm.age_b[j],&cm.logl_b[j],&gravity,&metal,&cm.Rs_b[j],&MB,
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.cl_b[j],&cm.type_b[j]);
    RomanFilter(w, metal, cm.age_b[j], double(cm.Mab_b[2][j]) , double(cm.Mab_b[1][j]) );
    cm.Mab_b[3][j]=w.MH; 
    cm.Mab_b[4][j]=w.MJ;
   
    if(fabs(cm.Mab_b[3][j]-cm.Mab_b[2][j])>1.5 or fabs(cm.Mab_b[4][j]-cm.Mab_b[1][j])>1.5  or 
    cm.mass_b[j]<0.0 or cm.mass_b[j]==0.0 or cm.tef_b[j]<0.0 or  cm.age_b[j]>10 or metal>0.9 or cm.cl_b[j]==6 or
    cm.type_b[j]>=8.0 or (cm.cl_b[j]==5 and int(cm.type_b[j])>8.0) or (cm.cl_b[j]==6 and int(cm.type_b[j])!=9)){
    cout<<"ERROR: Mab_b(y-band): "<<cm.Mab_b[4][j]<<"\t Mab_b(z-band): "<<cm.Mab_b[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<endl;  cin>>yye;}
   
    j++;} 
    fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
    cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;






////=================================== THICK DISK =============================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c%c.dat",'C','M','D','T','k','W','s');//sorted by mass
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTkW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.tef_t[j],&cm.age_t[j],&cm.logl_t[j],&gravity,&metal,&cm.Rs_t[j],&MB,
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.cl_t[j],&cm.type_t[j]);
    RomanFilter(w, metal, cm.age_t[j], double(cm.Mab_t[2][j]) , double(cm.Mab_t[1][j]) );
    cm.Mab_t[3][j]=w.MH; 
    cm.Mab_t[4][j]=w.MJ;
    
    if(fabs(cm.Mab_t[3][j]-cm.Mab_t[2][j])>1.5 or fabs(cm.Mab_t[4][j]-cm.Mab_t[1][j])>1.5 or 
    cm.mass_t[j]<0.0||  cm.mass_t[j]==0.0 or cm.tef_t[j]<0.0 or metal>0.2||cm.cl_t[j]==6|| cm.type_t[j]>=8.0 or
    (cm.cl_t[j]==5 and int(cm.type_t[j])>8) or (cm.cl_t[j]==6 and int(cm.type_t[j])!=9)){
    cout<<"ERROR: Mab_t(y-band): "<<cm.Mab_t[4][j]<<"\t Mab_t(z-band): "<<cm.Mab_t[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;




////=================================== STELLAR HALO ===========================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','h','W','s');//sorted by mass
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDhW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.tef_h[j],&cm.age_h[j], &cm.logl_h[j],&gravity,&metal,&cm.Rs_h[j],&MB,
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.cl_h[j],&cm.type_h[j]);
    RomanFilter(w, metal, cm.age_h[j] , double(cm.Mab_h[2][j]) , double(cm.Mab_h[1][j])  );
    cm.Mab_h[3][j]=w.MH; 
    cm.Mab_h[4][j]=w.MJ;
    
   
    if(fabs(cm.Mab_h[3][j]-cm.Mab_h[2][j])>1.5  or fabs(cm.Mab_h[4][j]-cm.Mab_h[1][j])>1.5 or 
    cm.mass_h[j]<0.0 or cm.mass_h[j]==0.0 or  cm.age_h[j]<0 or cm.cl_h[j]<0  or cm.cl_h[j]==6  or  cm.tef_h[j]<0.0 or
    metal>0.1 or cm.cl_h[j]>7 or cm.type_h[j]>9 or (cm.cl_h[j]==5 and int(cm.type_h[j])>8) or (cm.cl_h[j]==6 and int(cm.type_h[j])!=9) or (cm.cl_h[j]<5 and int(cm.type_h[j])==9)){
    cout<<"ERROR: Mab_h(y-band): "<<cm.Mab_h[4][j]<<"\t Mab_h(z-band): "<<cm.Mab_h[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
   cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
   cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s, int numt)
{
   int    flagf=0;
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
