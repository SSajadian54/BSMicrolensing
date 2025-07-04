### This file was made on 8/4/1404
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib.figure import Figure
from IPython.display import display, Math
import emcee
from numpy import conj
import pylab
from matplotlib import cm
from matplotlib import colors
cm=colors.ListedColormap(['purple', 'blue', 'darkgreen','yellowgreen', 'orange', 'red'])
from matplotlib import gridspec
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)



cmap=plt.get_cmap('viridis')
mpl.rc('font',**{'family':'serif','serif':['Palatino']})
mpl.rc('text', usetex=True)
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams["font.size"] = 11.5
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Computer Modern Sans"]
mpl.rcParams["text.usetex"] = True
mpl.rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
mpl.rcParams['text.usetex'] = True

import VBMicrolensing
vbb = VBMicrolensing.VBMicrolensing()
vbb.RelTol = 1e-04
vbb.Tol=1e-04
vbb.LoadESPLTable("./ESPL.tbl"); 
vbb.a1=0.0

###############################################################################

labels=['u0', 't0', 'tE', 'fb', 'ros']
nam=[r"$\rm{Nsimulation}$", r"$\rm{Ndetection}$",r"$\rm{Galactic}~\rm{longitude}$", r"$\rm{Galactic}~\rm{Latitude}$",
     r"$\rm{Structure}_{\star}$", r"$D_{\rm{s}}(\rm{kpc})$", r"$v_{\star}(km/s)$" , r"$M_{\star, 1}(M_{\odot})$", 
     r"$M_{\star, 2}(M_{\odot})$",r"$\mathcal{A}$",  r"$T_{\star, 1}(K)$", r"$T_{\star, 2}(K)$", r"$R_{\star, 1}(R_{\odot})$", 
     r"$R_{\star, 2}(R_{odot})$", r"$\rm{CL}_1$", r"$\rm{CL}_2$", 
     r"$\rm{Type}_{\star, 1}$", r"$\rm{Type}_{\star, 2}$", r"$a/R_{\odot}$", r"$T(\rm{days})$", r"$\epsilon$", 
     r"$q$",r"$i(\rm{deg})$",r"$\theta(\rm{deg})$", r"$\tau \times 10^{6}$",  
     r"$\rm{Structure}_{\rm{l}}$", r"$D_{\rm{l}}(\rm{kpc})$", r"$v_{\rm{l}}(km/s)$",  
     r"$M_{\rm{l}}(M_{\odot})$", r"$v_{\rm{rel}}(km/s)$", r"$R_{\rm{E}}(\rm{AU})$", r"$t_{\rm{E}}(\rm{days})$", 
     r"$\rho_{\star, 1}$", r"$\rho_{\star, 2}$", r"$u_{0}$", r"$\xi(\rm{deg})$", r"$t_{0}$",  
     r"$t_{\rm{min}}(\rm{days})$", r"$t_{\rm{max}}(\rm{days})$", r"$t_{\rm{p}}$", r"$\delta_{t}(\rm{days})$",  
     r"$\rm{peak}~\rm{flag}$", r"$\Delta \chi^{2}_{\rm{L}}$", r"$\mathcal{F}$", r"$\chi^{2}_{\rm{base}}$", 
     r"$\chi^{2}_{\star,1}$", r"$\chi^{2}_{\star,2}$",r"$\chi^{2}_{\rm{real}}$",r"$\chi^{2}_{\rm{CM}}$",r"$f_{b}$", r"$n_{\rm{data}}$", r"$Flag_{3}$"]#52


nam2=[r"$Nsim$", r"$\rm{filter}$" , r"$m_{\rm{base}}(\rm{mag})$", r"$N_{\rm{blend}}$", r"$f_{\rm{b}}$", r"$M_{\rm{ab},~1}(\rm{mag})$", r"$m_{\rm{ap}, 1}(\rm{mag})$", r"$M_{\rm{ab}, 2}(\rm{mag})$", r"$m_{\rm{ap}, 2}(\rm{mag})$", r"$\rm{Extinction}(\rm{mag})$"]#10
    

###############################################################################

sm=1.0e-50
Rsun=6.957*pow(10.0,8.0); #solar radius [meter]
AU= 1.4960*pow(10.0,11.0);
Tobs =float(7*30.0)
LsstAAC=float(0.7*1000.0)## mili arcsec LSST 

ndim =int(5)
Nwal =int(60)
nstep=int(10000)
nc   =int(52)
p0=np.zeros((Nwal,ndim))
pb=np.zeros((ndim,3))
b =np.zeros((ndim*3+2))
prt=np.zeros((ndim))

###############################################################################
def likelihood(p, tim,Magni,erra, prt):
    lp=prior(p,prt)
    if(lp<0.0):
        return(-np.inf)
    return lnlike2(p,tim,Magni,erra)
###############################################################################
def prior(p, prt): 
    u01, t01, tE1, fb1, rho1 =p[0],   p[1],   p[2],   p[3],   p[4] 
    u00, t00, tE0, fb0, rho0 =prt[0], prt[1], prt[2], prt[3], prt[4] 
    if(u01>sm and u01<1.0 and t01>0.0 and t01<Tobs and tE1>0.0 and tE1<250.0 and fb1>sm and fb1<1.0 and rho1>sm and rho1<10.0
    and abs(t01-t00)<tE0 and abs(tE1-tE0)<tE0 and abs(fb1-fb0)<0.8 and abs(rho1-rho0)<float(7.0*rho0)): 
        return(0); 
    return(-1.0);     
###############################################################################
def lnlike2(p, tim,Magni,erra):
    u01, t01, tE1, fb1, rho1=p[0], p[1], p[2], p[3],  p[4]
    vbb.a1=0.0;
    As=np.zeros(len(tim))
    for d in range(len(tim)):  
        u=np.sqrt((tim[d]-t01)*(tim[d]-t01)/tE1/tE1 + u01*u01); 
        if(u>float(20.0*rho1)):  As[d]=(u*u+2.0)/np.sqrt(u*u*(u*u+4.0))
        else:                    As[d]=vbb.ESPLMag2(u,rho1)
    As=fb1*As+1.0-fb1 
    return(-1.0*np.sum((As-Magni)**2.0/(erra*erra)))
###############################################################################
def mcmcp(prt, tim, Magni, erra, NUM, u0FF, t0FF, tEFF):

    p0[:,0]=np.abs(np.random.normal(u0FF,0.2,     Nwal))#u0
    p0[:,1]=np.abs(np.random.normal(t0FF,0.2*tEFF,Nwal))#t0
    p0[:,2]=np.abs(np.random.normal(tEFF,0.2*tEFF,Nwal))#tE
    p0[:,3]=np.abs(np.random.normal(prt[3],0.4,   Nwal))#fb
    p0[:,4]=np.abs(np.random.normal(prt[4],prt[4]*0.6,Nwal))#ros
      
    sampler=emcee.EnsembleSampler(Nwal,ndim,likelihood,args=(tim,Magni,erra, prt), threads=8)
    sampler.run_mcmc(p0,  nstep, progress=True)
    print("END OF MCMC ***********************************")
    ####################################################
    Chain=sampler.get_chain(flat=True)
    for w in range(ndim):
        mcmc=np.percentile(Chain[:,w], [16, 50, 84])
        pb[w,0]= mcmc[1]
        pb[w,1:] = np.diff(mcmc)
        print(str(labels[w]), pb[w,0],  pb[w,1], pb[w,2])
    chi2s=sampler.get_log_prob(flat=True)
    tt=int(np.argmin(-1.0*chi2s))
    Chibf=abs(lnlike2(Chain[tt,:],tim, Magni, erra))
    ffd=  abs(chi2s[tt])
    print("Chi2s, argmin,  param_min, chis: ",chi2s, tt, Chain[tt,:], ffd, Chibf)
    if(abs(Chibf-ffd)>0.1):  
        print("BIG ERROR: chi2s:  ",  Chibf, ffd)
        print("Best-fitted parameters:  ",  Chain[tt, :])
        print("Initial parameters:  ", prt)
        input("Enter a number ")
    u0f, t0f, tEf, fbf, rosf=Chain[tt,0], Chain[tt,1], Chain[tt,2], Chain[tt,3], Chain[tt,4]
    du0=0.5*(abs(pb[0,1])+abs(pb[0,2]));     
    dt0=0.5*(abs(pb[1,1])+abs(pb[1,2]));
    dtE=0.5*(abs(pb[2,1])+abs(pb[2,2]));     
    dfb=0.5*(abs(pb[3,1])+abs(pb[3,2]));
    dro=0.5*(abs(pb[4,1])+abs(pb[4,2]));
    b=np.array([NUM,u0f,pb[0,1],pb[0,2],t0f,pb[1,1],pb[1,2],tEf,pb[2,1],pb[2,2],fbf,pb[3,1],pb[3,2],rosf,pb[4,1],pb[4,2],Chibf])#17
    save=open("./files/BestfitL.dat","a")
    np.savetxt(save,b.reshape(-1,ndim*3+2),fmt='%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.7f %.7f %.7f %.1f')
    save.close()
    ##################################################
    timt=np.arange(t0f-3.0*tEf, t0f+3.0*tEf, 0.005*tEf)
    us=  np.sqrt(u0f**2.0+(timt-t0f)**2.0/tEf/tEf)
    mfit=np.zeros((len(timt)))
    for ii in range(len(timt)):  
        if(us[ii]>float(20.0*rosf)): mfit[ii]=float(us[ii]*us[ii]+2.0)/np.sqrt(us[ii]*us[ii]*(us[ii]*us[ii]+4.0))
        else:                        mfit[ii]=vbb.ESPLMag2(us[ii] , rosf)
    mfit=fbf*mfit+1.0-fbf  
    pn=np.array([tEf, dtE, rosf, dro])   
    return(timt, mfit, Chibf, pn) 
###############################################################################

def Asymetry(array, tim): 
    N=len(array)
    N0=np.argmax(array)
    asym=0.0; 
    numb=0.0
    t0=tim[N0]
    ts=0.0
    for i in range(N-N0):
        n1=int(N0+i)
        n2=int(N0-i);  dtt=0.0
        if(n1<N and n2>=0):
            if(i==0 or i==1): dtt=abs(t0-tim[n1])#days
            else:             dtt=abs(tim[n1]-tim[n1-1])
            ts+=dtt
            if(ts>3.0):#cadence_LSST
                ts=ts-3.0
                numb+=1.0  
                asym+=(array[n1]-array[n2])**2.0
    return( np.sqrt(asym/(numb+0.0000174236)/1.0) )        
###############################################################################

def funcini(tim , Atot, tE):
    ##Full Width of Half Maximum #
    N  =int(len(Atot))
    Am =abs(np.max(Atot-1.0)*0.5)
    dis=np.abs(Atot-(Am+1.0))
    n1 =np.argmin(dis[:int(N/2)])
    t1 =float(tim[n1])
    n2 =np.argmin(dis[int(N/2):])
    t2 =float(tim[n2+int(N/2)])
    FWHM=abs(t2-t1)
    print("FWHM: " ,  t1, t2, n1, n2+int(N/2), FWHM)
    #################################
    ####  U0_fitted && t0_fitted  ##
    u0=np.arange(0.00001, 1.0, 0.0001)
    AA=(u0*u0+2.0)/np.sqrt(u0*u0*(u0*u0+4.0))
    n3=np.argmin(np.abs(AA-np.max(Atot)))
    u0f=u0[n3]
    t0f=tim[np.argmax(Atot)]
    nt0=int(np.argmax(Atot))
    print("u0f, t0f, N: ", u0f, t0f, n3, AA[n3], np.max(Atot), nt0 )
    #################################
    ### tE_fitted  ##
    tEr=np.arange(tE-7.0, tE+7.0 , 0.001)
    u=np.sqrt(((t2-t0f)/tEr)**2.0+ u0f*u0f)
    A=np.abs((u*u+2.0)/np.sqrt(u*u*(u*u+4.0))-(Am+1.0))
    n4=np.argmin(A)
    tEf=tEr[n4] 
    print("tEf: ", tEf, n4, A[n4], Am+1.0)   
    return(u0f, t0f, tEf)
###############################################################################

 
save2=open("./files/BestfitL.dat","a+")
save2.close()

f1=open("./files/resGBLsstA.dat","r")
nf= sum(1 for line in f1)  
par=np.zeros((nf,nc))
detec=np.zeros((nf,nc))
par2=np.zeros((nf*6,10))
par=np.loadtxt("./files/resGBLsstA.dat")  
par2=np.loadtxt("./files/resGBLsstB.dat")  


frac1=0.0
frac2=0.0
simple=int(0)
nlens=int(0)
magb=np.zeros((6)); Map1=np.zeros((6)); Map2=np.zeros((6)); 
################################################################################
for i in range(nf):
    ei, ndet, lon, lat, strucs, Ds =int(par[i,0]), par[i,1], par[i,2], par[i,3], par[i,4], par[i,5]  
    vs, Ms1, Ms2, age, tef1,tef2, Rs1= par[i,6], par[i,7], par[i,8], par[i,9], par[i,10],par[i,11],par[i,12]
    Rs2, cl1, cl2, typ1, typ2, Semi, period=par[i,13],par[i,14],par[i,15],par[i,16],par[i,17],par[i,18],par[i,19]
    ecen, q, inc, tet, opt, strucl, Dl=  par[i,20],par[i,21],par[i,22],par[i,23],par[i,24],par[i,25],par[i,26]
    vl, Ml, Vt, RE,tE, ros1, ros2=           par[i,27],par[i,28],par[i,29],par[i,30],par[i,31],par[i,32],par[i,33]
    u0, ksi,t0, tmin, tmax, tp, dto=par[i,34],par[i,35],par[i,36],par[i,37],par[i,38],par[i,39],par[i,40]
    peak, dchiL, f21, chiB, chi1= par[i,41], par[i,42],par[i,43],par[i,44], par[i,45]
    chi2, chiR, chiCM, fb, ndata, threef  =  par[i,46], par[i,47], par[i,48], par[i,49], int(par[i,50]), par[i,51]

    for k in range(6): 
        cc=int(i*6+k)
        nsim, filt, magb[k], nsbl, blend, Mab1, Map1[k], Mab2, Map2[k], Ext=par2[cc,:]  
   
    Dmag =abs(Map1[2] - Map2[2])#r-band
    sproj=float(Semi*Rsun*np.sin(inc*np.pi/180.0)/AU/Ds)##mas
    cri2=sproj/LsstAAC
    prt=np.array([u0, t0, tE, fb, 0.5*(ros1+ros2)])   
    print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH", ei)
    print("Masses:  ", Ms1, Ms2, Dmag, Map1[2], Map2[2], blend,   nsbl, magb[2])
    #input(" Enter a number ")       
    ############################################################################
    if(dchiL>=300.0 and peak>0 and threef>0 and ndata>1 and ei>3172):
        detec[int(nlens),:]=par[i,:]
        nlens+=1
        f2=open('./files/L{0:d}.dat'.format(ei),"r")
        nd= sum(1 for line in f2)  
        f2=open('./files/M{0:d}.dat'.format(ei),"r")
        nm= sum(1 for line in f2)  
        dat=np.zeros((nd,4))
        dat=np.loadtxt('./files/L{0:d}.dat'.format(ei)) 
        mod=np.zeros((nm,10))
        mod=np.loadtxt('./files/M{0:d}.dat'.format(ei)) 
        
        if(Dmag<2.0):   frac1+=1.0
        if(cri2<=1.0):  frac2+=1.0 
        flag=0;  
        chi2m=0.0
        tEf=0.0;  dtE=0.0;  rosf=0.0;  dro=0.0;  
        if(Dmag<2.0 and cri2<=1.0):  
            flag=1
            u0FF, t0FF, tEFF= funcini(mod[:,0],mod[:,9], tE)
            print("paramo, fitted:  ",  tE,  tEFF, t0,  t0FF, u0, u0FF,  fb,  ros1,  ros2)
            timt, mfit, chi2m, pp= mcmcp(prt, dat[:,0], dat[:,1], dat[:,2]+1.0e-10,  ei, u0FF, t0FF, tEFF)
            tEf, dtE, rosf, dro=pp[0], pp[1], pp[2], pp[3]
            if(abs(chi2m-chiR)<100): simple+=1; 
        print("***************************************************************")
        print("Counter:  ",           ei,frac1,frac2,simple, chi2m, chiR, flag)
        print("t0,  u0,  tE:       ", prt , Dl, Ds, RE, ksi, Ml,  Vt)
        print("Orbit_parameters:   ", Semi,  q,  period, inc,  tet,  ecen)
        print("dchi,  peak, f21:   ", dchiL,   peak, f21,  threef, cri2, Dmag)
        
        
        n1=np.argmin(np.sqrt(mod[:,3]*mod[:,3]+mod[:,4]*mod[:,4]))
        n2=np.argmin(np.sqrt(mod[:,5]*mod[:,5]+mod[:,6]*mod[:,6]))
        t01=mod[n1,0]
        t02=mod[n2,0]
        Dt0=abs(t01-t02)
        u01=np.sqrt(mod[n1,3]*mod[n1,3]+mod[n1,4]*mod[n1,4])
        u02=np.sqrt(mod[n2,5]*mod[n2,5]+mod[n2,6]*mod[n2,6])
        asym=Asymetry(mod[:,9], mod[:,0])/np.mean(np.abs(dat[:,2]))#normalized 
        #######################################################################
        
        ymax1,ymin1= np.max(mod[:,7:9]), np.min(mod[:,7:9])
        ymax2,ymin2= np.max(dat[:,1]),   np.min(dat[:,1])
        if(flag>0):  ymax3, ymin3= np.max(mfit), np.min(mfit)
        
        plt.clf()
        plt.cla()
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        plt.plot(mod[:,0],mod[:,7],'r--',label=r"$\rm{Star}~1$", lw=1.5)
        plt.plot(mod[:,0],mod[:,8],'b:', label=r"$\rm{Star}~2$", lw=1.5)
        plt.plot(mod[:,0],mod[:,9],'k-', label=r"$\rm{Binary}~\rm{System}$",lw=1.7)
        if(flag>0): plt.plot(timt,mfit,'g-.', label=r"$\rm{ESPL}~\rm{Fitted}~\rm{model}$",lw=1.85)        
        for jq in range(nd): 
            nvv=abs(int(dat[jq,3]))
            plt.errorbar(dat[jq,0],dat[jq,1],yerr=dat[jq,2],fmt=".",markersize=7.7,color=cm(nvv),ecolor=cm(nvv),elinewidth=0.2, capsize=0,alpha=0.75)
        plt.ylabel(r"$\rm{Magnification}$",fontsize=18)
        plt.xlabel(r"$\rm{time}(\rm{days})$",fontsize=18)
        plt.title(
        r"$t_{\rm{E},r}\rm{(days)}=$"+str(round(tE,1))+
        r"$,~\log_{10}[\rho_{\star,r}]=$"+str(round(np.log10(ros1),2))+",~"+str(round(np.log10(ros2),2))+
        r"$,~\Delta t_{0}(\rm{days})=$"+str(round(Dt0,2))+r"$,~\Delta m(\rm{mag})=$"+str(round(Dmag,1))+
        r"$,~\mathcal{A}=$"+str(round(asym,1))+"\n"+
        r"$t_{\rm{E},b}\rm{(days)}=$"+str(round(tEf,1))+r"$\pm$"+str(round(dtE,1))+
        r"$,~\log_{10}[\rho_{\star,b}]=$"+str(round(np.log10(rosf),2))+r"$\pm$"+str(round(np.log10(dro),2))+
        r"$,~\chi^{2}_{\rm{r}}=$"+str(int(chiR))+r"$,~\chi^{2}_{\rm{b}}=$"+str(int(chi2m))+
        r"$,~\rm{dof}=$"+str(int(nd-5)),fontsize=13.5,color='darkblue')
        if(flag>0): pylab.ylim([ np.min(np.array([ymin1, ymin2, ymin3])),  np.max(np.array([ymax1, ymax2, ymax3]))   ])
        else:       pylab.ylim([ np.min(np.array([ymin1, ymin2])),  np.max(np.array([ymax1, ymax2]))   ])
        pylab.xlim([ np.min(mod[:,0]),  np.max(mod[:,0])])
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        ax1.legend(prop={"size":13.0})
        ###################################################
        dx=abs(np.max(mod[:,0]) - np.min(mod[:,0]) )
        dy=abs(np.max(np.array([ymin1 , ymin2]))-np.min(np.array([ymax1 , ymax2])) )
        x1=np.min(mod[:,0])
        y1=np.min(np.array([ymax1, ymax2]))
        ###################################################
        left, bottom, width, height=[ 0.1,  0.5,  abs(0.35),  abs(0.35) ]
        ax4 = fig.add_axes([ left, bottom, width, height ])
        ax4.plot(mod[:,1], mod[:,2], "k-",lw=1.5)#, label=r"$\rm{CM}-\rm{Lens}~\rm{Trajectory}$")
        ax4.plot(mod[:,3], mod[:,4], "r--", lw=1.5)#, label=r"$\rm{Star1}-\rm{Lens}~\rm{Trajectory}$")
        ax4.plot(mod[:,5], mod[:,6], "b:", lw=1.5)#, label=r"$\rm{Star2}-\rm{Lens}~\rm{Trajectory}$")
        ax4.plot(0.0,0.0, "k*", markersize=5.2)
        ax4.set_aspect('equal', adjustable='box')
        ax4.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax4.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        plt.axis('off')
        #####################################################
        fig=plt.gcf()
        fig.tight_layout()
        fig.savefig("./FigsLsst/LBinaryL{0:d}_{1:d}.jpg".format(int(ei),int(flag)), dpi=200)
        #######################################################################
        print("************ Light curve was plotted ***********************")    


print("faction of simplw events with Dm<2 && Sp<LsstA:", float(simple*100.0/nlens)) 
print("BS Microlensing events with Delta m<2:  ",         float(frac1*100.0/nlens))    
print("BS Microlensing events with sp<LSST_Accuracy:  ", float(frac2*100.0/nlens))    
print(nf, nlens, frac1, frac2, simple )
  



