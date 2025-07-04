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


sm=1.0e-50
Rsun=6.957*pow(10.0,8.0); #solar radius [meter]
AU= 1.4960*pow(10.0,11.0);
Tobs =float(7*30.0)
LsstAAC=float(0.7*1000.0)## mili arcsec LSST 

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
     r"$\chi^{2}_{\star,1}$", r"$\chi^{2}_{\star,2}$",r"$\chi^{2}_{\rm{real}}$",r"$\chi^{2}_{\rm{CM}}$",r"$f_{b}$", r"$n_{\rm{data}}$", r"$Flag_{3}$"]


nam2=[r"$Nsim$", r"$\rm{filter}$" , r"$m_{\rm{base}}(\rm{mag})$", r"$N_{\rm{blend}}$", r"$f_{\rm{b}}$", r"$M_{\rm{ab},~1}(\rm{mag})$", r"$m_{\rm{ap}, 1}(\rm{mag})$", r"$M_{\rm{ab}, 2}(\rm{mag})$", r"$m_{\rm{ap}, 2}(\rm{mag})$", r"$\rm{Extinction}(\rm{mag})$"]#10
    
###############################################################################

frac1=0.0
frac2=0.0
frac3=0.0
simple=int(0)
nlens=int(0)
nc  =int(52)
f2=open("./files/BestfitL.dat","r")
na= sum(1 for line in f2)  
ana=np.zeros((na, 17))
f1=open("./files/resGBLsstA.dat","r")
nf= sum(1 for line in f1)
par=  np.zeros((nf,nc))
detec=np.zeros((nf,nc))
par2=  np.zeros((nf*6,10))
detec2=np.zeros((nf*6,10))
par= np.loadtxt("./files/resGBLsstA.dat")  
par2=np.loadtxt("./files/resGBLsstB.dat")  
ana= np.loadtxt("./files/BestfitL.dat")  
nn=int(0)
good=np.zeros((3000,16))
magb=np.zeros((6)); Map1=np.zeros((6)); Map2=np.zeros((6)); 


save2=open("./Table2.dat","w")
save2.close()
###############################################################################

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
    print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH", ei)

    ###########################################################################

    if(dchiL>=300.0 and peak>0 and threef>0 ):#and magb[2]<24.3 and magb[2]>16.0
        detec[int(nlens),:]=par[i,:]
        for k in range(6): 
            detec2[int(nlens*6+k),:]=par2[int(i*6+k),:]
        nlens+=1
        if(Dmag<2.0):   frac1+=1.0
        if(cri2<=1.0):  frac2+=1.0 
        if(Dmag<2.0 and cri2<=1.0):  
            frac3+=1.0
            nn=-1
            for k in range(na):   
                if(ana[k,0]==ei): 
                    nn=k; 
                    break; 
            if(nn<0):                 
                print("It seems that there is no fitted model for this event ,  ", ei)
            if(nn>=0):         
                num,u0f,eu01, eu02, t0f, et01, et02, tEf, etE1, etE2, fbf, efb1, efb2, rosf, eros1, eros2, chi2m= ana[nn,:]
                etE= 0.5*(abs(etE1)+abs(etE2))
                eros=0.5*(abs(eros1)+abs(eros2))
                if(num!=ei or tEf<=0.0 or rosf<=0.0  or u0f<=0.0  or fbf<=0.0 or etE<0.0 or eros<=0.0 or chi2m<0.0 or fb>1.0): 
                    print("Error, num, ei",   num, ei ,  tEf, u0f, fbf,  rosf, eros, etE, chi2m, fb)
                    input("Enter a number  ")
                if(abs(chi2m-chiR)<100 and abs(tE-tEf)<50.0): 
                    if(f21<1.0): ros=ros1
                    else:        ros=ros2
                    good[simple,:]=np.array([tE,fb,u0,ros,Dmag,sproj, magb[2], tEf,rosf,fbf,u0f,etE, eros, ei,2.0*abs(tE-tEf)/tE, q])#16
                    simple+=1; 
                    if(abs(tE-tEf)>100.0): 
                        print (tE, tEf, ei, chi2m, chiR) 
                        input("Enter a number   ")

###############################################################################
print ("fraction of simulated Lening events with Dmag<2:  ",  frac1*100.0/nlens)
print ("fraction of simulated lensing events with Sp<LSSTAA: ", frac2*100.0/nlens)
print ("fraction of simulated lensing events with Sp<LSSTAA && Dmag<2: ", frac3*100.0/nlens)
print ("fraction of them which looks simple:  ",  simple*100.0/frac3)
print ("****************************")


print("<tE>: ",   np.mean(good[:simple,0]))
print("<fb>: ",   np.mean(good[:simple,1]))
print("<u0>: ",   np.mean(good[:simple,2]))
print("<ros>: ",  np.mean(good[:simple,3]))
print("<Dmag>: ", np.mean(good[:simple,4]))
print("<sp>: ",   np.mean(good[:simple,5]))
print("<mbase>: ",np.mean(good[:simple,6]))
print("<tEf>: ",  np.mean(good[:simple,7]))
print("<rosf>: ", np.mean(good[:simple,8]))
print("<fbf>: ",  np.mean(good[:simple,9]))
print("<u0f>: ",  np.mean(good[:simple,10]))
print("<etE>: ",  np.mean(good[:simple,11]))
print("<eros>: ", np.mean(good[:simple,12]))
print("\Delta tE: ", np.mean( good[:simple,0]-good[:simple,7] ) )
print("<sigma_m_l>: ", np.mean(good[:simple,14]))
print("<q>: ", np.mean(good[:simple,15]) )

sav=np.array([frac3*100.0/nlens, simple*100.0/nlens, 
np.mean(good[:simple,0]), np.mean(good[:simple,3]), np.mean(good[:simple,2]),
np.mean(good[:simple,7]), np.mean(good[:simple,8]), np.mean(good[:simple,10]), np.mean(good[:simple,4]), np.mean(good[:simple,14])])#10
save2=open("./Table2.dat","a+") 
np.savetxt(save2,sav.reshape((-1,10)),
fmt="LSST & $%.2f$ & $%.2f$ & $%.2f$ & $%.4f$ & $%.2f$ & $%.2f$ & $%.4f$ & $%.2f$ & $%.2f$ & $%.3f$")    
save2.close()
    
# epsilon_1, epsilon_2,  tE, ros, u0,  tEf, rosf, u0f,  Dmag, sigma_ml
    
###############################################################################

inip=[r"$t_{\rm{E}}(\rm{days})$", r"$f_{\rm{b}}$", r"$u_{0}$", r"$\rho_{\star}$", r"$\Delta m(\rm{mag})$", r"$s_{\rm{p}}(\rm{mas})$", r"$m_{\rm{base}}(\rm{mag})$"]##7

for i in range(7): 
    plt.cla()
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    plt.plot(good[:int(simple),i], good[:int(simple),0]-good[:int(simple),7], "r*", markersize=5.0)
    plt.xlabel(str(inip[i]),fontsize=18)
    plt.ylabel(r"$\Delta t_{\rm{E}}(\rm{days})$",fontsize=18)
    #plt.xlim([ np.min(par[:int(simple),i]) , np.max(par[:int(simple),i]) ])
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    if(i==3): plt.xscale('log')
    plt.grid(True)
    plt.grid(linestyle='dashed')  
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./FigsLsst/Histo/DtE{0:d}.jpg".format(i), dpi=200)             
    ###########################################################################                      
    plt.cla()
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    plt.plot(good[:int(simple),i], np.abs(good[:int(simple),3]-good[:int(simple),8]), "r*", markersize=5.0)
    plt.xlabel(str(inip[i]),fontsize=18)
    plt.ylabel(r"$|\Delta \rho_{\star}|$",fontsize=18)
    #plt.xlim([ np.min(par[:int(simple),i]) , np.max(par[:int(simple),i]) ])
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.yscale('log')
    if(i==3): plt.xscale('log')
    plt.grid(True)
    plt.grid(linestyle='dashed')  
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./FigsLsst/Histo/Dros{0:d}.jpg".format(i), dpi=200) 


print("tE",  np.min(good[:int(simple),0]), np.min(good[:int(simple),7]), np.max(good[:int(simple),0]), np.max(good[:int(simple),7]))
print("ros", np.min(good[:int(simple),3]), np.min(good[:int(simple),8]), np.max(good[:int(simple),3]), np.max(good[:int(simple),8]))
plt.cla()
plt.clf()
fig= plt.figure(figsize=(8,6))
plt.scatter(np.log10(good[:int(simple),0]), np.log10(good[:int(simple),3]),color='r',marker='o',edgecolors='darkred',s=27.0, alpha=0.65)
plt.scatter(np.log10(good[:int(simple),7]), np.log10(good[:int(simple),8]),color='g',marker='s',edgecolors='darkgreen',s=24.0, alpha=0.6)
plt.xlabel(r"$\log_{10}[t_{\rm{E}}(\rm{days})]$",fontsize=18)
plt.ylabel(r"$\log_{10}[\rho_{\star}]$",fontsize=18)
plt.text(1.85,-1.,r"$\rm{Real}~\rm{values}$", fontsize=18, color="r", alpha=1.0)
plt.text(1.85,-1.3,r"$\rm{Best}-\rm{fitted}~\rm{values}$", fontsize=18, color="g", alpha=1.0)
plt.title(r"$\rm{Simple}-\rm{like}~\rm{and}~\rm{BS}~\rm{Microlensing},~\rm{LSST}$", fontsize=18, color="k")
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.xlim([np.log10(3.8),     np.log10(275.0)])
plt.ylim([np.log10(1.5e-5), np.log10(0.15)])
#plt.yscale('log')
#plt.xscale('log')
#plt.grid(True)
#plt.grid(linestyle='dashed')  
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./FigsLsst/Histo/DmapL.jpg", dpi=200)                   
            
###############################################################################
for i in range(nc):
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    ax= plt.gca()              
    plt.hist(par[:,i],30,histtype='bar',ec='darkgreen',facecolor='green',alpha=0.5,rwidth=1.5,label=r"$\rm{Simulated}~\rm{events}$")
    plt.hist(detec[:int(nlens),i],30,histtype='bar',ec='darkred',facecolor='red',alpha=0.5,rwidth=1.5,label=r"$\rm{Detected}~\rm{events}$")
    y_vals = ax.get_yticks()
    ax.set_yticks(y_vals)
    ax.set_yticklabels(['{:.2f}'.format(float(1.0*x*(1.0/nf))) for x in y_vals]) 
    y_vals = ax.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
    ax.set_xlabel(str(nam[i]),fontsize=19,labelpad=0.1)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.legend(prop={"size":12.5})
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig=plt.gcf()
    fig.savefig("./FigsLsst/Histo/histL{0:d}.jpg".format(i),dpi=200)
print ("****  All histos are plotted ***********************  ")   
################################################################################

################################################################################
for j in range(6): 
    for i in range(10):
        plt.cla()
        plt.clf()
        fig= plt.figure(figsize=(8,6))
        ax= plt.gca()              
        plt.hist(  par2[j:int(nf*5):5,i],30,histtype='bar',ec='k',facecolor='green',alpha=0.5,rwidth=1.5,label=r"$\rm{Sim}~\rm{events}$")
        plt.hist(detec2[j:int(nlens*5):5,i],30,histtype='bar',ec='darkred',facecolor='red',alpha=0.5,rwidth=1.5,label=r"$\rm{Detected}~\rm{events}$")
        y_vals = ax.get_yticks()
        ax.set_yticks(y_vals)
        ax.set_yticklabels(['{:.2f}'.format(float(1.0*x*(1.0/nf))) for x in y_vals]) 
        y_vals = ax.get_yticks()
        plt.ylim([np.min(y_vals), np.max(y_vals)])
        ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
        ax.set_xlabel(str(nam2[i])+r"$\rm{Filter}=$"+str(j), fontsize=19,labelpad=0.1)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        plt.legend(prop={"size":12.5})
        plt.grid("True")
        plt.grid(linestyle='dashed')
        fig=plt.gcf()
        fig.savefig("./FigsLsst/Histo/histLB{0:d}_{1:d}.jpg".format(i,j),dpi=200)
print ("****  All histos(2) are plotted ***********************  ")   
################################################################################

            
            
            
            
            
            
            
            
            
