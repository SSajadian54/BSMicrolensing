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
from matplotlib.ticker import FormatStrFormatter
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
#import VBMicrolensing
#vbb = VBMicrolensing.VBMicrolensing()
#vbb.RelTol = 1e-04
#vbb.Tol=1e-04
#vbb.LoadESPLTable("./ESPL.tbl"); 
#vbb.a1=0.0

t01=0.0
sigm=0.004## for WFIRST and for a star with the magnitude 19.5 
siga=abs(pow(10.0,-0.4*sigm)-1.0)
###############################################################################
def model(param, flag, NM):
    tE, t02, u01, u02, a= param[0], param[1], param[2], param[3], param[4]
    tim= np.arange(-tE*2.1 , tE*2.1, float(15.0/60.0/24.0) )
    N=int(len(tim))
    dis=np.zeros((N))
    u1= np.zeros((N))
    u2= np.zeros((N))
    As1=np.zeros((N))
    As2=np.zeros((N))
    Atot=np.zeros((N))
    u1=np.sqrt(u01*u01+(tim-t01)*(tim-t01)/tE/tE)
    u2=np.sqrt(u02*u02+(tim-t02)*(tim-t02)/tE/tE)
    As1= (u1*u1+2.0)/np.sqrt(u1*u1*(u1*u1+4.0)) #+alfa)/(1.0+alfa) 
    As2= (u2*u2+2.0)/np.sqrt(u2*u2*(u2*u2+4.0))#*alfa)/(1.0+alfa)
    Atot=(As1+As2*a)/(1.0+a)
    ###########################################################################
    ####Full Width of Half Maximum  #####
    Am=abs(np.max(Atot-1.0)*0.5)
    dis=np.abs(Atot-(Am+1.0))
    n1=np.argmin(dis[:int(N/2)])
    t1= float(tim[n1])
    n2=np.argmin(dis[int(N/2):])
    t2= float(tim[n2+int(N/2)])
    FWHM=abs(t2-t1)
    print("FWHM: " ,  t1, t2, n1, n2+int(N/2), FWHM)
    ###########################################################################
    ####  U0_fitted && t0_fitted ####
    u0=np.arange(0.00001, 1.0, 0.0001)
    AA=(u0*u0+2.0)/np.sqrt(u0*u0*(u0*u0+4.0))
    n3=np.argmin(np.abs(AA-np.max(Atot)))
    u0f=u0[n3]
    t0f=tim[np.argmax(Atot)]
    nt0=int(np.argmax(Atot))
    print("u0f, t0f, N: ", u0f, t0f, n3, AA[n3], np.max(Atot), nt0 )
    ###########################################################################
    ### tE_fitted ####
    tEr=np.arange(tE-7.0, tE+7.0 , 0.001)
    u=np.sqrt( ((t2-t0f)/tEr)**2.0+ u0f*u0f)
    A=np.abs((u*u+2.0)/np.sqrt(u*u*(u*u+4.0))-(Am+1.0))
    n4=np.argmin(A)
    tEf=tEr[n4] 
    print("tEf: ", tEf, n4, A[n4], Am+1.0)   
    ###########################################################################
    ####Fitted_model####
    ufit =np.sqrt(u0f*u0f+(tim-t0f)*(tim-t0f)/tEf/tEf)
    Afit =(ufit*ufit+2.0)/np.sqrt(ufit*ufit*(ufit*ufit+4.0))
    chi2n=np.sum((Afit-Atot)**2.0/siga/siga)/N
    Delta=np.sqrt(np.sum((Afit-Atot)**2.0)/N)
    print("chi2n: ", chi2n, siga , Delta )
    asym=0.0; num=0.0
    Nh=int(N-nt0)
    for i in range(Nh):
        nr=int(nt0+i)
        nl=int(nt0-i)
        if(nr<N and nl>=0):  
            asym+=(Atot[nr]-Atot[nl])**2.0;
            num +=1.0 
    asym=np.sqrt(asym/num)
    print("Asymmetry:  ",  asym,    siga,  num)
    ###########################################################################
    if(flag==1):  
        plt.cla()
        plt.clf()
        fig=plt.figure(figsize=(8,6))
        spec3=gridspec.GridSpec(nrows=3, ncols=1, figure=fig)
        ax1= fig.add_subplot(spec3[:2,0])
        ax1.plot(tim, (As1+a)/(1.0+a), "r--", lw=1.4, label=r"$\rm{Star}~1$")
        ax1.plot(tim, (1.0+a*As2)/(1.0+a), "b:", lw=1.4, label=r"$\rm{Star}~2$")
        ax1.plot(tim, Atot, "k-", lw=1.1,label=r"$\rm{Binary}~\rm{Source}$")
        ax1.plot(tim, Afit, "g-." , lw=1.4, label=r"$\rm{Best}-\rm{Fitted}~\rm{Model}$")
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax1.set_ylabel(r"$\rm{Magnification}$", fontsize=17, labelpad=0.1)
        ax1.set_title(
        r"$t_{\rm{E}, r}(\rm{days})=$"+str(round(tE,2))+ 
        r"$,~\Delta t_{0}=$"+str(round(t02-t01,2))+
        r"$,~u_{0, 1}=$"+str(round(u01,1))+r"$,~u_{0, 2}=$"+str(round(u02,1))+
        r"$,~\mathcal{F}=$"+str(round(a,2))+"\n"+
        r"$t_{\rm{E}, f}(\rm{days})=$"+str(round(tEf,2))+
        r"$,~t_{0, f}=$"+str(round(t0f,2))+
        r"$,~u_{0, f}=$"+str(round(u0f,2))+
        r"$,~\mathcal{A}=$"+str(round(asym,3))+
        r"$,~\Delta=$"+str(round(Delta,3)), fontsize=16, color="darkblue")
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        ax1.set_xlim([ np.min(tim), np.max(tim) ])
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.grid("True")
        ax1.grid(linestyle='dashed')
        plt.legend()
        plt.legend(loc=1,fancybox=True, shadow=True, fontsize=15, framealpha=0.85)
        #######################################
        ax2= fig.add_subplot(spec3[2,0],sharex=ax1) 
        ax2.plot(tim, Atot-Afit ,"k-", lw=1.5)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax2.set_xlabel(r"$\rm{time}(\rm{days})$",fontsize=17, labelpad=0.1)
        ax2.set_ylabel(r"$\Delta A$",  fontsize=17,labelpad=0.1)
        ax2.set_xlim([np.min(tim), np.max(tim)])
        ax2.grid("True")
        ax2.grid(linestyle='dashed')
        plt.subplots_adjust()
        #######################################    
        fig=plt.gcf()
        fig.savefig("./Figs/test{0:d}.jpg".format(NM),dpi=200)
    param[5:]=np.array([tEf, u0f, t0f, Delta, asym, chi2n, FWHM ])    
    return(param)
###############################################################################

NM=1
### parameters 
save2=open("./Output.txt","w")
save2.close()
param=np.zeros((12))

for i in range(500000): 
    tE= float(np.random.rand(1)*45.0)+5.0
    a=  float(np.random.rand(1))
    u01=float(np.random.rand(1))
    u02=float(np.random.rand(1))
    t02=float(np.random.rand(1)*0.1*tE)
    NM+=1
    param[:5]=np.array([tE, t02, u01, u02, a])
    flag=0
    if(i<100): flag=1
    param=model(param,flag,NM);      #tEf, u0f, t0f, Delta, asym, chi2n, FWHM
    save2=open("./Output.txt","a+")
    np.savetxt(save2,param.reshape(-1,12),fmt='%.4f %.4f  %.4f %.4f  %.4f  %.6f  %.6f  %.6f  %.6f   %.7f   %.3f   %.3f') 
    save2.close()
    print("*********************************:  ", NM)
    
    
'''    
for i in range(5):
    tE=np.random.rand(1)*0.1*tE)float(10.0+i*10.0)
    for j in range(5):  
        a=float(0.05+j*0.23)
        for k in range(5):  
            u01=float(0.1+k*0.21)
            for l in range(5): 
                u02=float(0.1+l*0.21) 
                for ii in range(10):  
                    NM+=1
                    t02=float(np.random.rand(1)*0.1*tE)
                    print("Initial parameters:  ", tE, a, u01, u02, t02)
'''                    
##############################################





