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
Error=0.004 
ndim=12
nam=[r"$t_{\rm{E}}(\rm{days})$", r"$\Delta t_{0}(\rm{days})$", r"$u_{0, 1}$", r"$u_{0, 2}$", r"$\mathcal{F}$", 
r"$t_{\rm{E}, f}(\rm{days})$", r"$u_{0, f}$", r"$t_{0, f}(\rm{days})$", r"$\Delta$", r"$\mathcal{A}$", r"$\chi^{2}_{\rm{n}}$", 
r"$\rm{FWHM}(\rm{days})$", r"$\Delta t_{\rm{E}}(\rm{days})$"]#13


nq=10
par= np.zeros((nq,4,5))
#Dt0=np.zeros((10,3))
#u01=np.zeros((10,3))
#u02=np.zeros((10,3))
#fra=np.zeros((10,3))

par[:,0,0]=np.array([5.0, 10, 15, 20,  25, 30,  35, 40,  45, 50])#10  tE,
par[:,0,1]=np.array([0.0,0.5, 1.0,1.5, 2.0,2.5, 3.0,3.5, 4.0, 5.0])#10,   Dt0
par[:,0,2]=np.array([0.0,0.15 ,0.3, 0.4,0.5,0.6,0.7,0.8, 0.9, 1.0])#10,   u01
par[:,0,3]=np.array([0.0,0.15 ,0.3, 0.4,0.5,0.6,0.7,0.8, 0.9, 1.0])#10,   U02
par[:,0,4]=np.array([0.0,0.2 ,0.3, 0.4,0.5,0.6,0.7,0.8, 0.9, 1.0])#10,   flux_ratio

#Effi= np.zeros((nq))
###############################################################################
def func(val, par, k, DtE, flag): 
    nt=-1
    if(val<par[0,0,k]):            nt=0
    elif(val>par[int(nq-1),0,k]):  nt=int(nq-1)
    else: 
        for j in range(nq-1): 
            if(float((val-par[j,0,k])*(val-par[j+1,0,k]))<=0.0):
                nt=j
                break 
    if(nt<0):  
        print("Big Error nt is negative:  ",  nt)
        input("Enter a number ")
    par[nt,1,k]+= DtE*flag# variation of tE
    par[nt,2,k]+= 1.0*flag# PSPL-like ones    
    par[nt,3,k]+= 1.0     # All simulated events
###############################################################################
#tE, t02, u01, u02, a, tEf, u0f, t0f, Delta, asym, chi2n, FWHM 
#0   1     2    3   4   5     6   7     8     9     10     11  

f2=open('./Output.txt')
nr=int(len(f2.readlines()))
res=np.zeros(( nr,ndim))
res=np.loadtxt('./Output.txt')


save2=open("./Table1.txt","w")
save2.close()

sigm= np.array([
0.001000787266916262, 0.0010504179299388695, 0.0013164612912327737,  0.0024082589099197684, 0.005983255916403701,  0.010087137650944121, 0.017062880789085626,  0.062256016296805036, 0.14414985937062885,  0.34839492983686626,  0.01, 0.02, 0.04, 0.1])#10+4
mapp= np.array([12.155057835273517,  14.217933920645809, 16.28356516962735, 18.216814396782496, 20.282202170979044 , 21.31509841431101, 22.219089102671028, 24.020895488076746, 25.01330910780286, 26.01505999779293, 21.0, 22.0, 23.0, 24.0])#10+4
siga=np.abs(np.power(10.0,-0.4*sigm)-1.0)

#sigm=np.array([])
#mapp=np.array([])
#siga=np.abs(np.power(10.0,-0.4*sigm)-1.0)


for i in range(len(siga)): 
    Error=float(siga[i])  
    res2=np.zeros((nr , ndim))
    k=0
    for j in range(nr):  
        if(res[j,8]<Error and res[j,9]<Error): 
            res2[k,:]=res[j,:]
            k+=1
    arr=np.array([ np.log10(Error), np.log10(sigm[i]),  mapp[i],  float(k*100.0/nr),  np.mean(res2[:k,0]),  np.mean(res2[:k,2]), np.mean(res2[:k,3]), np.mean(res2[:k,1]), np.mean(res2[:k,4]), np.log10(np.mean(res2[:k,8])), np.log10(np.mean(res2[:k,9])), np.mean(np.abs(res2[:k,0]-res2[:k,5])), 2.0*np.mean(np.abs(res2[:k,0]-res2[:k,5])/res2[:k,0]) ])
    save2=open("./Table1.txt","a+")       
    np.savetxt(save2, arr.reshape((-1,13)),
    fmt="$%.3f$ & $%.3f$ & $%.1f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$")  
    save2.close()      
            
            
###############################################################################            
res2=np.zeros((nr , ndim)) 
Error=float(siga[3])             
k=0
for i in range(nr): 
    flag=0
    if(res[i,8]<Error and res[i,9]<Error ): 
        res2[k,:]=res[i,:]
        k+=1
        flag=1
    func(float(res[i,0]),par, 0, abs(res[i,0]-res[i,5]), flag )
    func(float(res[i,1]),par, 1, abs(res[i,0]-res[i,5]), flag )
    func(float(res[i,2]),par, 2, abs(res[i,0]-res[i,5]), flag )
    func(float(res[i,3]),par, 3, abs(res[i,0]-res[i,5]), flag )
    func(float(res[i,4]),par, 4, abs(res[i,0]-res[i,5]), flag )
for i in range(5):
    for j in range(nq): 
        par[j,1,i]=float(par[j,1,i]/(par[j,2,i]+0.0001))# <Delta tE>  
        par[j,2,i]=float(par[j,2,i]*100.0/(par[j,3,i]+0.0001))#Efficiency    
print("number of unresolvable events:  ", k, float(k*100.0/nr))

###############################################################################
# <Delta tE> versus 5 parameters
for i in range(5): 
    plt.cla()
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    plt.plot(par[:(nq-1),0,i], par[:(nq-1),1,i], "ro")
    plt.xlabel(str(nam[i]),fontsize=18)
    plt.ylabel(str(nam[12]),fontsize=18)
    plt.xlim([ np.min(par[:,0,i]) , np.max(par[:,0,i]) ])
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.grid(True)
    plt.grid(linestyle='dashed')  
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./scatterp/StepDT{0:d}.jpg".format(i), dpi=200)
    
###############################################################################    
# fraction f PSPL-like events to toatl simulated ones     
for i in range(5): 
    plt.cla()
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    plt.plot(par[:int(nq-1),0,i], par[:int(nq-1),2,i], "ro")
    plt.xlabel(str(nam[i]),fontsize=18)
    plt.ylabel(r"$\varepsilon[\%]$",fontsize=18)
    plt.xlim([ np.min(par[:,0,i]) , np.max(par[:,0,i]) ])
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.grid(True)
    plt.grid(linestyle='dashed')  
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./scatterp/Effi{0:d}.jpg".format(i), dpi=200)    
###############################################################################
# Histograms 
for i in range(ndim):
    plt.cla()
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    ax= plt.gca()  
    plt.hist( res[:,i],40,histtype='bar',ec='darkgreen',facecolor='g',alpha=0.6, rwidth=1.0,label=r"$\rm{All}~\rm{BSPL}~\rm{events}$")
    plt.hist(res2[:k,i],40,histtype='bar',ec='darkred',  facecolor='red',alpha=0.5,rwidth=1.0,label=r"$\rm{PSPL}-\rm{like}~\rm{events}$")
    y_vals =ax.get_yticks()
    ax.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/nr)) for x in y_vals]) 
    y_vals =ax.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    ax.set_xlabel(str(nam[i]),fontsize=19,labelpad=0.1)
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.legend()
    plt.legend(loc='upper right',fancybox=True, shadow=True)
    plt.legend(prop={"size":16})
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./scatterp/Histo{0:d}.jpg".format(i), dpi=200)    
    print("Histogram *************, ", i, str(nam[i]), "\n"  )
###############################################################################
# Scatter plot 
for i in range(5):  # Delta, asym, chi2n,  DtE
    for j in range(5):# tE, t02, u01, u02, a  
        plt.clf()
        plt.cla()
        fig=plt.figure(figsize=(8,6))    
        if(i!=4): plt.scatter(res2[:k,j] , res2[:k,i+8], s=10.0,marker="o", color="g")
        if(i==4): plt.scatter(res2[:k,j], np.abs(res2[:k,0]-res2[:k,5]), s=10.0,marker="o", color="g" )
        plt.ylabel(str(nam[i+8]),fontsize=18)
        plt.xlabel(str(nam[j]),fontsize=18)
        #plt.xlim([0.0,10.0])
        #plt.yscale('log')
        plt.xticks(fontsize=18, rotation=0)
        plt.yticks(fontsize=18, rotation=0)
        plt.grid(True)
        plt.grid(linestyle='dashed')  
        fig=plt.gcf()
        #fig.tight_layout()
        fig.savefig("./scatterp/splot{0:d}_{1:d}.jpg".format(i,j), dpi=200)
    print("Scatter plot ************, ", i, str(nam[i+8]) )    

###############################################################################    
    
    


