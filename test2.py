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
import VBMicrolensing
vbb = VBMicrolensing.VBMicrolensing()
vbb.RelTol = 1e-04
vbb.Tol=1e-04
vbb.LoadESPLTable("./ESPL.tbl"); 
vbb.a1=0.0

u0=0.05
t0=0.0
tE=20.0
rho=0.1

tim=np.arange(-tE*1.2,   tE*1.2,  0.01)
u= np.sqrt( (tim-t0)*(tim-t0)/tE/tE + u0*u0 )
As=(u*u+2.0)/np.sqrt( u*u*(u*u+4.0) )
vbb.a1=0.0;
Atot=np.zeros(( len(tim) ))
for i in range(len(tim)):              
    Atot[i]=vbb.ESPLMag2(u[i],rho)

###############################################################################
####Full Width of Half Maximum  #####
N=int(len(tim))
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

##############################################
NM=1

plt.cla()
plt.clf()
fig=plt.figure(figsize=(8,6))
spec3=gridspec.GridSpec(nrows=3, ncols=1, figure=fig)
ax1= fig.add_subplot(spec3[:2,0])
ax1.plot(tim, As, "r--", lw=1.4, label=r"$\rm{PSPL}$")
ax1.plot(tim, Atot, "k-", lw=1.1,label=r"$\rm{ESPL}$")
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_ylabel(r"$\rm{Magnification}$", fontsize=17, labelpad=0.1)
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
ax2.plot(tim, As-Atot ,"k-", lw=1.5)
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
fig.savefig("./Figs/trhos{0:d}.jpg".format(NM),dpi=200)







