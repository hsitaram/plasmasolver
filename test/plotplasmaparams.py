import matplotlib.pyplot as plt
import numpy as np
from sys import argv
#===================================================================
def getspecinfo():

	str1="#NUMBER DENSITIES"
	infile=open('plasma_inputs','r')

	infile.readline()
	line=infile.readline()

	l=float(line.split()[1])

	for line in infile:
			if(line == "#NUMBER DENSITIES\n"):
				break
	specnames=[]

	for line in infile:
		specnames.append(line.split()[0])

	nions=0
	for i in range(len(specnames)):
		if(len(specnames[i].split('+')) > 1):
			nions=nions+1
		if(len(specnames[i].split('-')) > 1):
			nions=nions+1
	
	infile.close()

	nneutrals=len(specnames)-nions-1

	return(specnames,nions,nneutrals,l)
#===================================================================
def readfile(filename,nions,nspecies):

	infile=open(filename,'r')

	especnum=4
	nneutrals=nspecies-nions-1

	x        = np.array([])
	elecden  = np.array([])

	ionden   = np.array([])
	nden     = np.array([])
	potential= np.array([])
	efield   = np.array([])
	electemp = np.array([])
	jheating = np.array([])
	elasticol = np.array([])
	inelasticol = np.array([])
	electroncurr = np.array([])
	ioncurr = np.array([])

	for line in infile:

		spltline=line.split()
		x = np.append(x,float(spltline[0]))
		potential=np.append(potential,float(spltline[1]))
		efield=np.append(efield,float(spltline[2]))
		elecden=np.append(elecden,float(spltline[especnum]))

		for i in range(nions):
			ionden=np.append(ionden,float(spltline[especnum+i+1]))

		for i in range(nneutrals):
			nden=np.append(nden,float(spltline[especnum+nions+i+1]))

		electemp = np.append(electemp,float(spltline[especnum+nspecies+1]))
		jheating = np.append(jheating,float(spltline[especnum+nspecies+2]))
		elasticol = np.append(elasticol,float(spltline[especnum+nspecies+3]))
		inelasticol = np.append(inelasticol,float(spltline[especnum+nspecies+4]))
		electroncurr = np.append(electroncurr,float(spltline[especnum+nspecies+5]))
		ioncurr = np.append(ioncurr,float(spltline[especnum+nspecies+6]))

	return(x,elecden,ionden,nden,potential,efield,electemp,jheating,elasticol,inelasticol,electroncurr,ioncurr)
#===================================================================

fname=argv[1]

(specnames,nions,nneut,gapd)=getspecinfo()

(x,eden,iden,nden,pot,efield,etemp,ejheat,elcol,inelcol,electroncurr,ioncurr)=readfile(fname,nions,len(specnames))
x=x/gapd

xlen=6.0
ylen=3.0
me = 5
m=2
n=2
lt=3


npoints=len(x)
iondens=np.zeros((nions,npoints))
ndens=np.zeros((nneut,npoints))

for i in range(npoints):
	for j in range(nions):
		iondens[j][i]=iden[i*nions+j]

for i in range(npoints):
	for j in range(nneut):
		ndens[j][i]=nden[i*nneut+j]

(fig,ax)=plt.subplots(m,n,figsize=(n*xlen,m*ylen))
#fig.tight_layout()

i=0
j=0
ax[i][j].plot(x,eden,'k*-',markevery=me,markeredgecolor='black',linewidth=lt,label="E")
for n in range(nions):
	ax[i][j].plot(x,iondens[n],markevery=me+7,linewidth=lt,label=specnames[n+1])
ax[i][j].set_ylabel("number density (#/m3)")
lg=ax[i][j].legend(loc="best",fontsize=12)
lg.draw_frame(False)

i=0
j=1
ax[i][j].plot(x,pot,'r*-',markevery=me,markeredgecolor='red',linewidth=lt,label="Potential")
ax[i][j].set_ylabel("Voltage (V)")
ax_twin=ax[i][j].twinx()
ax_twin.plot(x,efield,'bo-',markevery=me+7,markeredgecolor='blue',linewidth=lt,label="Electric field")
ax_twin.set_ylabel("Electric field (V/m)")
lg=ax[i][j].legend(loc=1,fontsize=12)
lg.draw_frame(False)
lg=ax_twin.legend(loc=2,fontsize=12)
lg.draw_frame(False)

i=1
j=0
ax[i][j].plot(x,pot,'r*-',markevery=me,markeredgecolor='red',linewidth=lt,label="Potential")
ax[i][j].set_ylabel("Voltage (V)")
ax[i][j].set_xlabel("gap distance")
ax_twin=ax[i][j].twinx()
ax_twin.plot(x,etemp,'bo-',markevery=me+7,markeredgecolor='blue',linewidth=lt,label="Electron temperature")
ax_twin.set_ylabel("electron temperature (eV)")
lg=ax[i][j].legend(loc=1,fontsize=12)
lg.draw_frame(False)
lg=ax_twin.legend(loc=2,fontsize=12)
lg.draw_frame(False)

i=1
j=1
ax[i][j].plot(x,ejheat,'r*-',markevery=me,markeredgecolor='red',linewidth=lt,label="Joule heating")
ax[i][j].plot(x,elcol,'bo-',markevery=me,markeredgecolor='blue',linewidth=lt,label="Elastic collisions")
ax[i][j].plot(x,inelcol,'gs-',markevery=me,markeredgecolor='green',linewidth=lt,label="Inelastic collisions")
ax[i][j].set_xlabel("gap distance")
ax[i][j].set_ylabel("power density (W/m3)")
lg=ax[i][j].legend(loc="best",fontsize=12)
lg.draw_frame(False)

fig.tight_layout()
plt.savefig("plasmaparams.pdf")
plt.show()
