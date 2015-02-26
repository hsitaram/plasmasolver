import matplotlib.pyplot as plt
import numpy as np
from sys import argv
#===================================================================
def getspecinfo():

	str1="#NUMBER DENSITIES\n"
	infile=open('plasma_inputs','r')

	infile.readline()
	line=infile.readline()

	l=float(line.split()[1])

	for line in infile:
			if(line == str1):
				break
	specnames=[]

	for line in infile:
		specnames.append(line.split()[0])

	nbgspc=0
	for i in range(len(specnames)):
		if(specnames[i] == 'E'):
			break
	nbgspc=i
			
	nions=0
	for i in range(len(specnames)):
		if(specnames[i] == 'E'):
			especnum=i
		if(len(specnames[i].split('+')) > 1):
			nions=nions+1
		if(len(specnames[i].split('-')) > 1):
			nions=nions+1
	
	infile.close()

	nneutrals=len(specnames)-nions-nbgspc

	return(specnames,especnum,nbgspc,nions,nneutrals,l)
#===================================================================
def readfile(filename,enum,nbgspecies,nions,nneutrals,nspecies):

	infile=open(filename,'r')

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

	especnum=enum+3

	for line in infile:

		spltline=line.split()

		offset=0
		x = np.append(x,float(spltline[0]))
		
		offset=offset+1
		potential=np.append(potential,float(spltline[offset]))

		offset=offset+1
		efield=np.append(efield,float(spltline[offset]))
		
		for i in range(nbgspecies):
			offset=offset+1
			nden=np.append(nden,float(spltline[offset+i]))

		offset=offset+1
		elecden=np.append(elecden,float(spltline[offset]))
		
		for i in range(nions):
			offset=offset+1
			ionden=np.append(ionden,float(spltline[offset+i]))

		for i in range(nneutrals):
			offset=offset+1
			nden=np.append(nden,float(spltline[offset+i]))

		offset=offset+1
		electemp = np.append(electemp,float(spltline[offset]))
		
		offset=offset+1
		jheating = np.append(jheating,float(spltline[offset]))

		offset=offset+1
		elasticol = np.append(elasticol,float(spltline[offset]))
		
		offset=offset+1
		inelasticol = np.append(inelasticol,float(spltline[offset]))

		offset=offset+1
		electroncurr = np.append(electroncurr,float(spltline[offset]))

		offset=offset+1
		ioncurr = np.append(ioncurr,float(spltline[offset]))


	return(x,elecden,ionden,nden,potential,efield,electemp,jheating,elasticol,inelasticol,electroncurr,ioncurr)
#===================================================================

(specnames,enum,nbgspecies,nions,nneut,gapd)=getspecinfo()

fnames=[]

for i in range(len(argv)-1):
	fnames.append(argv[i+1])

xlen=6.0
ylen=3.0
me = 5
m=2
n=2
lt=3

(fig,ax)=plt.subplots(m,n,figsize=(n*xlen,m*ylen))

for fnum in range(len(fnames)):
	fname=fnames[fnum]

	(x,eden,iden,nden,pot,efield,etemp,ejheat,elcol,inelcol,electroncurr,ioncurr)=readfile(fname,enum,nbgspecies,nions,nneut,len(specnames))
	x=x/gapd

	npoints=len(x)
	iondens=np.zeros((nions,npoints))

	for i in range(npoints):
		for j in range(nions):
			iondens[j][i]=iden[i*nions+j]

	i=0
	j=0
	ax[i][j].plot(x,eden,markevery=me,linewidth=lt,label=fnames[fnum])
	ax[i][j].set_ylabel("electron density (#/m3)")
	lg=ax[i][j].legend(loc="best",fontsize=12)
	lg.draw_frame(False)

	i=0
	j=1
	ax[i][j].plot(x,iondens[0],markevery=me,linewidth=lt,label=fnames[fnum])
	#ax[i][j].plot(x,nden,markevery=me,linewidth=lt,label=fnames[fnum])
	ax[i][j].set_ylabel("ion density (#/m3)")

	i=1
	j=0
	ax[i][j].plot(x,pot,markevery=me,linewidth=lt,label=fnames[fnum])
	ax[i][j].set_ylabel("Voltage (V)")
	
	i=1
	j=1
	ax[i][j].plot(x,etemp,markevery=me,linewidth=lt,label=fnames[fnum])
	ax[i][j].set_ylabel("Electron temp (eV)")

	print "total current:",electroncurr[-1]+ioncurr[-1]

fig.tight_layout()
#plt.savefig("plasmaparams.pdf")
plt.show()
