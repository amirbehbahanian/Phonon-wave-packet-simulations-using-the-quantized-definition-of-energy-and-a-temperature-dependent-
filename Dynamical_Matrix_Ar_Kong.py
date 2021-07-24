    # -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 12:07:56 2019

@author: A02197412
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import integrate
import copy

f=open('Amir.log')
x=f.read()
# =============================================================================
x = x.split('\n')
startingindex = [i for i,a in enumerate(x) if a=='# qx\t qy \t qz \t\t Phi(q)']

y1 = x[startingindex[-1]+1:startingindex[-1]+1+64000]


y1 = [y1[a].split() for a in range(len(y1))]


y1 = [[float(y1[a][b]) for b in range(len(y1[a]))] for a in range(len(y1))]


df = pd.DataFrame(y1)


del(x, y1)

Kpoints = df.index[(df[0]-df[1]<=0) & (df[0]+df[1]+df[2]-1.5<=0) & (df[1]-df[2]>=0) & (df[0]-1<=0)]
# =============================================================================
col=int((len(df.columns)-3)/2)

nacell=1    #number of atoms per cell
dim=3
matdim=nacell*dim
Numatoms=1
convling = 15.63312493
conv = (1.60218e-19/(1e-20*(Numatoms*1e-3)/(6.022e23)))

for j in range(len(df)):
    z=[]
    for i in range(col):
        z.append(complex(df.iloc[j,2*i+3],df.iloc[j,2*i+4]))    #Creates the complex numbers from the imaginary and the real part 
    z=np.transpose(z)  
    D=np.vstack((z[0:matdim],z[matdim:2*matdim],z[2*matdim:3*matdim])) #Stacks the values to creat a 3*3 Matrix
    D=D* conv       #converting the ev/(A^2*(grams/mol))  to J/(m^2*Kg) 
    
    globals()['w%s' %(j)]=np.linalg.eig(D)
# =============================================================================
# 
GX=df.loc[(df[0]<=0.5) & (df[1]==0) & (df[2]<=0.5) & (df[0]==df[2])].index                         # GAMMA TO X
XW=df.loc[(df[0]==0.5) & (df[1]<=0.25) & (df[2]<=0.75) & (df[2]>=0.5) & (df[2]-df[1]==0.5)].index
KG=df.loc[(df[1]<=1) & (df[0]==df[2]) & (df[2]<=0.5) & (df[2]+df[0]==df[1])].index
GL=df.loc[(df[0]==df[1]) & (df[1]==df[2]) & (df[2]<=0.5)].index
  
 
freq_GX=[]
for i in range(round(int(len(GX)))):
    freq_GX.append(np.sqrt(globals()['w%s' %(GX[i])][0])/(2*np.pi))
freq_XW=[]
for i in range(round(int(len(XW)))):
    freq_XW.append(np.sqrt(globals()['w%s' %(XW[i])][0])/(2*np.pi))
freq_KG=[]
for i in range(round(int(len(KG)))):
    freq_KG.append(np.sqrt(globals()['w%s' %(KG[i])][0])/(2*np.pi))
freq_GL=[]
for i in range(round(int(len(GL)))):
    freq_GL.append(np.sqrt(globals()['w%s' %(GL[i])][0])/(2*np.pi))

# ############### Separates different branches of the Disp. Rel. ################
 
for j in range(3):
    globals()['freq_GX%s' %(j+1)]=[]
    for i in range(len(GX)):
        globals()['freq_GX%s' %(j+1)].append(freq_GX[i][j])

for j in range(3):
    globals()['freq_XW%s' %(j+1)]=[]
    for i in range(len(XW)):
        globals()['freq_XW%s' %(j+1)].append(freq_XW[i][j])

for j in range(3):
    globals()['freq_KG%s' %(j+1)]=[]
    for i in range(len(KG)):
        globals()['freq_KG%s' %(j+1)].append(freq_KG[i][j])
for j in range(3):
    globals()['freq_GL%s' %(j+1)]=[]
    for i in range(len(GL)):
        globals()['freq_GL%s' %(j+1)].append(freq_GL[i][j])

# ========================= Branch Data Correction ============================
"""KG"""
a = copy.copy(freq_KG1[10:16])
b = copy.copy(freq_KG3[10:16])
freq_KG1[10:16] = copy.copy(b)
freq_KG3[10:16] = copy.copy(a)

a = copy.copy(freq_KG1[16:21])
b = copy.copy(freq_KG2[16:21])
freq_KG1[16:21] = copy.copy(b)
freq_KG2[16:21] = copy.copy(a)

# ========================= Density of States =================================
"""
Dont forget that the frequencies are in angular form and they are aquared then f=sqrt(w)/(2*pi)
"""
freq = []
for i in range(len(df)):
    freq.append(globals()['w%s' %i][0])
    del(globals()['w%s' %i])
freq = pd.DataFrame(freq)
freq = freq.apply(lambda x:np.real(np.sqrt(x)/(2*np.pi)))
freq = freq.iloc[Kpoints]

Freq=np.arange(0,2.001e12,2e12/100)
data = []
for i in range(len(Freq)):
    try:
        data.append(sum(freq[(freq>=Freq[i]) & (freq<Freq[i+1])].count()))
    except IndexError:
        data.append(0)
data = pd.DataFrame(data)
data.insert(0,'Freq',Freq)
data.columns = ['Freq', 'DOS']
x = np.array(data['Freq'])
y = np.array(data['DOS'])
integ = integrate.simps(y,x)
data['DOS'] = data['DOS'].apply(lambda x: x/(integ/3))

# ========================= Density of States END =============================

fig = plt.figure('Disp. Relation and DOS',figsize=(7.2, 4))
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.98, top=0.98, wspace=0, hspace=0)

grid = plt.GridSpec(3, 4, hspace=0.05, wspace=0.05)

disp_rel = fig.add_subplot(grid[0:3, 0:3])
DOS = fig.add_subplot(grid[0:3, 3:], sharey=disp_rel)

color='^k'
markersize=2.5
fontsize=12
FN= 'Times New Roman'

disp_rel.plot(range(0,len(GX)),[10**-12*x for x in freq_GX1],color, label= '10K', ms=markersize)
disp_rel.plot(range(0,len(GX)),[10**-12*x for x in freq_GX2],color, ms=markersize)
disp_rel.plot(range(0,len(GX)),[10**-12*x for x in freq_GX3],color, ms=markersize)
disp_rel.plot(range(len(GX), len(GX)+len(XW)), [10**-12*x for x in freq_XW], color, ms=markersize)
disp_rel.plot(range(len(GX)+len(XW), len(GX)+len(XW)+len(XW)), [10**-12*x for x in freq_XW[::-1]], color, ms=markersize)
disp_rel.plot(range(len(GX)+2*len(XW),len(GX)+2*len(XW)+len(KG)), [10**-12*x for x in list(reversed(freq_KG))],color, ms=markersize)
disp_rel.plot(range(len(GX)+2*len(XW)+len(KG), len(GX)+2*len(XW)+len(KG)+len(GL)), [10**-12*x for x in freq_GL],color, ms=markersize)


Label=[r'$\Gamma$','X','W','X','K',r'$\Gamma$', 'L']
disp_rel.set_xticks([0,len(GX), len(GX)+len(XW)-0.5,len(GX)+2*len(XW)-0.5
, len(GX)+2*len(XW)+len(KG)/4-0.5,len(GX)+2*len(XW)+len(KG)-0.5,len(GX)+2*len(XW)+len(KG)+len(GL)-0.5])

disp_rel.set_xticklabels(Label, minor=False)
disp_rel.axvline(x=len(GX), color='k')
disp_rel.axvline(x=len(GX)+len(XW)-0.5, color='k', linestyle=':')
disp_rel.axvline(x=len(GX)+2*len(XW)-0.5, color='k')
disp_rel.axvline(x=(len(GX)+2*len(XW)+len(KG)/4)-0.5, color='k', linestyle=':')
disp_rel.axvline(x=len(GX)+2*len(XW)+len(KG)-0.5, color='k')
disp_rel.set_xlabel('Wave Number' , fontsize=fontsize, fontname=FN)
disp_rel.set_ylabel('Frequency' r' $THz$', fontsize=fontsize, fontname=FN)
disp_rel.legend(loc='upper left')

DOS.plot(10**12*data['DOS'],10**-12*data['Freq'],'-k', label= '10K')
plt.setp(DOS.get_yticklabels(), visible=False)
DOS.legend(loc='lower right')

#== Compare with the experiment DOI:https://doi.org/10.1103/PhysRevB.10.364 ===

# fig = plt.figure('Disp. Relation',figsize=(3.5, 4))
# plt.subplots_adjust(left=0.205, bottom=0.15, right=0.98, top=0.98, wspace=0, hspace=0)

# grid = plt.GridSpec(3, 4, hspace=0.05, wspace=0.05)

# disp_rel = fig.add_subplot(grid[0:4, 0:4])

# color='^k'
# markersize=2.5
# fontsize=12
# FN= 'Times New Roman'

# disp_rel.plot(range(0,len(GX)),[10**-12*x for x in freq_GX1],color, label= '10K', ms=markersize)
# disp_rel.plot(range(0,len(GX)),[10**-12*x for x in freq_GX2],color, ms=markersize)
# disp_rel.plot(range(0,len(GX)),[10**-12*x for x in freq_GX3],color, ms=markersize)

# disp_rel.plot(range(len(GX),len(GX)+len(KG)), [10**-12*x for x in list(reversed(freq_KG))],color, ms=markersize)
# disp_rel.plot(range(len(GX)+len(KG), len(GX)+len(KG)+len(GL)), [10**-12*x for x in freq_GL],color, ms=markersize)

# Label=[r'$\Gamma$','X','K',r'$\Gamma$', 'L']
# disp_rel.set_xticks([0,len(GX), len(GX)+len(KG)/4-0.5,len(GX)+len(KG)-0.5,len(GX)+len(KG)+len(GL)-0.5])
# disp_rel.set_xticklabels(Label, minor=False)
# disp_rel.set_yticks([0,0.5, 1, 1.5, 2])

# disp_rel.axvline(x=len(GX), color='k')
# disp_rel.axvline(x=(len(GX)+len(KG)/4)-0.5, color='k', linestyle=':')
# disp_rel.axvline(x=len(GX)+len(KG)-0.5, color='k')
# disp_rel.set_xlabel('Wave Number' , fontsize=fontsize, fontname=FN)
# disp_rel.set_ylabel('Frequency' r' $THz$', fontsize=fontsize, fontname=FN)

# k00 = [i*21 for i in [0 , 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9]]
# f00 =  [i*0.24180 for i in [0, 0.9, 1.79, 2.5, 2.67, 3.65, 3.46, 4.79, 4.18, 5.85, 4.76, 6.84, 5.25, 7.44, 
#                             5.61, 7.97, 5.82, 8.27]]
# disp_rel.plot(k00,f00,'ro',  label= 'exp_10K',fillstyle='none')

# k10 = [(abs(i-1)+1)*20 for i in [0 , 0.1, 0.1, 0.1,
#                                  0.2, 0.2,0.2,
#                                  0.3, 0.3, 0.3,
#                                  0.4, 0.4, 0.4,
#                                  0.5, 0.5, 0.5, 
#                                  0.6, 0.6,
#                                  0.7, 0.7, 0.7,
#                                  0.8, 0.8, 0.8,
#                                  0.9, 0.9, 0.9]]
# f10 =  [i*0.24180 for i in [0, 0.85, 1.312, 2.02,
#                             1.65, 2.56, 3.87,
#                             2.46, 3.78, 5.47,
#                             3.24, 4.89, 6.51,
#                             3.95, 5.91, 7.48,
#                             4.63, 7.44,
#                             5.24, 7.45, 7.22,
#                             5.63, 8.03, 6.55,
#                             5.79, 8.27, 6.15]]
# disp_rel.plot(k10, f10,'ro', fillstyle='none')

# k11 = [(i+2)*20 for i in [0, 0.2, 0.2 , 0.4, 0.4, 0.6, 0.6, 0.8 , 0.8,  1, 1]]
# f11 =  [i*0.24180 for i in [0, 1.19, 2.54, 2.31, 5.05, 3.2, 6.83, 3.81, 8.07, 3.97, 8.35]]
# disp_rel.plot(k11, f11,'ro', fillstyle='none')
# disp_rel.legend(loc='upper left',  prop={'size': 5.5})

## ============================================================================
# 2e10 is the frequency interval and has to be implemented in E equation
# T = 10
# hbar = (6.626*10**-34)/(2*np.pi)
# kb = 1.380649*10**-23
# Bose_Einstein = 1/(np.exp(hbar*(2*np.pi*Freq)/(kb*T)))
# E = hbar*2*np.pi*Freq*Bose_Einstein
# data.insert(2,'Energy',E*data['DOS'])


# x = np.array(data['Freq'])
# y = np.array(data['Energy'])
# integ = integrate.simps(y,x)
# print('DOSresult: %s \n3/2kbT: %s' %(integ,3/2*kb*T) )
#fontsize=12
#FN= 'Times New Roman'
#plt.plot(temp, KbT,'.', label ='kbT results')
#plt.plot(temp,dos_res,'o', label='Dos reluts')
#plt.legend(loc='bottom right')
#plt.xlabel('Temperature',fontsize=fontsize, fontname=FN)
#plt.ylabel('Energy' r'  $(Jouls)$',fontsize=fontsize, fontname=FN)
#plt.xticks(temp)

## ============================================================================

#np.savetxt('energy_per_freqinterval_%sK.txt' %T, data,delimiter=" ", header='Freq Energy(Jouls)', fmt='%g %g %g', comments='')
#
#names = ['GX', 'KG' , 'GL']
#for i in names:
#    for j in range(3):
#        globals()['freq_%s%s' %(i,j+1)] = pd.Series(globals()['freq_%s%s' %(i,j+1)])
#Fdata = pd.concat([freq_GX3, freq_GX2, freq_GX1, freq_KG3, freq_KG2, freq_KG1, freq_GL3, freq_GL2, freq_GL1], axis=1)
#Fdata.columns = ['freq_GX1', 'freq_GX2', 'freq_GX3', 'freq_KG1', 'freq_KG2', 'freq_KG3', 'freq_GL1', 'freq_GL2', 'freq_GL3']
#head = 'freq_GX1 freq_GX2 freq_GX3 freq_KG1 freq_KG2 freq_KG3 freq_GL1 freq_GL2 freq_GL3'
#np.savetxt('Dispersion_Frequencies_%sK.txt' %T, np.real(Fdata), delimiter=" ", header= head, comments='')