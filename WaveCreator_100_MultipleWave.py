# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 12:20:09 2020

@author: A02197412
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 18:03:01 2019

@author: A02197412
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Ampfinder_100 import Ampfinder

from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
from decimal import Decimal


""" SET THE Constants"""
###################
T = 10
AtomCoord = '%sK\\atomcord_%sK.0' %(T,T)

### Lattice Constant
if T==10:
    lam = 5.28936         
elif T==30:
    lam = 5.33567
elif T==50:
    lam = 5.39197
elif T==70:
    lam = 5.46562 

dist=lam/2
sig= lam*10        ### Standard Deviation of the Gaussian Envelope

###################

""" SET THE FREQUENCIES"""
###################
Branch_TA = 'freq_GX1'
Branch_LA = 'freq_GX3'

fromFREQ = 2
FREQNUM = 10
FREQSTEP = 0.1

counter = 0
for i in np.arange(fromFREQ,FREQNUM,FREQSTEP):
    counter = counter + 1
    globals()['freq%s' %(counter)] = (i)*1e11

Dir = 100
###################

# =============================================================================

f=open('DATA\energy_per_freqinterval_%sK.txt' %T)
x=f.read()
x = x.split('\n')
columns = x[0:1]
columns = columns[0].split()
x = x[1:-1]
x = [x[a].split() for a in range(len(x))]
x = [[float(x[a][b]) for b in range(len(x[a]))] for a in range(len(x))]
energy = pd.DataFrame(x)
energy.columns = columns
del(x, columns)

# =============================================================================

f=open('DATA\Dispersion_Frequencies_%sK.txt' %T)
x=f.read()
x = x.split('\n')
columns = x[0:1]
columns = columns[0].split()
x = x[1:-1]
x = [x[a].split() for a in range(len(x))]
x = [[float(x[a][b]) for b in range(len(x[a]))] for a in range(len(x))]
Freq = pd.DataFrame(x)
Freq.columns = columns
Freq.insert(0,'%K',np.arange(0,1+1/20,1/20))
del(x, columns)

# =============================================================================

"""Get the K value for the required Frequency"""

def interpolk(freq,BranchTA):
   for i in range(len(Freq)):
       if freq<Freq['%s' %(BranchTA)][i] and freq>Freq['%s' %(BranchTA)][i-1]:
           fit = np.polyfit(Freq['%K'],Freq[BranchTA],4)
           P  = np.poly1d(fit)
           roots = np.roots(P-freq)
           k = [abs(x) for x in roots if x>=0 and x<=1]
           k = float(np.array(k))
           return k
#       elif freq==0:
#           return 0
   return print('The Frequency is not in range')  

"""Get the LA Frequency"""    

def interpolLA(freq,BranchTA,BranchLA):
   for i in range(len(Freq)):
       if freq<Freq['%s' %(BranchTA)][i] and freq>Freq['%s' %(BranchTA)][i-1]:
           lowrange = Freq['%s' %(BranchLA)][i-1]
           highrange = Freq['%s' %(BranchLA)][i]
           slope = (highrange - lowrange)/(Freq['%K'][i] - Freq['%K'][i-1])
           LA_freq = slope*(interpolk(freq,BranchTA)- Freq['%K'][i-1])+ Freq['%s' %(BranchLA)][i-1]
           return LA_freq
   return print('The Frequency is not in range') 

"""Get the energy for the required Frequency"""    
       
def interpolEng(freq,Branch):
   if freq<=max(energy['Freq']):
       ENERGIES = []
       kpos = interpolk(freq,Branch)+ 3/(sig)
       if interpolk(freq,Branch)- 3/(sig)>=0:
           kneg = interpolk(freq,Branch)- 3/(sig)
#       elif interpolk(freq,Branch)- 3/(sig)<0:
#           kneg = 0
       P =np.polyfit(Freq['%K'], Freq[Branch],4)
       ENERGIES = energy.loc[(energy['Freq']>np.polyval(P,kneg)) & (energy['Freq']<np.polyval(P,kpos))]
       def G(k):
           sigma = 1/(sig)
           const = (1/((sigma)*np.sqrt(2*np.pi)))
           G = const * np.exp(-(k-(interpolk(freq,Branch)))**2/(2*sigma**2))
           return G
       ENERGIES.insert(3,'k',[(interpolk(x,Branch)) for x in ENERGIES['Freq']])	
       ENERGIES = ENERGIES.reset_index(drop=True)
       ENERGIES.insert(4,'delk',0)
       for i in range(len(ENERGIES)):
           if i<len(ENERGIES)-1:
               ENERGIES.iloc[i,4] = (ENERGIES['k'].loc[i+1] - ENERGIES['k'].loc[i])
           else:
               ENERGIES.iloc[i,4] = (kpos - ENERGIES['k'].loc[i])
       ENERGIES.insert(5,'Gaussian',G(ENERGIES['k']))
       ENERGIES.insert(6,'TotEng',ENERGIES['Energy(Jouls)']*ENERGIES['delk']*ENERGIES['Gaussian'])
       return sum(ENERGIES['TotEng'])
   return print('The Frequency is not in range')
# =============================================================================

""" Branch1 is the Transvese and Branch2 is just the Longinudinal in the same direction
    The following function takes the Ta freq. and LA and TA branches of the same freq
    including the direction of travel and provides the energy in the wave packet"""

# ===================                                         =================
    
def ENG(freq,Branch1,Branch2,LDdirection):
    
    if LDdirection==100:
        LD = 0.33
    elif LDdirection==110:
        LD = 0.47
    else:
        LD = 0.2
    
    P = np.polyfit(Freq['%K'], Freq[Branch2],4)
    k = interpolk(freq,Branch1)
    if np.polyval(P,k)>max(Freq[Branch1]):
        LAcons = 1
    else:
        LAcons = 1/3
    
    total_TA_LA = LD*( ((1/3)*interpolEng(freq,Branch1)) + (LAcons*interpolEng(np.polyval(P,k),Branch2)) )
    return total_TA_LA


# =============================================================================
    ############# The Frequency Entries kick in ########################
counter = 0
for i in np.arange(fromFREQ,FREQNUM,FREQSTEP):
    counter = counter + 1
    globals()['Amplitude%s' %(counter)] = Ampfinder( ENG(globals()['freq%s' %(counter)],Branch_TA,Branch_LA, Dir) , interpolk(globals()['freq%s' %(counter)],Branch_TA))

#### ((interpolEng(freq1)/3)+(interpolEng(freq2)/3))*960*6.242e+18

# =============================================================================

num_unitcell_iny=60
f = open(AtomCoord)
x = f.read()
x = x.split("\n",9)[9]
x = x.replace('\n','/n')
x=' '.join(x.split())
x = x.replace('/n', '\n')
x = x.split('\n')
x = [x[i].split() for i in range(len(x))]
y = [[float(x[i][j]) for j in range(len(x[i]))] for i in range(len(x))]
df=pd.DataFrame(y)
del(x)
del(y)
df=df.drop(columns=0)
df=df.drop([len(df)-1])
df.columns=['x','y','z']

counter = 0
for i in np.arange(fromFREQ,FREQNUM,FREQSTEP):
    counter = counter + 1
    globals()['eta%s' %(counter)] = interpolk(globals()['freq%s' %(counter)],Branch_TA)
    globals()['k%s' %(counter)] = globals()['eta%s' %(counter)]*np.pi/lam


mean=max(df['y'])/2
expden=2*sig**2

def delz(x,k,Amp):
    displace = (Amp)*np.cos((x*k)-mean)*np.exp(-(x-mean)**2/expden)
    return displace

def vel1(x, Amp, freq, k):
    velocity = (-1)*(Amp)*(freq*2*np.pi)*np.sin((x*k)-mean)*np.exp(-(x-mean)**2/expden)
    velocity = velocity*10**-12
    return velocity


dffinal = df  
#disp = []
for i in range(len(df)):
    xlayer = int(round(dffinal.iloc[i,0]/dist))
    ylayer = int(round(dffinal.iloc[i,1]/dist))
    zlayer = int(round(dffinal.iloc[i,2]/dist))
    DISP = 0
    counter = 0
    for j in np.arange(fromFREQ,FREQNUM,FREQSTEP):
        counter = counter + 1
        DISP = DISP + delz(ylayer * dist,globals()['k%s' %(counter)], globals()['Amplitude%s' %(counter)]) 
    dffinal.iloc[i,2] = (zlayer*dist) + DISP
    dffinal.iloc[i,1] = ylayer * dist
    dffinal.iloc[i,0] = xlayer * dist

# =========================Structuring the data file===========================

dffinal.insert(0,'Atom Kind',1)
dffinal.insert(0,'ID',dffinal.index+1)

########################### Velocities ########################################

Velocities = pd.DataFrame(0, index=np.arange(1,len(dffinal)+1), columns = ['vx', 'vy', 'vz'])
for i in range(len(dffinal)):
    Velocities.iloc[i,0] = 0
    Velocities.iloc[i,1] = 0
    Velocities.iloc[i,2] = 0#vel(dffinal.iloc[i,3], Amplitude1, freq1, k1) + vel2(dffinal.iloc[i,3], Amplitude2, freq2, k2)
    
# =========================Check the DATA======================================

#plt.plot(dffinal.loc[(dffinal['x']>1.5)& (dffinal['x']<4) & (dffinal['z']>1.5) & (dffinal['z']<4)]['z'])
    
#fig=plt.figure()
#ax=fig.add_subplot(111,projection='3d')
#ax.scatter(dffinal['x'],dffinal['y'],dffinal['z'],c='black')
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z') 
#ax.set_yticks(np.arange(0,num_unitcell_iny*lam,num_unitcell_iny*lam/3))   

# =========================== Out put file ====================================

Filename = 'atomcord_Freq_%sFreqs_from%.2eto%.2e.txt' %((FREQNUM-fromFREQ)/FREQSTEP,globals()['freq%s' %(1)],globals()['freq%s' %(int((FREQNUM-fromFREQ)/FREQSTEP))])

f = open(AtomCoord)
x = f.read()
x = x.split("\n",8)
head="LAMMPS data file via write_data, version 17 %s, timestep = 0\n\n\
%s atoms\n\
%s atom types\n\n\
%s xlo xhi\n\
%s ylo yhi\n\
%s zlo zhi\n\n\
Masses\n\n\
1 39.948\n\n\
Atoms # atomic\n" %(datetime.today().strftime('%Y-%m-%d'), len(df), max(df['Atom Kind']), x[5], x[6], x[7])
np.savetxt(Filename, df,delimiter=" ", header=head, fmt='%d %d %.6f %.6f %g', comments='')


with open(Filename, 'a') as file:
    file.write('\nVelocities\n\n')
with open(Filename, 'a') as file:
    file.write(Velocities.to_string(header = False))

