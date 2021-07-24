# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:55:59 2019

@author: A02197412
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

test = 'off'
FREQUENCY = 8e11
T = 70
Branch = 'freq_GX1'
step = 10
fftstep = 253  
### for 50 and 30: 20, 251|251|253|251
# =============================================================================

f=open('F:\Argon\Create the Signal\Phase 1\Amplitude\DATA\Dispersion_Frequencies_%sK.txt' %T)
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
   return print('The Frequency is not in range')  

# =============================================================================


def section(i):
    freq = '%.2e' %FREQUENCY
    df1=[]
    f=open('F:/Argon/Analysis/Phase 1/Main/Single Wave Interaction/k-Energy-Vg_Check/%sK/%s/dump.PE.allsections.%s' %(T,freq, i*200))
    x=f.read()
    x = x.split('\n')
    x = x[9:-1]
    x = [x[a].split() for a in range(len(x))]
    x = [[float(x[a][b]) for b in range(len(x[a]))] for a in range(len(x))]
    df1=pd.DataFrame(x)
    del(x)
    df1.columns=['ID', 'c_PE', 'x', 'y', 'z', 'vx', 'vy', 'vz']
    df1=df1.reset_index(drop=True)
    df1=df1.drop([len(df1)-1])
    df1=df1.sort_values(by=['y'])
    
    F = np.fft.fft(df1.loc[(df1['x']>1.5)& (df1['x']<4) & (df1['z']>1.5) & (df1['z']<4)]['vz'])
    F = abs(np.real(F))
    return F

ffk_index = np.where(section(step)==max(section(step)))[0][0]

ffk = []
for i in range(0, fftstep):
    ffk.append( section(i)[ffk_index])
    
# =============================================================================
# ##ffk = pd.DataFrame(np.zeros(fftstep),columns=['ffk(k,t)'])
# ##for i in range(0, fftstep):
# ##    ffk.loc[i] = section(i)[ffk_index]
# =============================================================================

def c(delt):
    c = 0
    for i in range(len(ffk)):
        try:
            c = c + ffk[i]*ffk[i+delt]
        except IndexError:
            break
    c = c/len(ffk)
    return c

c0=0
for i in range(len(ffk)):
    c0 = c0 + ffk[i]**2
c0 = c0/len(ffk)


R = [abs(c(i)) for i in range(len(ffk))]
R = [x-np.mean(R) for x in R]

# =============================================================================
# ##ffk.insert(1,'ffk(k,0).ffk(k,t)',0)
# ##ffk0 = ffk.iloc[1,0]
# ##for i in range(len(ffk)):
# ##    ffk.iloc[i,1] = ffk0 * ffk['ffk(k,t)'].loc[i]
# ##
# ##ffk0 = ffk.iloc[1,1]    
# ##ffk.insert(2,'<ffk(k,0).ffk(k,t)>',0)
# ##for i in range(len(ffk)):
# ##    if i==0:
# ##        ffk.iloc[i,2] = 0
# ##    else:
# ##        ffk.iloc[i,2] = ffk['ffk(k,0).ffk(k,t)'].loc[0:i].sum()/((i)*ffk0)
# =============================================================================

FFT = np.fft.fft(R)
FFT = np.real(FFT)
FFT = abs(FFT)
Fs = 1/(200*10**-15)
freq = np.arange(0,(Fs/2),Fs/len(FFT))
fontsize=12
FN= 'Times New Roman'

# =============================================================================
if test=='off':
    if FREQUENCY==2e11:
        fig = plt.figure('K and W FFT at %sK' %T,figsize=(4,4))
        grid = plt.GridSpec(4, 4, hspace=0.05, wspace=0.05)
        plt.subplots_adjust(left=0.18, bottom=0, right=0.98, top=0.98, wspace=0, hspace=0)
        kfft = fig.add_subplot(grid[:-1, 0:2])
        wfft = fig.add_subplot(grid[:-1, 2:], sharey=kfft)
if test=='on':
    fig = plt.figure('K and W FFT at %sK' %T,figsize=(7.2, 6))
    grid = plt.GridSpec(4, 4, hspace=0.05, wspace=0.15)
    plt.subplots_adjust(left=0.09, bottom=0, right=0.98, top=0.98, wspace=0, hspace=0)
    kfft = fig.add_subplot(grid[:-1, 0:2])
    wfft = fig.add_subplot(grid[:-1, 2:], sharey=kfft)


### Lattice Constant
if T==10:
    lam = 5.28936  
    zfreq = [0.8e11, 3.4e11, 5.836e11, 8.801e11]       
elif T==30:
    lam = 5.33567
    zfreq = [1.698e11, 4.203e11, 6.603e11, 9.489e11]
elif T==50:
    lam = 5.39197
    zfreq = [1.713e11, 4.229e11, 7.364e11, 10.168e11]
elif T==70:
    lam = 5.46562 
    zfreq = [1.728e11, 5.040e11, 8.797e11, 11.933e11]
        
dist=lam/2
sig1= lam*10        ### Standard Deviation of the Gaussian Envelope
k = np.arange(0,0.5,0.001)
sig2 = 1/(sig1)
F = section(step)
K=interpolk(FREQUENCY,Branch)

if FREQUENCY==2e11:
    W=0.2
    MARKER1='^y'
    MARKER2='^-y'
    wfft.arrow(zfreq[0]/(10**12),0,0,-0.2,length_includes_head=True, facecolor='yellow', edgecolor='black',
          head_width=0.08, head_length=0.04)
    wfft.axvline(x=W, color='k',linestyle=':',label=r'Exp.Freq.')
if FREQUENCY==4e11:
    W=0.4
    MARKER1='vg'
    MARKER2='v-g'
    wfft.arrow(zfreq[1]/(10**12),0,0,-0.2,length_includes_head=True, facecolor='green', edgecolor='black',
          head_width=0.08, head_length=0.04)
    wfft.axvline(x=W, color='k',linestyle=':')
elif FREQUENCY==6e11:
    W=0.6
    MARKER1='^c'
    MARKER2='^-c'
    wfft.arrow(zfreq[2]/(10**12),0,0,-0.2,length_includes_head=True, facecolor='cyan', edgecolor='black',
          head_width=0.08, head_length=0.04)
    wfft.axvline(x=W, color='k',linestyle=':')
elif FREQUENCY==8e11:
    W=0.8
    MARKER1='<r'
    MARKER2='<-r'
    wfft.arrow(zfreq[3]/(10**12),0,0,-0.2,length_includes_head=True, facecolor='red', edgecolor='black',
          head_width=0.08, head_length=0.04)
    wfft.axvline(x=W, color='k',linestyle=':')

    
Fg = np.exp(-((k-K/2)**2)/(2*sig2**2))

fontsize=12
FN= 'Times New Roman'

wfft.plot(freq[1:]/(10**12),FFT[1:len(freq)]/max(FFT[1:len(freq)]),MARKER2, label= '%s THz'%(FREQUENCY/(10**12)) )
wfft.set_xlabel('Frequency (THz)',fontsize=fontsize, fontname=FN)
wfft.set_xticks(np.arange(0,2.8,0.6))
plt.setp(wfft.get_yticklabels(), visible=False)
wfft.legend(loc ='upper right', prop={'size': 7})
wfft.set_ylim(-0.2,1.05)


kfft.plot(np.arange(1.0/int(len(F)),0.5,1.0/int(len(F))), F[1:int(len(F)/2)]/max(F[1:int(len(F)/2)]),MARKER1
          , markersize = 4, label ='%s THz' %(FREQUENCY/(10.0**12)))
kfft.plot(k,Fg,'k')
kfft.set_xlabel('Wave number  ' r'(units of $\frac{2\pi}{a}$)',fontsize=fontsize, fontname=FN)
kfft.set_ylabel('Normalized FFT value',fontsize=fontsize, fontname=FN)
kfft.set_xticks(np.arange(0,0.6,0.1),(0,0.2,0.4,0.6,0.8,1))
kfft.set_xlim(0,0.5)
kfft.set_ylim(-0.2,1.05)
kfft.legend(loc ='upper right', prop={'size': 6.5})


# =============================================================================
# fig=plt.figure()
# ax=fig.add_subplot(111,projection='3d')
# ax.scatter(df1['x'],df1['y'],df1['z'],c='black')
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z') 
# =============================================================================
