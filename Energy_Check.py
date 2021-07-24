# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:39:18 2019

@author: A02197412
"""
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np


for T in np.arange(10,90,20):
#T = 70
    std_data = []
    mean_data = []
    Branch1 = 'freq_GX1'
    Branch2 = 'freq_GX3'
    
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
    
    # =============================================================================
    
    f=open('F:\Argon\Create the Signal\Phase 1\Amplitude\DATA\energy_per_freqinterval_%sK.txt' %T)
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
       return print('The Frequency is not in range')  
    
    def interpolEng(freq,Branch):
       if freq<=max(energy['Freq']):
           ENERGIES = []
           kpos = interpolk(freq,Branch)+ 3/(sig)
           kneg = interpolk(freq,Branch)- 3/(sig)
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
    frequencies = [2e11,4e11,6e11,8e11,1e12]
    CalcEng = [1000*ENG(frequencies[0],Branch1,Branch2, 100)*960*6.242e+18, 1000*ENG(frequencies[1],Branch1,Branch2, 100)*960*6.242e+18,
               1000*ENG(frequencies[2],Branch1,Branch2, 100)*960*6.242e+18, 1000*ENG(frequencies[3],Branch1,Branch2, 100)*960*6.242e+18,
               1000*ENG(frequencies[4],Branch1,Branch2, 100)*960*6.242e+18]
    
    formatc='hk'
    if T==10:
        formatshade = 'lightblue'
    elif T==30:
        formatshade = 'lightsteelblue'
    elif T==50:
        formatshade = 'cornflowerblue'
    elif T==70:
        formatshade = 'royalblue'
    
    for f in frequencies:
        FREQUENCY = f
        freq = '%.2e' %FREQUENCY
        f=open('%sK\%s\log.lammps' %(T,freq))
        x=f.read()
        x = x.split('\n')
        index = x.index('Step Temp KinEng PotEng TotEng ')
        columns = x[index].split(' ')[0:5]
        x = x[index+1:index+52]
        x = [x[a].split() for a in range(len(x))]
        x = [[float(x[a][b]) for b in range(len(x[a]))] for a in range(len(x))]
        data = pd.DataFrame(x)
        data.columns = columns
        del(x, columns,freq,FREQUENCY)
        data = data.iloc[:,2]
        std_data.append(np.std(data))
        mean_data.append(np.mean(data))
    
    
    
    if T==10:
        fontsize=12
        FN= 'Times New Roman'
        fig = plt.figure('Energy Oscillation of 10_30_50_70K' ,figsize=(6, 7))
        plt.subplots_adjust(left=0.10, bottom=0.12, right=0.98, top=0.98, wspace=0, hspace=0)
        grid = plt.GridSpec(2, 2, hspace=0.15, wspace=0.11)
        ax  = fig.add_subplot(grid[0:2, 0:2])
        ax1 = fig.add_subplot(grid[0, 0])
        ax2 = fig.add_subplot(grid[0, 1])
        ax3 = fig.add_subplot(grid[1, 0])
        ax4 = fig.add_subplot(grid[1, 1])
        ax1.errorbar([x*10**-12 for x in frequencies], CalcEng, fmt=formatc,markersize=5, capsize=2, label='Expected Energy at %sK' %T)
        ax1.fill_between([x*10**-12 for x in frequencies], [1000*(m+s) for m, s in zip(mean_data,std_data)],
                         [1000*(m-s) for m, s in zip(mean_data,std_data)],color= formatshade, label='$\sigma$ KE at %sK MD' %T)
        ax1.legend(loc ='upper left', prop={'size': 7})
        ax1.set_title('(a)', y=-0.16)
        
    if T==30:
        ax2.errorbar([x*10**-12 for x in frequencies], CalcEng, fmt=formatc,markersize=5, capsize=2, label='Expected Energy at %sK' %T)
        ax2.fill_between([x*10**-12 for x in frequencies], [1000*(m+s) for m, s in zip(mean_data,std_data)],
                         [1000*(m-s) for m, s in zip(mean_data,std_data)],color= formatshade, label='$\sigma$ KE at %sK MD' %T)
        ax2.legend(loc ='upper left', prop={'size': 7})
        ax2.set_title('(b)', y=-0.16)
        
    elif T==50:
        ax3.errorbar([x*10**-12 for x in frequencies], CalcEng, fmt=formatc,markersize=5, capsize=2, label='Expected Energy at %sK' %T)
        ax3.fill_between([x*10**-12 for x in frequencies], [1000*(m+s) for m, s in zip(mean_data,std_data)],
                         [1000*(m-s) for m, s in zip(mean_data,std_data)],color= formatshade, label='$\sigma$ KE at %sK MD' %T)
        ax3.legend(loc ='upper left', prop={'size': 7})
        ax3.set_title('(c)', y=-0.16)
        
    elif T==70:
        ax4.errorbar([x*10**-12 for x in frequencies], CalcEng, fmt=formatc,markersize=5, capsize=2, label='Expected Energy at %sK' %T)
        ax4.fill_between([x*10**-12 for x in frequencies], [1000*(m+s) for m, s in zip(mean_data,std_data)],
                         [1000*(m-s) for m, s in zip(mean_data,std_data)],color= formatshade, label='$\sigma$ KE at %sK MD' %T)
        ax4.legend(loc ='upper left', prop={'size': 7})
        ax4.set_title('(d)', y=-0.16)
    
#    plt.errorbar([x*10**-12 for x in frequencies], CalcEng, fmt=formatc,markersize=5, capsize=2, label='Expected Energy at %sK' %T)
#    plt.fill_between([x*10**-12 for x in frequencies], [1000*(m+s) for m, s in zip(mean_data,std_data)], [1000*(m-s) for m, s in zip(mean_data,std_data)],
#                                   color= formatshade, label='$\sigma$ KE at %sK' %T)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax.set_xlabel('Frequency ' r'$(THz)$',fontsize=fontsize, fontname=FN)
ax.set_ylabel('Energy ' r'$(mev)$',fontsize=fontsize, fontname=FN)
#    plt.ylabel('Energy ' r'$(mev)$',fontsize=fontsize, fontname=FN)
    

