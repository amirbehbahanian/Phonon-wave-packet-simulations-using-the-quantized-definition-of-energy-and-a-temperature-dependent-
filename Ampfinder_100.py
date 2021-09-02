# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:12:04 2019

@author: A02197412
"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import sympy as sy

#eta1=0.1
#energy=2.2e-25

def Ampfinder(energy,eta1):

    num_unitcell_iny = 60
    lam=5.2562*10**-10
    eta2=eta1
    k1=eta1*np.pi/lam
    k2=eta2*np.pi/lam    
    sig=52.562*10**-10
    mean=num_unitcell_iny*lam/2
    expden=2*sig**2
    
    y = np.arange(0,60*lam,lam)
    z = np.zeros(len(y))
    x = np.zeros(len(y))
    
    df = np.column_stack((x,y,z))
    df = pd.DataFrame(df)
    df.columns = ['x','y','z']
    
    
    # ===============================Signal shape==================================
    """
    Caculation of the Wave amplitude based on the Kinetic energy at each temperture
    """
    ep=1.6699999994802e-21#0.0104233206
    sig=3.4
    r=sy.Symbol('r')
    U=4*ep*((sig/r)**12-(sig/r)**6)
    dif=sy.diff(U,r)
    ans=sy.solve(dif,r)
    ans=[complex(x) for x in ans]
    ans=np.array(ans)
    for i in range(len(ans)):
        if (ans[i].imag==0 and ans[i].real>0):
            min_r=ans[i]
            break
    del(ans)
    min_r=min_r.real
    min_U=sy.limit(U,r,min_r)
    #                        ==========================
    KE= energy
    #                        ==========================
    ans=sy.solve(U-(min_U+KE),r)
    ans=[complex(x) for x in ans]
    ans=np.array(ans)
    Result=[]
    for i in range(len(ans)):
        if (ans[i].imag==0 and ans[i].real>0):
            Result.append(ans[i])
    #del(ans)
    Result=np.real(Result)
    if len(Result)!=2:
        atom_disp = 'Wrong energy'
    else:
        #                       ===Amplitude Calculation===
        osc_amp = max(Result)-min(Result)
        #del(Result,ep,i,min_r,sig,min_U,r,U)
        """
                            End of Amplitude Calculation
        """
        
        atom_disp=0.000
        control = 0
        while control<osc_amp:
            dft = df
            atom_disp = atom_disp+0.005
            def delz(pos):
                delz1 = (atom_disp)*np.cos(pos*k2)*np.exp(-(pos-mean)**2/expden)
                return delz1
        # =============================================================================
        # 
        # =============================================================================
            dft=dft.loc[(dft['x']==0) & (dft['z']==0)]
            dft.sort_values('y', axis = 0, ascending = True, inplace = True, na_position ='last')
            for i in range(len(dft)):
                dft['z'][i]= dft['z'][i]+delz(dft['y'][i])
            dffinal=dft
            dffinal=pd.DataFrame(dffinal)
            
            dffinal.sort_values('y', axis = 0, ascending = True, inplace = True, na_position ='last') 
            
            
            dffinal['y']=[lam*x for x in range(len(dffinal))]
            
            dffinal.insert(3,'Dx',abs(dffinal.iloc[:,0].diff(periods=1)))
            dffinal.insert(4,'Dy',abs(dffinal.iloc[:,1].diff(periods=1)))
            dffinal.insert(5,'Dz',abs(dffinal.iloc[:,2].diff(periods=1)))
            
            dffinal.insert(6,'At-At Distance',np.sqrt(np.square(dffinal.iloc[:,3:6]).sum(axis=1)))   
            
            dffinal['Deviation from Eq']=abs(dffinal['At-At Distance'])-lam
            control = dffinal.iloc[1:-1,7].mean()
    return atom_disp

#plt.figure('Wave')
#plt.plot(dffinal['y'],dffinal['z'],'.k', dffinal.iloc[1:-1,1],dffinal.iloc[1:-1,7],'.r')