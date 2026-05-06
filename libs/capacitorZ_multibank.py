import libs.append_path
from add_python_libraries import *

import cmath
import numpy as np
import matplotlib.pyplot as plt

class Zcapbanks:  #single frequency functions
    def __init__(self,caplist:list):  #list of dictionaries
        self.listofbanks = caplist
        
    def Zbank(self,capdict,f):  #capdict keys = label,type,c,esr,esl,n
        resistance = capdict['ESR']
        reactance = capdict['ESL']*2*np.pi*f -1/(2*np.pi*f*capdict['C'])
        z_1cap = complex(resistance,reactance)
        if 'series' in capdict:
            series = int(capdict['series'])
        else:
            series = 0
        n = int(capdict['n'])
        if series and n%2:
            print("series config quantities must be multiples of 2")
        zmult = 2**series
        zgroups = n/zmult
        return z_1cap*zmult/zgroups    #complex(resistance,reactance)/df_cap['qnty']


    def Zeq(self,caplist,f):
        return sum([self.Zbank(capparamdict,f)**-1 for capparamdict in caplist])**-1

    def Z_eq_allbanks(self,f):
        return self.Zeq(self.listofbanks,f)

    def plot_Zbanks(self):
        def Z_mag_vect(Zfunc,freqdata):
            def Z_mag(f):
                return abs(Zfunc(f))
            Z_vect =  np.vectorize(Z_mag)
            return Z_vect(freqdata)
        ff=np.logspace(2, 8, 1000)
        Zeq_mag_vect = Z_mag_vect(self.Z_eq_allbanks,ff)
        for bankidx,capparamdict in enumerate(self.listofbanks):
            def Zfunc(f):
                return self.Zbank(capparamdict,f)
            plt.loglog(ff,Z_mag_vect(Zfunc,ff),label=f'bank {bankidx}')
        plt.loglog(ff,Zeq_mag_vect,label='Zeq allbanks')#,color='magenta')
        plt.legend(); plt.grid(True); plt.grid(True,which="minor",color='lime',linewidth=.2)

class Zcapbanks_df:  #single frequency functions
    def __init__(self,dfcaplist):  
        self.df_caplist = dfcaplist
        self.add_bank_label()
        self.listofbanks = self.df_caplist_to_listofbanks()

    #this function is just a patch to make this object compatible with current_capbanks_fft.py, which 
    #expects the dictlist listofbanks.   Need to refactor that library
    def df_caplist_to_listofbanks(self):
        df=self.df_caplist.copy()[['label','type','cap_uF(@vbias)','ESR_mOhm','ESL_nH','qnty','series']]
        df=df.rename(columns={'cap_uF(@vbias)':'C','ESR_mOhm':'ESR','ESL_nH':'ESL','qnty':'n'})
        df.C=df.C*1e-6;df.ESR=df.ESR*1e-3;df.ESL=df.ESL*1e-9
        lob = df.to_dict(orient='records')
        for bank in lob:
            if type(bank['series'])==str:
                bank['series']=0.0
        return lob #df.to_dict(orient='records')
    
    def add_bank_label(self):
        df = self.df_caplist
        for captype in df.type.unique(): #['mlcc','SP']:
            df_type=df[df.type==captype].copy();df_type=df_type.assign(label=range(len(df_type)))
            df_type['label']+=1
            df_type['label']=df_type['type']+df_type['label'].astype(str)
            df.loc[df.type==captype,'label'] = df_type['label']

    def Zbank(self,df_cap,f):  #capdict keys = label,type,c,esr,esl,n
        resistance = df_cap['ESR_mOhm']*1e-3
        reactance = df_cap['ESL_nH']*1e-9*2*np.pi*f -1/(2*np.pi*f*df_cap['cap_uF(@vbias)']*1e-6)
        z_1cap = complex(resistance,reactance)
        if type(df_cap['series'])==str:
            series = 0           
        else:
            series = int(df_cap['series'])            
        n = int(df_cap['qnty'])
        if series and n%2:
            print("series config quantities must be multiples of 2")
        zmult = 2**series
        zgroups = n/zmult
        return z_1cap*zmult/zgroups    #complex(resistance,reactance)/df_cap['qnty']

    def Zeq(self,dfcaplist,f):
        return sum([self.Zbank(dfcap,f)**-1 for _,dfcap in dfcaplist.iterrows()])**-1

    def Z_eq_allbanks(self,f):
        return self.Zeq(self.df_caplist,f)

    def plot_Zbanks(self,ax):
        def Z_mag_vect(Zfunc,freqdata):
            def Z_mag(f):
                return abs(Zfunc(f))
            Z_vect =  np.vectorize(Z_mag)
            return Z_vect(freqdata)
        ff=np.logspace(2, 8, 1000)
        Zeq_mag_vect = Z_mag_vect(self.Z_eq_allbanks,ff)
        for _,dfcap in self.df_caplist.iterrows():
            def Zfunc(f):
                return self.Zbank(dfcap,f)
            ax.loglog(ff,Z_mag_vect(Zfunc,ff),label=f'{dfcap.label} bank')
        ax.loglog(ff,Zeq_mag_vect,label='Zeq allbanks')#,color='magenta')
        ax.legend(loc='upper left',bbox_to_anchor=(0.0,-0.12))
        ax.grid(True); ax.grid(True,which="minor",color='lime',linewidth=.2)