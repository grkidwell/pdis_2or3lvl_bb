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
        return complex(resistance,reactance)/capdict['n']

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
        ff=np.logspace(2, 7, 1000)
        Zeq_mag_vect = Z_mag_vect(self.Z_eq_allbanks,ff)
        for bankidx,capparamdict in enumerate(self.listofbanks):
            def Zfunc(f):
                return self.Zbank(capparamdict,f)
            plt.loglog(ff,Z_mag_vect(Zfunc,ff),label=f'bank {bankidx}')
        plt.loglog(ff,Zeq_mag_vect,label='Zeq allbanks')#,color='magenta')
        plt.legend(); plt.grid(True); plt.grid(True,which="minor",color='lime',linewidth=.2)