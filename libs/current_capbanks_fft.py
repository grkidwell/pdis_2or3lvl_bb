import numpy as np
import scipy as sp
import scipy.signal

import libs.append_path
from add_python_libraries import *

class Current_capbanks_FFT:
    def __init__(self,currentsbuckfftobj,zcapbanksobj):
        self.i_buckfftobj = currentsbuckfftobj
        self.z_capbankobj = zcapbanksobj
        self.freq = self.i_buckfftobj.freq; self.freq[0] = 1
        self.Fs = self.freq[1]
        self.zinptotal = np.vectorize(self.z_capbankobj.Z_eq_allbanks)(self.freq)
        self.i_nbank_fft = [self.add_fft_data(bankdict) for bankdict in self.z_capbankobj.listofbanks]

    def Zbank_vectorize(self,capdict,f):
        def zbank_f(f):
            return self.z_capbankobj.Zbank(capdict,f)
        return np.vectorize(zbank_f)(f)

    def i_bank_fft(self,capdict):
        return self.i_buckfftobj.i_zin_fft*self.zinptotal/self.Zbank_vectorize(capdict,self.freq)

    def add_fft_data(self,bankdict):
        bankdict['fftdata'] = self.i_bank_fft(bankdict)
        return bankdict
    
    def plot_i_zbank_fft(self,icapfft):
        i_bank_fft_mag = np.abs(icapfft)
        ampscale=2/icapfft.size
        xticks = np.arange(8)*self.Fs
        i=self.freq>=0
        bar_width=300000
        fig,ax = plt.subplots(1,1,figsize=(8,4))
        rects1=ax.bar(self.freq[i],ampscale*i_bank_fft_mag[i],bar_width)
        ax.set_xlim(0,8*self.Fs)
        ax.set_xticks(xticks)
        plt.show()

    def plot_i_allbanks_fft(self):
        i_nbank_fft = self.i_nbank_fft
        f, ax = plt.subplots(2, 2, sharey=True)#,figsize=(12,4))
        freq = self.freq
        Fs = round(self.Fs,-4) #round(freq[1],-4)
        i_nbank_fft_mag = [np.abs(capdict['fftdata']) for capdict in i_nbank_fft]
        ampscale = 2/i_nbank_fft_mag[0].size
        i=freq>=0
        xticks = np.arange(8)*Fs
        bar_width=300000
        for n, ax_n in enumerate(ax.flat): 
            if (n < len(i_nbank_fft_mag)):
                rects = ax_n.bar(freq[i],ampscale*i_nbank_fft_mag[n][i],bar_width)
                ax_n.set_title(i_nbank_fft[n]['label'])
                ax_n.set_xlim(0,8*Fs)
                ax_n.set_xticks(xticks)#,labels=xlabels)
        f.tight_layout()
        plt.show()