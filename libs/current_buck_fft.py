import numpy as np
import scipy as sp
import scipy.signal

import libs.append_path
from add_python_libraries import *

class Current_buck_FFT:
    def __init__(self, current_waveforms_obj):
        self.i_waves_obj = current_waveforms_obj
        self.Fs = self.i_waves_obj.fs
        self.numpoints = self.i_waves_obj.numpoints
        self.timestep = self.i_waves_obj.time_datapoints[1]
        self.i_zin_fft = sp.fft.fft(self.i_waves_obj.i_zin_ac) #[:int(self.numpoints/2)]) if 4 cycles instead of 2
        self.freq = sp.fft.fftfreq(self.i_zin_fft.size,self.timestep)

    def plot_i_zin_fft(self):
        i_zin_fft_mag = np.abs(self.i_zin_fft)
        ampscale=2/self.i_zin_fft.size
        xticks = np.arange(8)*self.Fs
        i=self.freq>=0
        bar_width=300000
        fig,ax = plt.subplots(1,1,figsize=(8,4))
        rects1=ax.bar(self.freq[i],ampscale*i_zin_fft_mag[i],bar_width)
        ax.set_xlim(0,8*self.Fs)
        ax.set_xticks(xticks)
        plt.show()

# class Current_capbanks_FFT:
#     def __init__(self,currentsbuckfftobj,zcapbanksobj):
#         self.i_fftobj = currentsbuckfftobj
#         self.z_capbankobj = zcapbanksobj
#         self.freq = self.i_fftobj.freq; self.freq[0] = 1
#         self.Fs = self.freq[1]
#         self.zinptotal = np.vectorize(self.z_capbankobj.Z_eq_allbanks)(self.freq)

#     def Zbank_vectorize(self,capdict,f):
#         def zbank_f(f):
#             return self.z_capbankobj.Zbank(capdict,f)
#         return np.vectorize(zbank_f)(f)
        
#     def plot_i_zbank_fft(self,icapfft):
#         i_bank_fft_mag = np.abs(icapfft)
#         ampscale=2/icapfft.size
#         xticks = np.arange(8)*self.Fs
#         i=self.freq>=0
#         bar_width=300000
#         fig,ax = plt.subplots(1,1,figsize=(8,4))
#         rects1=ax.bar(self.freq[i],ampscale*i_bank_fft_mag[i],bar_width)
#         ax.set_xlim(0,8*self.Fs)
#         ax.set_xticks(xticks)
#         plt.show()