import libs.append_path
from add_python_libraries import *
import numpy as np
import scipy as sp
import scipy.signal

class Vripple_capbanks_waveforms:
    def __init__(self,i_cbank_wvfm_obj,bw = 10e6):
        self.i_caps_obj = i_cbank_wvfm_obj
        self.v_ripple_fft = self.get_vripple_fft()
        self.bw = bw
        self.volt_waveform = self.get_v_waveform()
        self.vpp = int((max(self.volt_waveform)-min(self.volt_waveform))*1000)
        
    def get_vripple_fft(self):
        capdict = self.i_caps_obj.i_capsfftobj.i_nbank_fft[0]
        i_zbank_fft = capdict['fftdata']
        freq = self.i_caps_obj.i_capsfftobj.freq
        self.sampling_rate = freq[1]
        zbank_fft=self.i_caps_obj.i_capsfftobj.Zbank_vectorize(capdict,freq)
        return i_zbank_fft*zbank_fft

    def get_v_waveform(self):
        def lp_filter(F, band_limit, sampling_rate):
            """Applies a simple lowpass filter using FFT in NumPy."""
            cutoff_index = int(band_limit/ sampling_rate)
            # Zero out frequencies above the cutoff for positive and negative frequencies
            F[cutoff_index + 1:-cutoff_index] = 0 
            return np.fft.ifft(F).real
        return lp_filter(self.v_ripple_fft,self.bw,self.sampling_rate)

    def plot_vripple(self,ax):
        fstep = self.sampling_rate #self.i_capsfftobj.freq[1]
        n = len(self.i_caps_obj.i_capsfftobj.freq)
        dt = 1/(n*fstep)
        time_axis = np.arange(0,n*dt,dt)
        #fig,ax = plt.subplots()
        ax.plot(time_axis,self.volt_waveform,label=f'{self.vpp}mVpp')
        ax.legend();ax.grid(True,color='lime',linewidth=.2)
        #plt.close(fig)
        #return fig, ax