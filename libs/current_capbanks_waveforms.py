import libs.append_path
from add_python_libraries import *
import numpy as np
import scipy as sp
import scipy.signal

class Current_capbanks_waveforms:
    def __init__(self,currents_capsfftobj):
        self.i_capsfftobj = currents_capsfftobj
        self.i_total_ac = sp.fft.ifft(self.i_capsfftobj.i_buckfftobj.i_zin_fft).real
        self.i_nbanks_ac = [self.do_ifft(bankdict) for bankdict in self.i_capsfftobj.i_nbank_fft]
        self.i_rms_per_cap = self.calculate_rms_currents()  #dict

    def do_ifft(self,bankdict):
        return sp.fft.ifft(bankdict['fftdata']).real

    def calculate_rms_currents(self):
        i_cap_dict = {}
        for idx,ibank in enumerate(self.i_nbanks_ac):
            bank_dict = self.i_capsfftobj.i_nbank_fft[idx]
            i_cap_dict[bank_dict['label']] = np.round(numpy_rms.rms(ibank)[0]/bank_dict['n']*2**bank_dict['series'],1)
        return i_cap_dict

    def plot_ibanks(self,ax):
        fstep = self.i_capsfftobj.freq[1]
        n = len(self.i_capsfftobj.freq)
        dt = 1/(n*fstep)
        time_axis = np.arange(0,n*dt,dt)
        for idx,ibank in enumerate(self.i_nbanks_ac):
            bank_dict = self.i_capsfftobj.i_nbank_fft[idx]
            irms_per_cap = self.i_rms_per_cap[bank_dict['label']]
            label = f'{bank_dict['label']} - {irms_per_cap} Arms/cap'
            ax.plot(time_axis,ibank,label=label)
        ax.set_xlim(0)
        ax.legend(loc='upper left',bbox_to_anchor=(0.0,-0.12))
        ax.grid(True,color='lime',linewidth=.2)
