import numpy as np
import scipy as sp
import scipy.signal

class Current_capbanks_waveforms:
    def __init__(self,currents_capsfftobj):
        self.i_capsfftobj = currents_capsfftobj
        self.i_total_ac = sp.fft.ifft(self.i_capsfftobj.i_buckfftobj.i_zin_fft).real
        self.i_nbanks_ac = [self.do_ifft(bankdict) for bankdict in self.i_capsfftobj.i_nbank_fft]

    def do_ifft(self,bankdict):
        return sp.fft.ifft(bankdict['fftdata']).real