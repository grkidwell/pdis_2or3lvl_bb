import numpy as np
import scipy as sp
import scipy.signal

class Current_buck_waveforms:
    def __init__(self, ind_obj):
        self.lobj = ind_obj
        self.level = self.lobj.ckt['ip']['lvl_config'] 
        self.fs = self.lobj.ckt['ip']['fs']
        self.ts13 = self.lobj.ckt['t_state13']; self.ts24 = self.lobj.ckt['t_state24']
        self.duty = self.ts13/(self.ts13+self.ts24)        
        self.dfactor = {'2 level':1,'3 level':2}[self.level]
        self.cycles=2; self.numpoints=1024; self.datapoints=np.linspace(0,1,self.numpoints) 
        self.time_datapoints=np.linspace(0,self.cycles/(self.fs*self.dfactor),self.numpoints)
        self.i_ind = self.create_ind_waveform()
        self.i_hs = self.create_hs_waveform()
        self.i_zin_ac  = self.create_i_zin_waveform()
        
    def create_ind_waveform(self):
        ipp = self.lobj.ipp; iout = self.lobj.idc      
        return sp.signal.sawtooth(2*np.pi*self.datapoints*self.cycles,self.duty)/2*ipp+iout

    def create_hs_waveform(self):
        d_hs = self.lobj.ckt['t_Qhs']/(self.ts13+self.ts24)/self.dfactor
        D_pulse = sp.signal.square(2*np.pi*self.datapoints*self.cycles/self.dfactor,d_hs)/2+.5
        return self.i_ind*D_pulse
        
    def create_i_zin_waveform(self):
        pin = self.lobj.ckt['ip']['pin']; vin = self.lobj.ckt['ip']['vin']
        Iin_avg = np.mean(self.i_hs) #pin/vin
        return self.i_hs-Iin_avg