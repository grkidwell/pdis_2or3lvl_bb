import math
import numpy as np
import pandas as pd

from flycap_cond_states import flycap_cond_states 

cap_filename = r'data/capacitor_data.xlsx'

def stripint(a):
    if type(a) == str:
        a = int(a[0])
    return a

class Caplosses:
    def __init__(self,lout_obj): #ckt_params:dict,idc,ipp,caps:dict):
        self.lout_obj = lout_obj; self.ckt_params = self.lout_obj.ckt
        self.ip = self.ckt_params['ip']; self.caps = self.ip['caps']

        # self.state_count = self.ckt_params['state count']
        # self.ts = {4:(2*self.ckt_params['t_state13']+2*self.ckt_params['t_state24']),
        #            2:self.ckt_params['t_state13']+self.ckt_params['t_state24']}[self.state_count]

        self.fs=self.lout_obj.fs_dcm  
        
        #self.idc=idc; self.ipp=ipp
        

        self.add_esr_to_capdict(cap_filename)
        self.summary = {'vin':self.p_inputcap(),
                        'flying':self.p_flyingcap(),
                        'vout':self.p_outputcap()}

    def add_esr_to_capdict(self,filename):
        df_caps = pd.read_excel(filename)
        esr_dict = df_caps.set_index('description').loc['ESR'][3:].to_dict()
        for cap,value in self.caps.items():
            self.caps[cap]['esr']=esr_dict[value['partnum']]
        

    def p_flyingcap(self):
        def rms_dcm_calc():
            ts = 1/(2*self.fs) #cap current pulses 2x per cycle to balance charge
            i_scaler = {'single':1,'series':1,'parallel':2}[self.ip['lout']['config']]
            cond_states = flycap_cond_states(self.ckt_params)
            ims = self.lout_obj.i_ms_bystate
            tstate = self.ckt_params['t_state']
            ms_tot = sum([0+tstate[stripint(state)]*ims[stripint(state)] for state in cond_states if state !=0])
            return (ms_tot/ts)**0.5*i_scaler
        cap = self.caps['flying']
        esr_equiv = cap['esr']/cap['n']
        # irms_4state = ((self.idc**2+self.ipp**2/2)*2*self.ckt_params['t_state24']/self.ts)**0.5
        # irms_2state = 0
        # irms = {4:irms_4state,
        #         2:irms_2state}[self.ckt_params['state count']]
        irms=rms_dcm_calc()
        return irms**2*esr_equiv
    
    def p_inputcap(self):
        cap = self.caps['vin']
        esr_equiv = cap['esr']/cap['n']
        #d = self.ckt_params['t_Qhs']/self.ts

        #need to re-derive this for both 4state and 2state and for dcm
        irms_4state = 0 #self.idc*(d*(1-d)+d*(1-d)**2*self.ipp**2/12)**0.5
        
        return irms_4state**2*esr_equiv

    def p_outputcap(self):
        cap = self.caps['vout']
        esr_equiv = cap['esr']/cap['n']
        irms = 0# self.ipp*(1/12)**0.5
        return irms**2*esr_equiv
        
