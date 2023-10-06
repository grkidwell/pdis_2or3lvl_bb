import math
import numpy as np
import pandas as pd

cap_filename = r'data/capacitor_data.xlsx'



class Caplosses:
    def __init__(self,ckt_params:dict,idc,ipp,caps:dict):
        self.ckt_params = ckt_params
        self.state_count = self.ckt_params['state count']
        self.ts = {4:(2*self.ckt_params['t_state13']+2*self.ckt_params['t_state24']),
                   2:self.ckt_params['t_state13']+self.ckt_params['t_state24']}[self.state_count]
        self.idc=idc; self.ipp=ipp
        self.caps = caps

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
        cap = self.caps['flying']
        esr_equiv = cap['esr']/cap['n']
        irms_4state = ((self.idc**2+self.ipp**2/2)*2*self.ckt_params['t_state24']/self.ts)**0.5
        irms_2state = 0
        irms = {4:irms_4state,
                2:irms_2state}[self.ckt_params['state count']]
        return irms**2*esr_equiv
    
    def p_inputcap(self):
        cap = self.caps['vin']
        esr_equiv = cap['esr']/cap['n']
        d = self.ckt_params['t_Qhs']/self.ts
        irms_4state = self.idc*(d*(1-d)+d*(1-d)**2*self.ipp**2/12)**0.5
        return irms_4state**2*esr_equiv

    def p_outputcap(self):
        cap = self.caps['vout']
        esr_equiv = cap['esr']/cap['n']
        irms = self.ipp*(1/12)**0.5
        return irms**2*esr_equiv
        
