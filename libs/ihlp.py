#in this version, will remove lparams dependent variable from Inductor_pdiss constructor 
#and create lparams dictionary within same constructor
#also change lset function 


import numpy as np
import pandas as pd
from winding_temp_dcm import dcr_temp

from circuit_4state import circuit_params as circuit_params_4state
from circuit_2state import circuit_params as circuit_params_2state

ihlp_filename = 'data/ihlp_core_data.xlsx'

#ckt_params = {'deltaV','duty','fs','Tamb','Idc'}

class Inductor_pdis:
    def __init__(self,inp_params): #ckt_params:dict,Lparams:dict):
        self.ip = inp_params.copy()
#        self.ckt = ckt_params
        self.ckt = {'2 level':circuit_params_2state(self.ip),
              '3 level':circuit_params_4state(self.ip)}\
              [self.ip['lvl_config']]

#        self.ipp = ckt_params['deltaV']*ckt_params['t_for_deltaV']/Lparams['Lout']/1e-6

        l_ip = self.ckt['ip']['lout']
        self.ind = lparams(l_ip['value(uH)'],create_ind_family_df(l_ip['family']))   #Lparams.copy()
        self.idc = self.ckt['Idc']
#        self.ind = Lparams

        self.do_ccm_stuff()
        
        self.ton_mult=self.ckt['ton_mult']
        
        #phasenode dcm and ccm frequency
        self.fs_dcm = round(self.dcm_ratio/(self.ckt['t_state13']+self.ckt['t_state24']),0)
        
        self.irms_dcm = self.i_rms_dcm()

        #self.ind['K1']=0  #no ac winding loss.  see winding_temp.py
        self.p_core    = self.pcore() #ckt_params,Lparams)
        self.tempco = 1/(234.45+25)
        #self.t_winding = round(dcr_temp(ckt_params,Lparams,self.p_core,self.tempco),1)
        self.t_winding = round(dcr_temp(self.irms_dcm,self.ckt,self.ind,self.p_core,self.tempco),1)
        self.DCR = self.ind['DCR']*(1+self.tempco*(self.t_winding-25))
        #self.p_ac = Lparams['K1']*self.ipp**2*ckt_params['fs_phase']**0.5*self.DCR
        self.p_ac = self.ind['K1']*self.ipp**2*self.fs_dcm**0.5*self.DCR
        self.p_dc = self.DCR*self.idc**2
        self.p_tot = self.p_dc+self.p_core+self.p_ac
        self.summary = {'dcr':self.p_dc,
                        'core':self.p_core,
                        'ipp':self.ipp,
                        'fs_dcm':self.fs_dcm,
                        'ton_mult':self.ton_mult,
                        'irms_dcm':self.irms_dcm}

    def do_ccm_stuff(self):
        #ton_mult  = self.ckt['ton_mult']
        ipp       = abs(self.ckt['deltaV']*self.ckt['t_for_deltaV']/self.ind['Lout']/1e-6) 
        dcm_ratio = self.idc/(ipp/2)
        fs_dcm = dcm_ratio/(self.ckt['t_state13']+self.ckt['t_state24'])
        ccm_hyst  = 1.1
        if fs_dcm > 2*self.ip['fs']:
            ton_mult=1
            ip=self.ckt['ip']
            ip['ton_mult']=ton_mult
            lvl_config=ip['lvl_config']
            self.ckt = {'2 level':circuit_params_2state(ip),
                        '3 level':circuit_params_4state(ip)}\
                        [lvl_config]
            ipp = abs(self.ckt['deltaV']*self.ckt['t_for_deltaV']/self.ind['Lout']/1e-6)
            dcm_ratio=1
            
        elif dcm_ratio > ccm_hyst:
            ton_mult = 1
            ip=self.ckt['ip']
            ip['ton_mult']=ton_mult
            lvl_config=ip['lvl_config']
            self.ckt = {'2 level':circuit_params_2state(ip),
                        '3 level':circuit_params_4state(ip)}\
                        [lvl_config]
            ipp = abs(self.ckt['deltaV']*self.ckt['t_for_deltaV']/self.ind['Lout']/1e-6)
            dcm_ratio = self.idc/(ipp/2)

        else:
            dcm_ratio = self.idc/(ipp/2)
            
        self.ipp = ipp
        self.dcm_ratio = min(1,dcm_ratio)

    def pcore(self): #,ckt_params:dict,Lparams:dict):
        #etckt = ckt_params['deltaV']*ckt_params['t_for_deltaV']/1e-6
        #bpk = etckt/Lparams['ET100']*100
        #d = ckt_params['duty']; fs = ckt_params['fs_phase']
        etckt = self.ckt['deltaV']*self.ckt['t_for_deltaV']/1e-6
        bpk = etckt/self.ind['ET100']*100
        d = self.ckt['duty']; fs = self.fs_dcm #ckt['fs_phase']
        fe = fs/(2*np.pi*(d-d**2))
        #return Lparams['K0']*fe**(Lparams['Kf']-1)*bpk**Lparams['Kb']*fs*1e-14
        return self.ind['K0']*fe**(self.ind['Kf']-1)*bpk**self.ind['Kb']*fs*1e-14
     
    def i_rms_dcm(self):
        d1 = self.ckt['t_state13']*self.fs_dcm
        d2 = self.ckt['t_state24']*self.fs_dcm
        return (self.idc**2+(d1+d2)/12*self.ipp**2)**0.5
    
    def losses(self):
        print(f'Total: {round(self.p_tot,3)}')
        print(f'DC: {round(self.p_dc,3)}')
        print(f'AC: {round(self.p_ac,3)}')
        print(f'core: {round(self.p_core,3)}')
        print(f'Temp: {round(self.t_winding,1)}')

def create_ind_family_df(familyname:str):
    return pd.read_excel(ihlp_filename,sheet_name = familyname,engine='openpyxl')
    
def closest_value(input_list, input_value):
    arr = np.asarray(input_list)
    i = (np.abs(arr - input_value)).argmin()
    return arr[i]

def lparams(lout:float,df_ind_family):
    lout_column = df_ind_family['Lout']
    cv = closest_value(lout_column.tolist(),lout)
    df_ind = df_ind_family[lout_column==cv]
    return df_ind.to_dict('records')[0]

#this next function should be unnecessary and depricated
def ind_pdis_obj(ckt_params,lout:float,df_ind_family):
    df_ind = df_ind_family[df_ind_family['Lout']==lout]
    lparams = df_ind.to_dict('records')[0]
    return Inductor_pdis(ckt_params,lparams)

def l_set(inp_params,family:str,inductor_list:list): #Ldict:dict):
    ip = inp_params.copy()
    ip['lout']['family'] = ind_family
    
    df_ind_family = create_ind_family_df(ind_family)
    def ip_replace_ind(ind):
        ip['lout']['value(uH)'] = ind
        return ip
    inductor_set = {ind:Inductor_pdis(ip_replace_ind(ind)) for ind in inductor_list if ind in df_ind_family['Lout'].tolist()}
    pset = [{'Lout':ind_obj.ind['Lout'],
             'Ptot':ind_obj.p_tot,
             'Pdc' :ind_obj.p_dc,
             'Pac' :ind_obj.p_ac,
             'Pcore':ind_obj.p_core,
             'Temp':ind_obj.t_winding
            } for ind_obj in inductor_set.values()
           ]    
    return pd.merge(df_ind_family,pd.DataFrame.from_dict(pset))